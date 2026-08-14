// Package capsule composes and reads /.cnt/provenance, the embedded closure of
// everything an artifact was built from.
//
// It is built by union, never by re-deriving anything: an artifact's capsule is
// its dependencies' records plus its dependencies' capsules, copied across
// unchanged. Two properties follow, and they are why this shape was chosen —
// completeness is inductive, so if every dependency's capsule is complete the
// union is complete, with no recursive resolution, no catalog access and no
// network at build time; and cycles cannot occur, because a dependency image
// existed before the artifact that mounts it.
package capsule

import (
	"errors"
	"fmt"
	"os"
	"path/filepath"
	"strings"

	"github.com/Justype/condatainer/internal/artifact/key"
	"github.com/Justype/condatainer/internal/artifact/meta"
	"github.com/Justype/condatainer/internal/image"
	"github.com/Justype/condatainer/internal/utils"
)

// DirName is the capsule directory inside /.cnt.
const DirName = "provenance"

// Path is where the capsule lives inside an image.
const Path = "/" + meta.DirName + "/" + DirName

// IdentityChars is how much of an identity an entry directory carries. It is
// addressing only: a reader regenerates the full digest from the entry sources
// rather than trusting twelve characters.
const IdentityChars = 12

// ErrInvalid reports a capsule that cannot be trusted as read.
var ErrInvalid = errors.New("invalid provenance capsule")

// Dep is one direct dependency to compose from: what it is called, what it says
// its identity is, and the image to read it out of.
type Dep struct {
	Name string
	// Identity is the complete key, empty when the image carries no records.
	// Composition compares both halves.
	Identity  meta.KeyRef
	ImagePath string
}

// EntryName is the directory one artifact's records live in: its name with
// slashes joined by `--`, exactly as image filenames and store entries already do
// so a dependency's own name never becomes directory levels, then `@` and the
// truncated identity. `@` cannot occur in a Conda name or version, so it never
// collides with the separator.
//
// Name and identity together address a record set; identity alone never does,
// because one solve published under two names has one identity.
func EntryName(name, identity string) string {
	short := strings.TrimPrefix(identity, "sha256:")
	if len(short) > IdentityChars {
		short = short[:IdentityChars]
	}
	return image.EncodeArtifactName(name) + "@" + short
}

// Compose stages the capsule for an artifact into metaDir, and reports whether
// the closure it produced is complete.
//
// Completeness is inherited, not recomputed: an artifact is complete when every
// dependency carried records *and* every dependency's own manifest said it was
// complete. One unrecorded image anywhere below makes everything above it
// incomplete, which is the honest answer.
func Compose(metaDir string, deps []Dep) (complete bool, err error) {
	complete = true
	if len(deps) == 0 {
		return complete, nil
	}

	dest := filepath.Join(metaDir, DirName)
	for _, dep := range deps {
		if dep.Identity.Empty() {
			// Nothing to copy and nothing to name it: an unrecorded dependency
			// is recorded in the manifest and the records, not here.
			complete = false
			continue
		}
		whole, err := composeOne(dest, dep)
		if err != nil {
			return false, err
		}
		if !whole {
			complete = false
		}
	}
	return complete, nil
}

// composeOne copies one dependency's manifest and rebuild sources into the
// capsule, then its own capsule entries across unchanged.
func composeOne(dest string, dep Dep) (complete bool, err error) {
	staging, err := os.MkdirTemp("", "cnt-capsule-")
	if err != nil {
		return false, fmt.Errorf("failed to create capsule staging dir: %w", err)
	}
	defer os.RemoveAll(staging) //nolint:errcheck

	// One extraction, not one read per file: every archive read spawns a process.
	if err := image.ExtractDir(dep.ImagePath, "/"+meta.DirName, staging); err != nil {
		return false, fmt.Errorf("cannot read provenance from %s: %w", dep.Name, err)
	}
	source := filepath.Join(staging, meta.DirName)
	manifest, err := readManifestDir(source)
	if err != nil {
		return false, fmt.Errorf("cannot read provenance manifest from %s: %w", dep.Name, err)
	}
	derived, err := key.VerifyDir(source, manifest)
	if err != nil {
		return false, fmt.Errorf("cannot verify provenance from %s: %w", dep.Name, err)
	}
	if manifest.Name != dep.Name || derived.Identity.Ref != dep.Identity {
		return false, fmt.Errorf("provenance from %s does not match the selected identity", dep.Name)
	}
	entry := filepath.Join(dest, EntryName(dep.Name, dep.Identity.Digest()))
	if err := utils.MkdirAllShared(entry); err != nil {
		return false, err
	}

	files := append([]string{meta.FileName}, manifest.Source.Files...)
	for _, name := range files {
		if err := copyFile(filepath.Join(source, name), filepath.Join(entry, name)); err != nil {
			return false, err
		}
	}

	// The dependency's own closure, copied across rather than re-derived.
	inherited, err := os.ReadDir(filepath.Join(source, DirName))
	if err == nil {
		for _, sub := range inherited {
			if !sub.IsDir() {
				continue
			}
			if err := copyTree(filepath.Join(source, DirName, sub.Name()), filepath.Join(dest, sub.Name())); err != nil {
				return false, err
			}
		}
	}

	return manifest.ProvenanceComplete == nil || *manifest.ProvenanceComplete, nil
}

// copyFile writes src to dst, creating dst's parent.
func copyFile(src, dst string) error {
	data, err := os.ReadFile(src)
	if err != nil {
		return fmt.Errorf("failed to read %s: %w", src, err)
	}
	if err := utils.MkdirAllShared(filepath.Dir(dst)); err != nil {
		return err
	}
	f, err := utils.CreateFileWritable(dst)
	if err != nil {
		return fmt.Errorf("failed to create %s: %w", dst, err)
	}
	if _, err := f.Write(data); err != nil {
		f.Close() //nolint:errcheck
		return fmt.Errorf("failed to write %s: %w", dst, err)
	}
	if err := f.Close(); err != nil {
		return fmt.Errorf("failed to close %s: %w", dst, err)
	}
	utils.ShareWithParentGroup(dst)
	return nil
}

// copyTree copies one capsule entry across. Deduplication is by directory name,
// which is (name, identity): a diamond stores the shared dependency once.
func copyTree(src, dst string) error {
	if _, err := os.Stat(dst); err == nil {
		return nil
	}
	entries, err := os.ReadDir(src)
	if err != nil {
		return fmt.Errorf("failed to list %s: %w", src, err)
	}
	if err := utils.MkdirAllShared(dst); err != nil {
		return err
	}
	for _, file := range entries {
		if file.IsDir() {
			continue // a capsule entry is flat
		}
		if err := copyFile(filepath.Join(src, file.Name()), filepath.Join(dst, file.Name())); err != nil {
			return err
		}
	}
	return nil
}
