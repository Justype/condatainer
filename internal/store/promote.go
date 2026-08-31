package store

import (
	"errors"
	"fmt"
	"os"
	"path/filepath"

	"github.com/Justype/condatainer/internal/artifact/compare"
	"github.com/Justype/condatainer/internal/artifactcache"
	"github.com/Justype/condatainer/internal/config"
	"github.com/Justype/condatainer/internal/image"
	"github.com/Justype/condatainer/internal/image/producer"
	"github.com/Justype/condatainer/internal/utils"
)

// ErrShadowed reports a promotion that would be overruled by a nearer root.
var ErrShadowed = errors.New("promotion would be shadowed by a nearer root")

// Promotion is what a promotion did, or would do.
type Promotion struct {
	Name string `json:"name"`
	Root string `json:"root"`
	// Promoted is the artifact that now answers to the bare name.
	Promoted Candidate `json:"promoted"`
	// Demoted is the artifact that answered to it before, now in the same
	// root's store/. Zero when the root held no flat copy.
	Demoted Candidate `json:"demoted,omitzero"`
	// DemotedTo is where it went, empty when nothing was demoted.
	DemotedTo string `json:"demoted_to,omitempty"`
	// Shadows is a flat copy in a farther root that this promotion now
	// overrules. It is untouched and still resolvable by identity.
	Shadows string `json:"shadows,omitempty"`
	// AlreadyFlat reports that the requested identity already answered to the
	// name and nothing was done.
	AlreadyFlat bool `json:"already_flat,omitempty"`
}

// Promote makes one identity the artifact a bare name resolves to.
//
// It only ever renames, and only inside the candidate's own root. That is what
// makes it safe for everyone else reading that root: nothing leaves the root, so
// no identity disappears from anyone's readable set whatever their configuration,
// and the demoted artifact is still found by identity in store/. A move across
// roots would take it out of the readable set of every user whose configuration
// does not include the destination.
//
// A candidate in a root that a nearer root shadows is refused rather than
// renamed into place, because promoting it there would change nothing that reads
// resolve.
func Promote(name string, q IdentityQuery, dirs []string) (Promotion, error) {
	return promote(name, q, dirs, compare.Read)
}

func promote(name string, q IdentityQuery, dirs []string, read artifactReader) (Promotion, error) {
	if dirs == nil {
		dirs = config.GetImageSearchPaths()
	}
	report := scan(ScanOptions{Dirs: dirs, Name: name}, read, nil)
	candidate, err := resolveReport(name, q, report)
	if err != nil {
		return Promotion{}, err
	}
	result := Promotion{Name: candidate.Name, Root: candidate.Root, Promoted: candidate}
	if candidate.Layout == LayoutFlat {
		result.AlreadyFlat = true
		return result, nil
	}

	// Reads are nearest-first, so only the nearest flat copy decides the name.
	incumbent, found := nearestFlat(report, dirs)
	if found {
		switch {
		case incumbent.Root == candidate.Root:
			result.Demoted = incumbent
		case nearer(dirs, incumbent.Root, candidate.Root):
			return Promotion{}, fmt.Errorf("%w: %s already answers to %s; copy this identity into that root first",
				ErrShadowed, incumbent.Path, candidate.Name)
		default:
			// The incumbent is farther out. Promoting here wins outright and
			// leaves it in place, still reachable by identity.
			result.Shadows = incumbent.Path
		}
	}

	flatPath := filepath.Join(candidate.Root, image.EncodeArtifactName(candidate.Name)+".sqf")
	if err := promoteLocked(&result, candidate, flatPath, read); err != nil {
		return Promotion{}, err
	}
	return result, nil
}

// promoteLocked performs the renames while holding the bare name's producer
// lock, so a build cannot claim the name in the window where it is unclaimed.
//
// The demotion happens first and completely: a crash between the two renames
// leaves both artifacts in store/ with the bare name free, which is a legal
// state that re-running repairs. Staging the incumbent as a .part instead would
// leave a file GC's staging sweep deletes, turning a crash into a lost artifact.
func promoteLocked(result *Promotion, candidate Candidate, flatPath string, read artifactReader) error {
	guard, err := producer.AcquireLocal(flatPath)
	if err != nil {
		return err
	}
	defer guard.Release() //nolint:errcheck

	if result.Demoted.Path != "" {
		demotedTo, err := demote(result.Demoted, read)
		if err != nil {
			return err
		}
		result.DemotedTo = demotedTo
	}
	if err := renameChecked(candidate.Path, flatPath); err != nil {
		return err
	}
	utils.ShareWithParentGroup(flatPath)
	artifactcache.Default().Forget(candidate.Path)

	result.Promoted.Path, result.Promoted.Layout = flatPath, LayoutFlat
	return nil
}

// demote moves the incumbent flat artifact into its own root's store/, and
// reports where it landed. Its identity is regenerated rather than trusted, so
// it is filed under the address it actually has.
func demote(incumbent Candidate, read artifactReader) (string, error) {
	if err := checkFree(incumbent.Path); err != nil {
		return "", err
	}
	artifact, err := read(incumbent.Path)
	if err != nil {
		return "", fmt.Errorf("cannot demote %s: %w", incumbent.Path, err)
	}
	identity := artifact.IdentityRef()
	storeDir := filepath.Join(incumbent.Root, DirName)
	if err := utils.MkdirAllShared(storeDir); err != nil {
		return "", fmt.Errorf("cannot create store directory: %w", err)
	}
	for chars := DefaultPrefixChars; chars <= len(identity.SHA256); chars++ {
		filename, err := Filename(artifact.Name, identity, chars)
		if err != nil {
			return "", err
		}
		target := filepath.Join(storeDir, filename)
		if _, err := os.Lstat(target); err == nil {
			continue
		} else if !os.IsNotExist(err) {
			return "", err
		}
		if err := renameChecked(incumbent.Path, target); err != nil {
			return "", err
		}
		artifactcache.Default().Forget(incumbent.Path)
		return target, nil
	}
	return "", fmt.Errorf("%w: no filename remains for %s", ErrIdentityCollision, artifact.Name)
}

// checkFree refuses an artifact that is protected or being read. The two stay
// distinct: reporting a pinned image as "in use" sends the reader hunting for a
// container that is not running.
func checkFree(path string) error {
	lock, err := image.AcquireLock(path, true)
	if err != nil {
		return err
	}
	return lock.Close()
}

// renameChecked refuses to overwrite anything: every target here is either free
// or was just vacated under the producer lock.
func renameChecked(from, to string) error {
	if _, err := os.Lstat(to); err == nil {
		return fmt.Errorf("cannot move %s to %s: it already exists", filepath.Base(from), to)
	} else if !os.IsNotExist(err) {
		return err
	}
	if err := os.Rename(from, to); err != nil {
		return fmt.Errorf("failed to move %s to %s: %w", filepath.Base(from), to, err)
	}
	return nil
}

// nearestFlat returns the flat copy that currently answers to the name.
func nearestFlat(report Report, dirs []string) (Candidate, bool) {
	best, found := Candidate{}, false
	for _, candidate := range report.Candidates {
		if candidate.Layout != LayoutFlat {
			continue
		}
		if !found || nearer(dirs, candidate.Root, best.Root) {
			best, found = candidate, true
		}
	}
	return best, found
}

// nearer reports whether root a is read before root b.
func nearer(dirs []string, a, b string) bool {
	return rootIndex(dirs, a) < rootIndex(dirs, b)
}

func rootIndex(dirs []string, root string) int {
	root = filepath.Clean(root)
	for i, dir := range dirs {
		if filepath.Clean(dir) == root {
			return i
		}
	}
	return len(dirs)
}
