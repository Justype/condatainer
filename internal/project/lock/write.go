package lock

import (
	"fmt"
	"os"
	"path/filepath"
	"strings"
	"time"

	"github.com/Justype/condatainer/internal/utils"
)

// stagingPrefix marks a directory or file mid-write. A dotted name keeps it out
// of the artifact listing, so an interrupted write is invisible rather than a
// half-formed artifact.
const stagingPrefix = "."

// stagingName builds a unique sibling name. Uniqueness only has to survive
// concurrent writers to one project, which pid and nanosecond give.
func stagingName(kind string) string {
	return fmt.Sprintf("%s%s-%d-%d", stagingPrefix, kind, os.Getpid(), time.Now().UnixNano())
}

// Publish writes a lock atomically and then removes artifact directories no
// published selection reaches.
//
// The order is the transaction. Marshalling validates, a complete document is
// written to a temporary sibling and renamed over the old one, and only then is
// anything deleted. A failure before the rename leaves the previous lock
// authoritative and at worst leaves harmless unreferenced directories behind;
// pruning never runs against a lock that was not published.
func Publish(root string, l *Lock) error {
	data, err := l.Marshal()
	if err != nil {
		return err
	}
	dir := Dir(root)
	if err := utils.MkdirAllShared(dir); err != nil {
		return err
	}

	staged := filepath.Join(dir, stagingName("lock")+".json")
	defer os.Remove(staged) //nolint:errcheck

	temp, err := utils.CreateFileWritable(staged)
	if err != nil {
		return fmt.Errorf("cannot stage the lock: %w", err)
	}
	if _, err := temp.Write(data); err != nil {
		temp.Close() //nolint:errcheck
		return fmt.Errorf("cannot write the lock: %w", err)
	}
	if err := temp.Sync(); err != nil {
		temp.Close() //nolint:errcheck
		return fmt.Errorf("cannot flush the lock: %w", err)
	}
	if err := temp.Close(); err != nil {
		return err
	}
	if err := os.Rename(staged, FilePath(root)); err != nil {
		return fmt.Errorf("cannot publish the lock: %w", err)
	}
	utils.ShareWithParentGroup(FilePath(root))

	return Prune(root, l)
}

// Prune removes artifact directories the published lock does not reach.
//
// Reachability is recomputed from the published file, and it is the closure
// rather than the selection set: a dependency directory is reached through its
// parent's manifest edges, so pruning against selections alone would delete
// exactly the entries a rebuild needs. Anything that fails to verify is left
// alone — an entry that cannot be read cannot be shown to be unreachable, and
// deleting on a read error would turn a corrupt file into data loss.
func Prune(root string, l *Lock) error {
	present, err := listArtifactDirs(root)
	if err != nil {
		return err
	}
	if len(present) == 0 {
		return nil
	}
	verified, _ := Verify(root, l)

	for _, relative := range present {
		if verified.Reachable[relative] {
			continue
		}
		if _, readable := verified.Entries[relative]; !readable {
			continue
		}
		if err := os.RemoveAll(filepath.Join(Dir(root), relative)); err != nil {
			return fmt.Errorf("cannot prune %s: %w", relative, err)
		}
	}
	return nil
}

// StageArtifact copies a prepared artifact directory into cnt-lock/artifacts/
// under its final name, atomically.
//
// The staged directory is built as a temporary sibling on the same filesystem
// and renamed into place, so a reader never sees a half-written artifact. An
// entry that already exists is left as it is: one (name, identity) has one set
// of records, and rewriting them could only replace them with different bytes
// claiming the same key.
func StageArtifact(root, entryName string, files map[string][]byte) (string, error) {
	if entryName != filepath.Base(entryName) || entryName == "." || entryName == ".." {
		return "", fmt.Errorf("%w: artifact entry %q is not a plain name", ErrInvalid, entryName)
	}
	relative := ArtifactPath(entryName)
	if err := validArtifactPath(relative); err != nil {
		return "", fmt.Errorf("%w: %v", ErrInvalid, err)
	}
	final := filepath.Join(Dir(root), relative)
	if _, err := os.Stat(final); err == nil {
		return relative, nil
	}
	if err := utils.MkdirAllShared(ArtifactsPath(root)); err != nil {
		return "", err
	}

	staging := filepath.Join(ArtifactsPath(root), stagingName("staging"))
	if err := utils.MkdirAllShared(staging); err != nil {
		return "", fmt.Errorf("cannot stage %s: %w", entryName, err)
	}
	defer os.RemoveAll(staging) //nolint:errcheck

	for name, data := range files {
		if name != filepath.Base(name) || name == "." || name == ".." {
			return "", fmt.Errorf("%w: source file %q is not a plain name", ErrInvalid, name)
		}
		if len(data) > MaxSourceBytes {
			return "", fmt.Errorf("%w: source %s is %d bytes, over the %d limit", ErrInvalid, name, len(data), MaxSourceBytes)
		}
		file, err := utils.CreateFileWritable(filepath.Join(staging, name))
		if err != nil {
			return "", err
		}
		if _, err := file.Write(data); err != nil {
			file.Close() //nolint:errcheck
			return "", err
		}
		if err := file.Close(); err != nil {
			return "", err
		}
	}

	if err := os.Rename(staging, final); err != nil {
		// A concurrent writer producing the same (name, identity) wins; its
		// bytes are ours by construction.
		if _, statErr := os.Stat(final); statErr == nil {
			return relative, nil
		}
		return "", fmt.Errorf("cannot publish %s: %w", entryName, err)
	}
	utils.ShareWithParentGroup(final)
	return relative, nil
}

// SweepStaging removes staging directories abandoned by an interrupted write.
// They are never readable as artifacts — the name is dotted and Verify only
// reads what listArtifactDirs returns — so this is housekeeping, not repair.
func SweepStaging(root string) error {
	entries, err := os.ReadDir(ArtifactsPath(root))
	if os.IsNotExist(err) {
		return nil
	}
	if err != nil {
		return err
	}
	for _, entry := range entries {
		if !entry.IsDir() || !strings.HasPrefix(entry.Name(), stagingPrefix) {
			continue
		}
		if err := os.RemoveAll(filepath.Join(ArtifactsPath(root), entry.Name())); err != nil {
			return err
		}
	}
	return nil
}
