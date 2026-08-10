package capsule

import (
	"fmt"
	"os"
	"path/filepath"
	"sort"
	"strings"

	"github.com/Justype/condatainer/internal/artifact/meta"
)

// Entry is one artifact in an embedded closure.
type Entry struct {
	// Dir is the entry directory name, <name>@<short identity>.
	Dir string
	// Name is the artifact's name, with its slashes restored.
	Name string
	// Identity is the truncated identity from the directory name. It addresses
	// the entry; the full digest is inside its identity.record.
	Identity string
	// Files are the entry's file names, sorted.
	Files []string
}

// Entries lists a capsule directory. A missing capsule is not an error: only
// data has dependencies, so most artifacts have none.
func Entries(dir string) ([]Entry, error) {
	found, err := os.ReadDir(dir)
	if os.IsNotExist(err) {
		return nil, nil
	}
	if err != nil {
		return nil, fmt.Errorf("cannot list capsule %s: %w", dir, err)
	}

	var out []Entry
	for _, node := range found {
		if !node.IsDir() {
			return nil, fmt.Errorf("%w: %s is not an entry directory", ErrInvalid, node.Name())
		}
		entry, err := readEntry(dir, node.Name())
		if err != nil {
			return nil, err
		}
		out = append(out, entry)
	}
	sort.Slice(out, func(i, j int) bool { return out[i].Dir < out[j].Dir })
	return out, nil
}

// readEntry parses one entry directory name and lists what it holds.
func readEntry(dir, name string) (Entry, error) {
	artifact, identity, ok := strings.Cut(name, "@")
	if !ok || artifact == "" || identity == "" {
		return Entry{}, fmt.Errorf("%w: %q is not <name>@<identity>", ErrInvalid, name)
	}
	entry := Entry{
		Dir:      name,
		Name:     strings.ReplaceAll(artifact, "--", "/"),
		Identity: identity,
	}

	files, err := os.ReadDir(filepath.Join(dir, name))
	if err != nil {
		return Entry{}, fmt.Errorf("cannot list capsule entry %s: %w", name, err)
	}
	for _, file := range files {
		if file.IsDir() {
			return Entry{}, fmt.Errorf("%w: %s nests a directory; entries are flat", ErrInvalid, name)
		}
		entry.Files = append(entry.Files, file.Name())
	}
	sort.Strings(entry.Files)
	return entry, nil
}

// Validate reports whether a capsule can be trusted as read, for an artifact of
// the given name and identity.
//
// Two failures matter. A **self-reference** means the closure claims the artifact
// was built from itself; a dependency image exists before the artifact that
// mounts it, so this cannot arise from a real build. A **truncated entry** — one
// with no manifest.json — cannot be rebuilt from, so it is worse than an absent
// one, which at least says so.
func Validate(dir, name, identity string) error {
	entries, err := Entries(dir)
	if err != nil {
		return err
	}
	self := EntryName(name, identity)
	for _, entry := range entries {
		if entry.Dir == self {
			return fmt.Errorf("%w: %s contains itself", ErrInvalid, name)
		}
		if !hasFile(entry.Files, meta.FileName) {
			return fmt.Errorf("%w: %s carries no %s", ErrInvalid, entry.Dir, meta.FileName)
		}
	}
	return nil
}

func hasFile(files []string, want string) bool {
	for _, f := range files {
		if f == want {
			return true
		}
	}
	return false
}
