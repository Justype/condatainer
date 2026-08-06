package catalog

import (
	"context"
	"io/fs"
	"os"
	"path/filepath"
	"strings"
)

// dirBackend reads a filesystem source directly. There is no index to keep
// honest and nothing to cache: names map to paths structurally, so resolving
// one is a stat, and only enumeration walks the tree.
type dirBackend struct{ src *Source }

func (d *dirBackend) read(_ context.Context, path string) ([]byte, error) {
	return os.ReadFile(filepath.Join(d.src.Base, filepath.FromSlash(path)))
}

// entries walks recipes/ and parses each file's headers.
func (d *dirBackend) entries(_ context.Context) (map[string]*Entry, error) {
	root := filepath.Join(d.src.Base, recipesDir)
	out := map[string]*Entry{}

	err := filepath.WalkDir(root, func(p string, entry fs.DirEntry, err error) error {
		if err != nil || entry.IsDir() || skipName(entry.Name()) {
			return err
		}
		rel, err := filepath.Rel(root, p)
		if err != nil {
			return err
		}
		f, err := os.Open(p)
		if err != nil {
			return err
		}
		defer f.Close()

		rec, err := ParseRecipe(recipesDir+"/"+filepath.ToSlash(rel), f)
		if err != nil {
			return err
		}
		out[rec.Name] = &rec.Entry
		return nil
	})
	if err != nil {
		return nil, err
	}
	return out, nil
}

// skipName drops files a collection carries but does not offer: documentation,
// and whatever an editor or the OS left behind.
func skipName(name string) bool {
	switch strings.ToLower(name) {
	case "readme.md", ".gitignore", ".gitkeep", ".ds_store":
		return true
	}
	return strings.HasPrefix(name, ".")
}
