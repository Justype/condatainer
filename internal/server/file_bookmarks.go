package server

import (
	"fmt"
	"os"
	"path/filepath"
	"time"

	"github.com/Justype/condatainer/internal/config"
)

// fileBookmark is one starred filesystem path shown in the dashboard's file
// tree, persisted server-side (not in browser localStorage) so it's shared
// across browsers/devices for the same user.
type fileBookmark struct {
	Path    string    `json:"path"`
	Label   string    `json:"label,omitempty"`
	AddedAt time.Time `json:"added_at"`
}

// fileBookmarksStore returns the jsonStore backing ~/.local/state/condatainer/file-bookmarks.json.
func fileBookmarksStore() (*jsonStore[[]fileBookmark], error) {
	path := config.GetFileBookmarksPath()
	if path == "" {
		return nil, fmt.Errorf("cannot determine file-bookmarks path")
	}
	return newJSONStore[[]fileBookmark](path), nil
}

// listFileBookmarks returns all file bookmarks (oldest first). Returns an
// empty slice, not an error, if the file-bookmarks file does not exist yet.
func listFileBookmarks() ([]fileBookmark, error) {
	store, err := fileBookmarksStore()
	if err != nil {
		return nil, err
	}
	list, err := store.Load()
	if err != nil {
		return nil, err
	}
	if list == nil {
		list = []fileBookmark{}
	}
	return list, nil
}

// addFileBookmark stars path, assigning label (defaults to
// filepath.Base(path) if empty). Re-adding an already-bookmarked path
// updates its label rather than duplicating the entry. path must currently
// exist on disk. Returns the full updated list.
func addFileBookmark(path, label string) ([]fileBookmark, error) {
	store, err := fileBookmarksStore()
	if err != nil {
		return nil, err
	}
	clean := filepath.Clean(os.ExpandEnv(path))
	abs, err := filepath.Abs(clean)
	if err != nil {
		return nil, fmt.Errorf("invalid path: %w", err)
	}
	if _, err := os.Stat(abs); err != nil {
		return nil, fmt.Errorf("path does not exist: %w", err)
	}
	if label == "" {
		label = filepath.Base(abs)
	}

	return store.Update(func(list *[]fileBookmark) error {
		if *list == nil {
			*list = []fileBookmark{}
		}
		for i := range *list {
			if (*list)[i].Path == abs {
				(*list)[i].Label = label
				return nil
			}
		}
		*list = append(*list, fileBookmark{Path: abs, Label: label, AddedAt: time.Now()})
		return nil
	})
}

// removeFileBookmark unstars path (a no-op if it was not bookmarked).
// Returns the full updated list.
func removeFileBookmark(path string) ([]fileBookmark, error) {
	store, err := fileBookmarksStore()
	if err != nil {
		return nil, err
	}
	clean := filepath.Clean(os.ExpandEnv(path))
	abs, err := filepath.Abs(clean)
	if err != nil {
		return nil, fmt.Errorf("invalid path: %w", err)
	}

	return store.Update(func(list *[]fileBookmark) error {
		if *list == nil {
			*list = []fileBookmark{}
		}
		filtered := (*list)[:0]
		for _, b := range *list {
			if b.Path != abs {
				filtered = append(filtered, b)
			}
		}
		*list = filtered
		return nil
	})
}
