package server

import (
	"fmt"

	"github.com/Justype/condatainer/internal/config"
)

// helperBookmarksStore returns the jsonStore backing ~/.local/state/condatainer/helper-bookmarks.json.
func helperBookmarksStore() (*jsonStore[[]string], error) {
	path := config.GetHelperBookmarksPath()
	if path == "" {
		return nil, fmt.Errorf("cannot determine helper-bookmarks path")
	}
	return newJSONStore[[]string](path), nil
}

// listHelperBookmarks returns all bookmarked helper names. Returns an empty
// slice, not an error, if the helper-bookmarks file does not exist yet.
func listHelperBookmarks() ([]string, error) {
	store, err := helperBookmarksStore()
	if err != nil {
		return nil, err
	}
	list, err := store.Load()
	if err != nil {
		return nil, err
	}
	if list == nil {
		list = []string{}
	}
	return list, nil
}

// addHelperBookmark stars name (idempotent — starring an already-bookmarked
// name is a no-op). Returns the full updated list.
func addHelperBookmark(name string) ([]string, error) {
	store, err := helperBookmarksStore()
	if err != nil {
		return nil, err
	}
	return store.Update(func(list *[]string) error {
		if *list == nil {
			*list = []string{}
		}
		for _, n := range *list {
			if n == name {
				return nil
			}
		}
		*list = append(*list, name)
		return nil
	})
}

// removeHelperBookmark unstars name (a no-op if it was not bookmarked).
// Returns the full updated list.
func removeHelperBookmark(name string) ([]string, error) {
	store, err := helperBookmarksStore()
	if err != nil {
		return nil, err
	}
	return store.Update(func(list *[]string) error {
		if *list == nil {
			*list = []string{}
		}
		filtered := (*list)[:0]
		for _, n := range *list {
			if n != name {
				filtered = append(filtered, n)
			}
		}
		*list = filtered
		return nil
	})
}
