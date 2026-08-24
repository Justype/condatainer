package cmd

import (
	"strings"
	"testing"

	"github.com/Justype/condatainer/internal/config"
)

// No scope means every writable root, which store.GC allows for a report and
// refuses for an apply.
func TestStoreGCDirsIsEmptyWithoutAScope(t *testing.T) {
	dirs, err := storeGCDirs("", "")
	if err != nil {
		t.Fatal(err)
	}
	if dirs != nil {
		t.Fatalf("dirs = %v, want nil so GC walks every writable root", dirs)
	}
}

// A layer nobody can spell is an error, not an empty selection.
func TestStoreGCDirsRejectsAnUnknownLayer(t *testing.T) {
	_, err := storeGCDirs("", "nowhere")
	if err == nil {
		t.Fatal("an unknown layer was accepted")
	}
}

// "Nothing to collect" and "you named a root that does not exist" are answers
// someone acts on differently, so a filter matching nothing is an error.
func TestStoreGCDirsRejectsADirMatchingNothing(t *testing.T) {
	_, err := storeGCDirs("definitely-not-an-image-root-xyzzy", "")
	if err == nil {
		t.Fatal("a --dir matching no image directory was accepted")
	}
	if !strings.Contains(err.Error(), "xyzzy") {
		t.Errorf("error does not name the filter: %v", err)
	}
}

// --dir narrows whatever --layer selected, rather than replacing it.
func TestStoreGCDirsAppliesLayerThenDir(t *testing.T) {
	dirs, err := storeGCDirs("", "u")
	if err != nil {
		t.Skipf("no user-layer image directory on this host: %v", err)
	}
	for _, dir := range dirs {
		if config.ClassifyDataDir(dir) != config.LayerUser {
			t.Errorf("%s is not in the user layer", dir)
		}
	}
	if _, err := storeGCDirs("definitely-not-here-xyzzy", "u"); err == nil {
		t.Error("a --dir matching nothing inside the layer was accepted")
	}
}
