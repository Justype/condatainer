package cmd

import (
	"os"
	"path/filepath"
	"testing"

	"github.com/Justype/condatainer/catalog"
	"github.com/Justype/condatainer/internal/artifact/meta"
)

func TestPresentListOverlayTrustsRecordedType(t *testing.T) {
	tests := []struct {
		name     string
		encoded  string
		runtime  meta.Runtime
		recorded bool
		want     listPresentation
	}{
		{
			name: "deep explicitly app", encoded: "team--tool--1.0", recorded: true,
			runtime: meta.Runtime{Name: "team/tool/1.0", Type: catalog.TypeApp},
			want:    listPresentation{name: "team/tool", version: "1.0"},
		},
		{
			name: "shallow explicitly data", encoded: "reference--1", recorded: true,
			runtime: meta.Runtime{Name: "reference/1", Type: catalog.TypeData},
			want:    listPresentation{name: "reference/1", data: true},
		},
		{
			name: "os", encoded: "ubuntu24--tool--1", recorded: true,
			runtime: meta.Runtime{Name: "ubuntu24/tool/1", Type: catalog.TypeOS},
			want:    listPresentation{name: "ubuntu24/tool/1", version: "(system app)"},
		},
		{
			name: "legacy depth fallback", encoded: "reference--genome--1", recorded: false,
			want: listPresentation{name: "reference/genome/1", data: true},
		},
	}
	for _, tc := range tests {
		t.Run(tc.name, func(t *testing.T) {
			if got := presentListOverlay(tc.encoded, tc.runtime, tc.recorded); got != tc.want {
				t.Fatalf("presentListOverlay() = %#v, want %#v", got, tc.want)
			}
		})
	}
}

// `list` is the by-name view of what is installed. Store entries are addressed
// by identity and reached through the `store` commands, so one name never
// appears here at several identities.
func TestScanOverlaysByDirDoesNotDescendTheStore(t *testing.T) {
	root := t.TempDir()
	for _, path := range []string{
		filepath.Join(root, "samtools--1.23.1.sqf"),
		filepath.Join(root, "store", "samtools--1.23.1@abcdef012345.sqf"),
	} {
		if err := os.MkdirAll(filepath.Dir(path), 0o775); err != nil {
			t.Fatal(err)
		}
		if err := os.WriteFile(path, nil, 0o664); err != nil {
			t.Fatal(err)
		}
	}

	dirs := scanOverlaysByDir([]string{root}, &SearchQuery{})
	if len(dirs) != 1 {
		t.Fatalf("dirs = %v, want one", dirs)
	}
	for name, path := range dirs[0].Paths {
		if filepath.Dir(path) != root {
			t.Errorf("%s listed from %s, want only the flat root", name, path)
		}
	}
	if len(dirs[0].Paths) > 1 {
		t.Errorf("paths = %v, want only the flat image", dirs[0].Paths)
	}
}
