package container

import (
	"path/filepath"
	"testing"

	"github.com/Justype/condatainer/internal/config"
)

// withImageDirs points the image search paths at dirs for one test.
func withImageDirs(t *testing.T, dirs ...string) {
	t.Helper()
	prev := config.GlobalDataPaths.ImagesDirs
	config.GlobalDataPaths.ImagesDirs = dirs
	t.Cleanup(func() { config.GlobalDataPaths.ImagesDirs = prev })
}

// A filename is an address only inside an image directory, so that is the only
// place a recorded name can disagree with one.
func TestNameDiagnostic(t *testing.T) {
	imageDir := t.TempDir()
	elsewhere := t.TempDir()
	withImageDirs(t, imageDir)

	if d := nameDiagnostic(filepath.Join(imageDir, "samtools--1.22.sqf"), "samtools/1.22"); d != nil {
		t.Errorf("agreeing name warned: %s", d.Message)
	}
	d := nameDiagnostic(filepath.Join(imageDir, "samtools--1.22.sqf"), "r/4.4.3")
	if d == nil {
		t.Fatal("a recorded name that the filename does not address went unreported")
	}
	if d.Level != "warn" {
		t.Errorf("level = %q, want warn", d.Level)
	}
	if d := nameDiagnostic(filepath.Join(elsewhere, "scratch.sqf"), "r/4.4.3"); d != nil {
		t.Errorf("a path outside the image directories addresses no name: %s", d.Message)
	}
	if d := nameDiagnostic(filepath.Join(imageDir, "anything.sqf"), ""); d != nil {
		t.Errorf("an image recording no name has nothing to disagree with: %s", d.Message)
	}
}
