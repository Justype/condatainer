package config

import (
	"os"
	"path/filepath"
	"strings"
	"testing"
)

// withImageDir points every image search at one scratch directory.
func withImageDir(t *testing.T) string {
	t.Helper()
	dir := filepath.Join(t.TempDir(), "images")
	if err := os.MkdirAll(dir, 0o755); err != nil {
		t.Fatal(err)
	}
	prev := GlobalDataPaths
	GlobalDataPaths.ImagesDirs = []string{dir}
	t.Cleanup(func() { GlobalDataPaths = prev })
	return dir
}

// withBase records a configured base for the duration of the test.
func withBase(t *testing.T, base string) {
	t.Helper()
	prev := Global.Base
	Global.Base = base
	t.Cleanup(func() { Global.Base = prev })
}

// GetBaseImage answers "which base is installed", not "where might one go".
// It used to fall back to the write path, handing every caller a path that did
// not exist and leaving Apptainer to report it.
func TestGetBaseImageRequiresAnInstalledFile(t *testing.T) {
	dir := withImageDir(t)
	withBase(t, "ubuntu24")

	if path, err := GetBaseImage(); err == nil {
		t.Fatalf("GetBaseImage = %q with nothing installed, want an error", path)
	} else if !strings.Contains(err.Error(), "ubuntu24/base") {
		t.Errorf("err = %v, want it to name the base it looked for", err)
	}

	installed := filepath.Join(dir, "ubuntu24--base.sqf")
	if err := os.WriteFile(installed, []byte("SQF"), 0o644); err != nil {
		t.Fatal(err)
	}
	got, err := GetBaseImage()
	if err != nil {
		t.Fatalf("GetBaseImage: %v", err)
	}
	if got != installed {
		t.Errorf("GetBaseImage = %q, want %q", got, installed)
	}
}

// With no base configured the failure is about configuration, not a lookup.
func TestGetBaseImageWithoutConfiguredBase(t *testing.T) {
	withImageDir(t)
	withBase(t, "")

	_, err := GetBaseImage()
	if err == nil {
		t.Fatal("an unconfigured base resolved to something")
	}
	if !strings.Contains(err.Error(), "no base configured") {
		t.Errorf("err = %v, want it to say no base is configured", err)
	}
}

// The name a build writes under (NewBuildObject's own "/" → "--" + ".sqf"
// transform of the resolved recipe name) and the name FindBaseImage looks for
// have to agree, or a rebuilt base is written where nothing will look for it.
func TestBaseWritePathMatchesFindPath(t *testing.T) {
	dir := withImageDir(t)
	withBase(t, "ubuntu24")

	fileName := BaseImageFileName()
	if fileName == "" {
		t.Fatal("BaseImageFileName is empty with a base configured")
	}
	if err := os.WriteFile(filepath.Join(dir, fileName), []byte("SQF"), 0o644); err != nil {
		t.Fatal(err)
	}
	if found := FindBaseImage(); found == "" {
		t.Errorf("a base written as %s is not found by FindBaseImage", fileName)
	}
}
