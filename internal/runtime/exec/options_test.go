package exec

import (
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/Justype/condatainer/internal/config"
)

// withImageDir points every image search at one scratch directory and records a
// configured base, so a test decides for itself whether one is installed.
func withImageDir(t *testing.T) string {
	t.Helper()
	dir := filepath.Join(t.TempDir(), "images")
	if err := os.MkdirAll(dir, 0o755); err != nil {
		t.Fatal(err)
	}
	prevPaths, prevBase := config.GlobalDataPaths, config.Global.DefaultDistro
	config.GlobalDataPaths.ImagesDirs = []string{dir}
	config.Global.DefaultDistro = "ubuntu24"
	t.Cleanup(func() { config.GlobalDataPaths, config.Global.DefaultDistro = prevPaths, prevBase })
	return dir
}

// There is no overlay-only execution: without a root there is no container to
// start, so a missing base fails here rather than inside Apptainer.
func TestResolveBaseImageRequiresAnInstalledBase(t *testing.T) {
	dir := withImageDir(t)

	if _, err := (Options{}).resolveBaseImage(""); err == nil {
		t.Fatal("execution was configured with no base image")
	}

	base := filepath.Join(dir, "ubuntu24--base.sqf")
	if err := os.WriteFile(base, []byte("SQF"), 0o644); err != nil {
		t.Fatal(err)
	}
	opts, err := (Options{}).resolveBaseImage("")
	if err != nil {
		t.Fatalf("resolveBaseImage: %v", err)
	}
	if opts.BaseImage != base {
		t.Errorf("BaseImage = %q, want %q", opts.BaseImage, base)
	}
}

// A root pulled out of the requested overlays wins over the configured
// default, without requiring one to be installed at all.
func TestResolveBaseImageRootWinsOverConfiguredDefault(t *testing.T) {
	dir := withImageDir(t) // no base installed — the configured default is not needed

	root := filepath.Join(dir, "myos.sqf")
	if err := os.WriteFile(root, []byte("SQF"), 0o644); err != nil {
		t.Fatal(err)
	}
	opts, err := (Options{}).resolveBaseImage(root)
	if err != nil {
		t.Fatalf("resolveBaseImage: %v", err)
	}
	if opts.BaseImage != root {
		t.Errorf("BaseImage = %q, want %q", opts.BaseImage, root)
	}
}

// A base named directly (by the caller, or the configured default) that does
// not exist is the mistake it looks like, not silently accepted.
func TestResolveBaseImageRejectsAMissingBase(t *testing.T) {
	dir := withImageDir(t)
	if err := os.WriteFile(filepath.Join(dir, "ubuntu24--base.sqf"), []byte("SQF"), 0o644); err != nil {
		t.Fatal(err)
	}

	missing := filepath.Join(dir, "nope.sqf")
	_, err := Options{BaseImage: missing}.resolveBaseImage("")
	if err == nil {
		t.Fatal("a base image that does not exist was accepted")
	}
	if !strings.Contains(err.Error(), missing) {
		t.Errorf("err = %v, want it to name %s", err, missing)
	}
}

// A .sif root is checked for /bin/bash; unlike the .sqf cases above, this
// exercises a real bash lookup, so a file that isn't a real SIF at all still
// has to fail resolveBaseImage rather than being waved through.
func TestResolveBaseImageChecksBashInASif(t *testing.T) {
	dir := withImageDir(t)
	root := filepath.Join(dir, "myos.sif")
	if err := os.WriteFile(root, []byte("not a real sif"), 0o644); err != nil {
		t.Fatal(err)
	}

	if _, err := (Options{}).resolveBaseImage(root); err == nil {
		t.Fatal("resolveBaseImage: want an error for a .sif with no readable bash")
	}
}

// ensureDefaults no longer touches BaseImage at all — that is
// resolveBaseImage's job, run separately after overlay setup.
func TestEnsureDefaultsFillsCommandOnly(t *testing.T) {
	opts, err := (Options{}).ensureDefaults()
	if err != nil {
		t.Fatalf("ensureDefaults: %v", err)
	}
	if opts.BaseImage != "" {
		t.Errorf("BaseImage = %q, want ensureDefaults to leave it untouched", opts.BaseImage)
	}
	if len(opts.Command) == 0 {
		t.Error("no default command")
	}
}
