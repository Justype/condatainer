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
	prevPaths, prevBase := config.GlobalDataPaths, config.Global.Base
	config.GlobalDataPaths.ImagesDirs = []string{dir}
	config.Global.Base = "ubuntu24"
	t.Cleanup(func() { config.GlobalDataPaths, config.Global.Base = prevPaths, prevBase })
	return dir
}

// There is no overlay-only execution: without a root there is no container to
// start, so a missing base fails here rather than inside Apptainer.
func TestEnsureDefaultsRequiresAnInstalledBase(t *testing.T) {
	dir := withImageDir(t)

	if _, err := (Options{}).ensureDefaults(); err == nil {
		t.Fatal("execution was configured with no base image")
	}

	base := filepath.Join(dir, "ubuntu24--base.sif")
	if err := os.WriteFile(base, []byte("SIF"), 0o644); err != nil {
		t.Fatal(err)
	}
	opts, err := (Options{}).ensureDefaults()
	if err != nil {
		t.Fatalf("ensureDefaults: %v", err)
	}
	if opts.BaseImage != base {
		t.Errorf("BaseImage = %q, want %q", opts.BaseImage, base)
	}
	if len(opts.Command) == 0 {
		t.Error("no default command")
	}
}

// An explicit -b naming nothing is the user's mistake, and is reported as that
// rather than falling back to the configured base.
func TestEnsureDefaultsRejectsMissingExplicitBase(t *testing.T) {
	dir := withImageDir(t)
	if err := os.WriteFile(filepath.Join(dir, "ubuntu24--base.sif"), []byte("SIF"), 0o644); err != nil {
		t.Fatal(err)
	}

	missing := filepath.Join(dir, "nope.sif")
	_, err := Options{BaseImage: missing}.ensureDefaults()
	if err == nil {
		t.Fatal("a base image that does not exist was accepted")
	}
	if !strings.Contains(err.Error(), missing) {
		t.Errorf("err = %v, want it to name %s", err, missing)
	}
}
