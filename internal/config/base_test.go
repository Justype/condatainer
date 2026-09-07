package config

import (
	"os"
	"path/filepath"
	"testing"

	"github.com/Justype/condatainer/catalog"
)

func TestEnsureDefaultDistroRecordsOnce(t *testing.T) {
	root := t.TempDir()
	if err := os.WriteFile(filepath.Join(root, "source.json"),
		[]byte(`{"schema":1,"default_distro":"ubuntu24"}`), 0o644); err != nil {
		t.Fatal(err)
	}
	if err := os.MkdirAll(filepath.Join(root, "recipes"), 0o755); err != nil {
		t.Fatal(err)
	}

	cfgDir := t.TempDir()
	t.Setenv("XDG_CONFIG_HOME", cfgDir)
	oldBase := Global.DefaultDistro
	Global.DefaultDistro = ""
	Global.Sources = []catalog.Spec{{Name: "cnt", Base: root}}
	ResetCatalog()
	t.Cleanup(func() { Global.DefaultDistro = oldBase; Global.Sources = nil; ResetCatalog() })

	cat, err := OpenCatalog(t.Context())
	if err != nil {
		t.Fatal(err)
	}

	if got := EnsureDefaultDistro(cat); got != "ubuntu24" {
		t.Fatalf("EnsureDefaultDistro = %q, want ubuntu24", got)
	}
	if Global.DefaultDistro != "ubuntu24" {
		t.Errorf("Global.DefaultDistro = %q", Global.DefaultDistro)
	}
	if BaseRecipeName() != "ubuntu24/base" {
		t.Errorf("BaseRecipeName = %q", BaseRecipeName())
	}

	// Sticky: a later upstream default does not revise what was recorded.
	Global.DefaultDistro = "ubuntu22"
	if got := EnsureDefaultDistro(cat); got != "ubuntu22" {
		t.Errorf("EnsureDefaultDistro overwrote a recorded distro: %q", got)
	}
	if def := SourceDefaultDistro(cat); def != "ubuntu24" {
		t.Errorf("SourceDefaultDistro = %q, want the source's newer recommendation", def)
	}
}
