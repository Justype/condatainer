package config

import (
	"os"
	"path/filepath"
	"testing"

	"github.com/Justype/condatainer/catalog"
)

func TestEnsureBaseRecordsOnce(t *testing.T) {
	root := t.TempDir()
	if err := os.WriteFile(filepath.Join(root, "source.json"),
		[]byte(`{"schema":1,"default_base":"ubuntu24"}`), 0o644); err != nil {
		t.Fatal(err)
	}
	if err := os.MkdirAll(filepath.Join(root, "recipes"), 0o755); err != nil {
		t.Fatal(err)
	}

	cfgDir := t.TempDir()
	t.Setenv("XDG_CONFIG_HOME", cfgDir)
	oldBase := Global.Base
	Global.Base = ""
	Global.Sources = []catalog.Spec{{Name: "cnt", Base: root}}
	ResetCatalog()
	t.Cleanup(func() { Global.Base = oldBase; Global.Sources = nil; ResetCatalog() })

	cat, err := OpenCatalog(t.Context())
	if err != nil {
		t.Fatal(err)
	}

	if got := EnsureBase(cat); got != "ubuntu24" {
		t.Fatalf("EnsureBase = %q, want ubuntu24", got)
	}
	if Global.Base != "ubuntu24" {
		t.Errorf("Global.Base = %q", Global.Base)
	}
	if BaseRecipeName() != "ubuntu24/base" {
		t.Errorf("BaseRecipeName = %q", BaseRecipeName())
	}

	// Sticky: a later upstream default does not revise what was recorded.
	Global.Base = "ubuntu22"
	if got := EnsureBase(cat); got != "ubuntu22" {
		t.Errorf("EnsureBase overwrote a recorded base: %q", got)
	}
	if def := SourceDefaultBase(cat); def != "ubuntu24" {
		t.Errorf("SourceDefaultBase = %q, want the source's newer recommendation", def)
	}
}
