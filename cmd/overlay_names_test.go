package cmd

import (
	"context"
	"os"
	"path/filepath"
	"testing"

	"github.com/Justype/condatainer/catalog"
	"github.com/Justype/condatainer/internal/config"
)

func TestResolveOverlayValuesOutsideAProject(t *testing.T) {
	dir := t.TempDir()
	for _, f := range []string{"samtools--1.21.sqf", "samtools--1.22.1.sqf", "ubuntu24--xfce4.sqf"} {
		if err := os.WriteFile(filepath.Join(dir, f), nil, 0o644); err != nil {
			t.Fatal(err)
		}
	}
	prevSources, prevDistro, prevPaths := config.Global.Sources, config.Global.DefaultDistro, config.GlobalDataPaths
	config.Global.Sources, config.Global.DefaultDistro = []catalog.Spec(nil), "ubuntu24"
	config.GlobalDataPaths.ImagesDirs = []string{dir}
	config.ResetCatalog()
	t.Cleanup(func() {
		config.Global.Sources, config.Global.DefaultDistro, config.GlobalDataPaths = prevSources, prevDistro, prevPaths
		config.ResetCatalog()
	})
	t.Chdir(t.TempDir())

	got, err := resolveOverlayValues(context.Background(),
		[]string{"samtools", "samtools/1.22:ro", "xfce4", "nothing/1.0"}, nil, false)
	if err != nil {
		t.Fatal(err)
	}
	want := []string{
		filepath.Join(dir, "samtools--1.22.1.sqf"),
		filepath.Join(dir, "samtools--1.22.1.sqf") + ":ro",
		filepath.Join(dir, "ubuntu24--xfce4.sqf"),
		"nothing/1.0",
	}
	for i := range want {
		if got[i] != want[i] {
			t.Errorf("overlay %d = %q, want %q", i, got[i], want[i])
		}
	}
}

func TestResolveOverlayValuesConstraint(t *testing.T) {
	if _, err := resolveOverlayValues(context.Background(), []string{"star>=2.7.0"}, nil, false); err == nil {
		t.Error("a version constraint was accepted")
	}
	got, err := resolveOverlayValues(context.Background(), []string{"star>=2.7.0"}, nil, true)
	if err != nil || got[0] != "star>=2.7.0" {
		t.Errorf("got %v, %v; want the value left as it was", got, err)
	}
}
