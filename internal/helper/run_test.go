package helper

import (
	"context"
	"os"
	"path/filepath"
	"testing"

	"github.com/Justype/condatainer/catalog"
	"github.com/Justype/condatainer/internal/config"
)

// checkAndInstallNamedOverlays resolves a bare, distro-unprefixed name
// against what's already installed via catalog.SolveName, with no build and
// no "condatainer create" subprocess when a local overlay already satisfies
// it.
func TestCheckAndInstallNamedOverlaysReusesInstalled(t *testing.T) {
	recipesRoot := t.TempDir()
	write := func(rel, body string) {
		p := filepath.Join(recipesRoot, filepath.FromSlash(rel))
		if err := os.MkdirAll(filepath.Dir(p), 0o755); err != nil {
			t.Fatal(err)
		}
		if err := os.WriteFile(p, []byte(body), 0o644); err != nil {
			t.Fatal(err)
		}
	}
	write("recipes/ubuntu24/rstudio-server.def",
		"#DESC:RStudio Server {version}\n#TARGET:ubuntu24/rstudio-server/{version}\n#PH:version:2025.05.0-160,2026.09.0-174\n")

	imagesDir := t.TempDir()
	installedName := "ubuntu24--rstudio-server--2026.09.0-174.sqf"
	if err := os.WriteFile(filepath.Join(imagesDir, installedName), []byte{}, 0o644); err != nil {
		t.Fatal(err)
	}

	prevSources, prevDistro := config.Global.Sources, config.Global.DefaultDistro
	prevPaths := config.GlobalDataPaths
	config.Global.Sources = []catalog.Spec{{Name: "test", Base: recipesRoot}}
	config.Global.DefaultDistro = "ubuntu24"
	config.GlobalDataPaths.ImagesDirs = []string{imagesDir}
	config.ResetCatalog()
	t.Cleanup(func() {
		config.Global.Sources, config.Global.DefaultDistro = prevSources, prevDistro
		config.GlobalDataPaths = prevPaths
		config.ResetCatalog()
	})

	paths, err := checkAndInstallNamedOverlays(context.Background(), []string{"rstudio-server"})
	if err != nil {
		t.Fatalf("checkAndInstallNamedOverlays: %v", err)
	}
	want := filepath.Join(imagesDir, installedName)
	if len(paths) != 1 || paths[0] != want {
		t.Errorf("paths = %v, want [%s] — should reuse the installed overlay with no build", paths, want)
	}
}
