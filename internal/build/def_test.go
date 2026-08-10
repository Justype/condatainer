package build

import (
	"bytes"
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/Justype/condatainer/catalog"
	"github.com/Justype/condatainer/internal/artifact/meta"
)

func TestSynthesizeDefFromURI(t *testing.T) {
	tmpDir := t.TempDir()

	t.Run("docker URI maps scheme and image", func(t *testing.T) {
		path, err := synthesizeDefFromURI("docker://ubuntu:22.04", tmpDir)
		if err != nil {
			t.Fatalf("unexpected error: %v", err)
		}
		data, err := os.ReadFile(path)
		if err != nil {
			t.Fatalf("read def: %v", err)
		}
		content := string(data)
		if !strings.Contains(content, "Bootstrap: docker") {
			t.Errorf("missing Bootstrap: docker, got:\n%s", content)
		}
		if !strings.Contains(content, "From: ubuntu:22.04") {
			t.Errorf("missing From: ubuntu:22.04, got:\n%s", content)
		}
		if !strings.Contains(content, "docker://ubuntu:22.04") {
			t.Errorf("header should record the source URI, got:\n%s", content)
		}
	})

	t.Run("rejects non-URI source", func(t *testing.T) {
		if _, err := synthesizeDefFromURI("/path/to/local.def", tmpDir); err == nil {
			t.Error("expected error for non-URI source, got nil")
		}
	})
}

func TestWriteRecordingDef(t *testing.T) {
	tmpDir := t.TempDir()
	cleanPath := filepath.Join(tmpDir, "clean.def")
	// No trailing newline, to exercise the newline-normalizing branch.
	if err := os.WriteFile(cleanPath, []byte("Bootstrap: docker\nFrom: alpine:3.19"), 0o644); err != nil {
		t.Fatalf("write clean def: %v", err)
	}

	metaDir := filepath.Join(tmpDir, meta.DirName)
	if err := meta.StageRuntime(metaDir, meta.Runtime{
		SchemaVersion: meta.SchemaVersion,
		Name:          "ubuntu24/base",
		Type:          catalog.TypeBase,
		Platform:      meta.NativePlatform(),
	}); err != nil {
		t.Fatalf("StageRuntime: %v", err)
	}
	if err := meta.StageManifest(metaDir, meta.Manifest{
		SchemaVersion: meta.SchemaVersion,
		Name:          "ubuntu24/base",
		Type:          catalog.TypeBase,
		BuildType:     "def",
		Platform:      meta.NativePlatform(),
	}); err != nil {
		t.Fatalf("StageManifest: %v", err)
	}

	recPath, err := writeRecordingDef(cleanPath, metaDir, tmpDir, "")
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}
	data, err := os.ReadFile(recPath)
	if err != nil {
		t.Fatalf("read recording def: %v", err)
	}
	content := string(data)

	// Original directives preserved.
	if !strings.Contains(content, "Bootstrap: docker") || !strings.Contains(content, "From: alpine:3.19") {
		t.Errorf("recording def dropped original directives:\n%s", content)
	}
	// Every staged document rides in on an appended %files section — a def build
	// has no packing step to add them during.
	if !strings.Contains(content, "%files") {
		t.Errorf("recording def has no %%files section:\n%s", content)
	}
	for src, dst := range map[string]string{meta.FileName: meta.Path, meta.RuntimeFileName: meta.RuntimePath} {
		want := "    " + filepath.Join(metaDir, src) + " " + dst
		if !strings.Contains(content, want) {
			t.Errorf("recording def missing embed line %q:\n%s", want, content)
		}
	}

	// The definition is not copied in: Apptainer already writes it to the rootfs
	// at /.singularity.d/Singularity, which survives extraction to .sqf.
	if strings.Contains(content, ".cnt-build-script") {
		t.Errorf("recording def embeds a redundant copy of the definition:\n%s", content)
	}
}

// A def build with nothing staged still has to produce a usable definition.
func TestWriteRecordingDefWithoutMetadata(t *testing.T) {
	tmpDir := t.TempDir()
	cleanPath := filepath.Join(tmpDir, "clean.def")
	if err := os.WriteFile(cleanPath, []byte("Bootstrap: docker\nFrom: alpine:3.19\n"), 0o644); err != nil {
		t.Fatal(err)
	}

	recPath, err := writeRecordingDef(cleanPath, "", tmpDir, "")
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}
	data, err := os.ReadFile(recPath)
	if err != nil {
		t.Fatal(err)
	}
	if strings.Contains(string(data), meta.Path) {
		t.Errorf("embedded a manifest that was never staged:\n%s", data)
	}
}

// A URI build has no recipe, so the synthesized def is its only chance to carry
// metadata: the headers have to survive the same parser a real recipe goes through.
func TestSynthesizedDefCarriesMetadata(t *testing.T) {
	tmpDir := t.TempDir()
	const uri = "docker://ubuntu:22.04"

	defPath, err := synthesizeDefFromURI(uri, tmpDir)
	if err != nil {
		t.Fatalf("synthesizeDefFromURI: %v", err)
	}
	data, err := os.ReadFile(defPath)
	if err != nil {
		t.Fatal(err)
	}

	recipe, err := catalog.ParseRecipe("myubuntu.def", bytes.NewReader(data))
	if err != nil {
		t.Fatalf("synthesized def does not parse as a recipe: %v\n%s", err, data)
	}
	if !strings.Contains(recipe.Description, uri) {
		t.Errorf("description = %q, want it to name %s", recipe.Description, uri)
	}
	if recipe.URL != uri {
		t.Errorf("url = %q, want %q", recipe.URL, uri)
	}
	if recipe.Type != catalog.TypeOS {
		t.Errorf("type = %q, want %q — a bare image is a root layer, not an app", recipe.Type, catalog.TypeOS)
	}

	// Apptainer still has to be able to build it.
	if !strings.Contains(string(data), "Bootstrap: docker") || !strings.Contains(string(data), "From: ubuntu:22.04") {
		t.Errorf("synthesized def lost its directives:\n%s", data)
	}
}
