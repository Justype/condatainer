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

	// Apptainer stores exactly these bytes at /.singularity.d/Singularity when
	// it builds the same URI, so an image imported from such a build derives the
	// same recipe digest. Hashing is over the comment-stripped def, which is why
	// the headers are absent here and the trailing blank line is not.
	t.Run("docker URI matches what Apptainer synthesizes", func(t *testing.T) {
		path, err := synthesizeDefFromURI("docker://ubuntu:22.04", tmpDir)
		if err != nil {
			t.Fatalf("unexpected error: %v", err)
		}
		data, err := os.ReadFile(path)
		if err != nil {
			t.Fatalf("read def: %v", err)
		}
		const want = "bootstrap: docker\nfrom: ubuntu:22.04\n\n"
		if got := string(catalog.StripComments(data)); got != want {
			t.Errorf("stripped def = %q, want %q", got, want)
		}
		if !strings.Contains(string(data), "docker://ubuntu:22.04") {
			t.Errorf("header should record the source URI, got:\n%s", data)
		}
	})

	t.Run("rejects non-URI source", func(t *testing.T) {
		if _, err := synthesizeDefFromURI("/path/to/local.def", tmpDir); err == nil {
			t.Error("expected error for non-URI source, got nil")
		}
	})
}

// The staged metadata reaches the image by being written into the sandbox, so
// what the pack sees at <sandbox>/.cnt is what ends up at /.cnt.
func TestCopyMetaIntoSandbox(t *testing.T) {
	dir := t.TempDir()
	metaDir := filepath.Join(dir, "staged")
	sandbox := filepath.Join(dir, "rootfs")
	if err := os.MkdirAll(metaDir, 0o755); err != nil {
		t.Fatal(err)
	}
	if err := os.MkdirAll(sandbox, 0o755); err != nil {
		t.Fatal(err)
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

	if err := copyMetaIntoSandbox(metaDir, sandbox); err != nil {
		t.Fatalf("copyMetaIntoSandbox: %v", err)
	}
	if _, err := os.Stat(filepath.Join(sandbox, meta.DirName, meta.FileName)); err != nil {
		t.Errorf("manifest did not land in the sandbox: %v", err)
	}
}

// A definition with nothing staged still builds.
func TestCopyMetaIntoSandboxWithoutMetadata(t *testing.T) {
	sandbox := t.TempDir()
	if err := copyMetaIntoSandbox("", sandbox); err != nil {
		t.Fatalf("unexpected error: %v", err)
	}
	if _, err := os.Stat(filepath.Join(sandbox, meta.DirName)); err == nil {
		t.Error("created a .cnt directory for metadata that was never staged")
	}
}
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

	// The directives still have to be readable — Apptainer parses a def header
	// case-insensitively, and so does parseBootstrap, which resolves the digest.
	if boot := parseBootstrap(data); boot.Agent != "docker" || boot.From != "ubuntu:22.04" {
		t.Errorf("parseBootstrap = %+v, want {docker ubuntu:22.04}:\n%s", boot, data)
	}
}

// The script's markers are read per line: Apptainer writes its own greetings and
// warnings onto the same streams, so a marker never arrives alone.
func TestReadBaseTools(t *testing.T) {
	noisy := "INFO:    Converting SIF file to temporary sandbox...\n" +
		"WARNING: group: unknown groupid 10295\n" +
		"MISSING: mksquashfs micromamba\n"
	report := readBaseTools(noisy)
	if got := strings.Join(report.missing, ","); got != "mksquashfs,micromamba" {
		t.Errorf("missing = %q, want mksquashfs,micromamba", got)
	}
	if report.noApptainer {
		t.Error("apptainer reported absent when the script never said so")
	}

	report = readBaseTools("NOAPPTAINER\n")
	if len(report.missing) > 0 {
		t.Errorf("missing = %v, want none", report.missing)
	}
	if !report.noApptainer {
		t.Error("the optional-tool marker was not read")
	}

	if report := readBaseTools("WARNING: nothing to say\n"); len(report.missing) > 0 || report.noApptainer {
		t.Errorf("a complete base was read as incomplete: %+v", report)
	}
}

// Every requirement has to reach the script, or the check silently stops
// covering one.
func TestBaseToolsScriptNamesEveryTool(t *testing.T) {
	for _, name := range append(append([]string{}, baseTools...), baseToolPaths...) {
		if !strings.Contains(baseToolsScript, name) {
			t.Errorf("the check never asks for %s", name)
		}
	}
}
