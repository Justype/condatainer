package build

import (
	"context"
	"encoding/json"
	"os"
	"path/filepath"
	"strings"
	"testing"
	"time"

	"github.com/Justype/condatainer/catalog"
	"github.com/Justype/condatainer/internal/config"
)

// writeRecipe creates a one-recipe collection and returns its root.
func writeRecipe(t *testing.T, nameVersion, content string) string {
	t.Helper()
	root := t.TempDir()
	path := filepath.Join(root, "recipes", filepath.FromSlash(nameVersion))
	if err := os.MkdirAll(filepath.Dir(path), 0o775); err != nil {
		t.Fatalf("failed to create recipes dir: %v", err)
	}
	if err := os.WriteFile(path, []byte(content), 0o664); err != nil {
		t.Fatalf("failed to write recipe: %v", err)
	}
	return root
}

// setTestSource points the catalog at a single filesystem collection.
func setTestSource(t *testing.T, root string) {
	t.Helper()
	oldSources := config.Global.Sources
	config.Global.Sources = []catalog.Spec{{Name: "test", Base: root}}
	config.ResetCatalog()
	t.Cleanup(func() {
		config.Global.Sources = oldSources
		config.ResetCatalog()
	})
}

func TestParseScriptMetadata_RequiresTTY(t *testing.T) {
	// Create a temp script with an #INPUT: declaration
	tmp, err := os.CreateTemp("", "script-*.sh")
	if err != nil {
		t.Fatalf("failed to create temp file: %v", err)
	}
	defer os.Remove(tmp.Name())

	content := "#INPUT:Please enter something\n"
	if _, err := tmp.WriteString(content); err != nil {
		t.Fatalf("failed to write to temp file: %v", err)
	}
	tmp.Close()

	base := &BuildObject{
		spec:        Spec{Image: ImageSpec{Name: "foo/bar"}},
		buildSource: tmp.Name(),
	}

	err = base.parseScriptMetadata(context.Background())
	if err == nil {
		t.Fatalf("expected error when interactive prompts are present but no TTY is available")
	}
}

func TestParseScriptMetadata_NoInteractive(t *testing.T) {
	// Create a temp script with no #INPUT: declaration
	tmp, err := os.CreateTemp("", "script-*.sh")
	if err != nil {
		t.Fatalf("failed to create temp file: %v", err)
	}
	defer os.Remove(tmp.Name())

	content := "#DEP:foo/1.0\necho hi\n"
	if _, err := tmp.WriteString(content); err != nil {
		t.Fatalf("failed to write to temp file: %v", err)
	}
	tmp.Close()

	base := &BuildObject{
		spec:        Spec{Image: ImageSpec{Name: "foo/bar"}},
		buildSource: tmp.Name(),
	}

	err = base.parseScriptMetadata(context.Background())
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}
	if len(base.spec.Dependencies) == 0 {
		t.Fatalf("expected dependencies to be parsed")
	}
}

func TestParseScriptMetadata_NcpusFromSlurm(t *testing.T) {
	// Create a temp script with SLURM cpus-per-task directive
	tmp, err := os.CreateTemp("", "script-*.sh")
	if err != nil {
		t.Fatalf("failed to create temp file: %v", err)
	}
	defer os.Remove(tmp.Name())

	content := "#!/bin/bash\n#SBATCH --cpus-per-task=8\necho hi\n"
	if _, err := tmp.WriteString(content); err != nil {
		t.Fatalf("failed to write to temp file: %v", err)
	}
	tmp.Close()

	base := &BuildObject{
		spec:        Spec{Image: ImageSpec{Name: "foo/bar"}},
		buildSource: tmp.Name(),
	}

	err = base.parseScriptMetadata(context.Background())
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}
	if base.effectiveNcpus() != 8 {
		t.Fatalf("expected effectiveNcpus=8 from SLURM directive, got %d", base.effectiveNcpus())
	}
}

func TestParseScriptMetadata_NcpusFromPBS(t *testing.T) {
	// Create a temp script with PBS ncpus directive
	tmp, err := os.CreateTemp("", "script-*.sh")
	if err != nil {
		t.Fatalf("failed to create temp file: %v", err)
	}
	defer os.Remove(tmp.Name())

	content := "#!/bin/bash\n#PBS -l select=1:ncpus=16\necho hi\n"
	if _, err := tmp.WriteString(content); err != nil {
		t.Fatalf("failed to write to temp file: %v", err)
	}
	tmp.Close()

	base := &BuildObject{
		spec:        Spec{Image: ImageSpec{Name: "foo/bar"}},
		buildSource: tmp.Name(),
	}

	err = base.parseScriptMetadata(context.Background())
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}
	if base.effectiveNcpus() != 16 {
		t.Fatalf("expected effectiveNcpus=16 from PBS directive, got %d", base.effectiveNcpus())
	}
}

func TestNewBuildObject_DoesNotParseInteractiveWhenInstalled(t *testing.T) {
	imagesDir := t.TempDir()

	// Point the catalog at a collection providing a recipe with an #INPUT:
	// prompt, so parsing it would block on a prompt if it were attempted.
	setTestSource(t, writeRecipe(t, "cellranger/8.0.1",
		"#!/bin/bash\n#INPUT:Please enter the license key\n"))

	// Create an overlay file in imagesDir to simulate already-installed overlay
	nameVersion := "cellranger/8.0.1"
	fileName := strings.ReplaceAll(catalog.Normalize(nameVersion), "/", "--") + ".sqf"
	overlayPath := filepath.Join(imagesDir, fileName)
	if err := os.WriteFile(overlayPath, []byte{}, 0o664); err != nil {
		t.Fatalf("failed to create overlay file: %v", err)
	}

	// Re-init data paths so the test's extra base dir is picked up
	config.InitDataPaths()

	// Call NewBuildObject: should return without attempting to parse the recipe inputs
	bo, err := NewBuildObject(context.Background(), nameVersion, false, imagesDir, false)
	if err != nil {
		t.Fatalf("NewBuildObject returned error: %v", err)
	}
	if bo == nil {
		t.Fatalf("expected non-nil BuildObject")
	}
}

func TestNewBuildObject_ErrorsWhenBuildLockExists(t *testing.T) {
	imagesDir := t.TempDir()

	// Point the catalog at a collection providing a recipe with an #INPUT:
	// prompt, so parsing it would block on a prompt if it were attempted.
	setTestSource(t, writeRecipe(t, "cellranger/8.0.1",
		"#!/bin/bash\n#INPUT:Please enter the license key\n"))

	// Simulate a build in progress: create a live JSON lock file next to the target image.
	// Use the current process PID and hostname so isBuildLockStale() treats it as active.
	// The build-in-progress guard checks base.tgt.Lock = target path + ".lock",
	// which lives in imagesDir.
	nameVersion := "cellranger/8.0.1"
	sqfName := strings.ReplaceAll(catalog.Normalize(nameVersion), "/", "--") + ".sqf"
	lockPath := filepath.Join(imagesDir, sqfName+".lock")
	liveLock := BuildLockInfo{
		Runner:    "local",
		Node:      shortHostname(),
		PID:       os.Getpid(), // this process is definitely alive
		CreatedAt: time.Now().Format(time.RFC3339),
	}
	lockData, err := json.Marshal(liveLock)
	if err != nil {
		t.Fatalf("failed to marshal lock: %v", err)
	}
	if err := os.WriteFile(lockPath, lockData, 0o664); err != nil {
		t.Fatalf("failed to create lock file: %v", err)
	}
	defer os.Remove(lockPath)

	// Re-init data paths so the test's extra base dir is picked up
	config.InitDataPaths()

	// Call NewBuildObject: should return an error with lock details.
	_, err = NewBuildObject(context.Background(), nameVersion, false, imagesDir, false)
	if err == nil {
		t.Fatal("expected error when build lock file exists, got nil")
	}
	if !strings.Contains(err.Error(), "build lock found") {
		t.Fatalf("unexpected error message: %v", err)
	}
}

// An external build's name decides its /cnt/<name> prefix and, through
// key.Role, which dependencies count toward its equivalence. #TARGET: is where
// that name comes from; -p only says where the file goes.
func TestExternalSourceTakesItsNameFromTarget(t *testing.T) {
	dir := t.TempDir()
	src := filepath.Join(dir, "build.sh")
	if err := os.WriteFile(src, []byte(
		"#!/usr/bin/env bash\n#TYPE: data\n#TARGET: star/2.7.11b/index\necho hi\n"), 0o644); err != nil {
		t.Fatal(err)
	}

	b, err := FromExternalSource(t.Context(), filepath.Join(dir, "idx"), src, false, dir, false)
	if err != nil {
		t.Fatalf("FromExternalSource: %v", err)
	}
	if got := b.spec.Image.Name; got != "star/2.7.11b/index" {
		t.Errorf("name = %q, want the declared #TARGET:", got)
	}
	if got := b.spec.Image.Prefix; got != "/cnt/star/2.7.11b/index" {
		t.Errorf("prefix = %q, want it derived from #TARGET:", got)
	}
	// The file keeps the -p name: the two namespaces are deliberately unrelated.
	if want := filepath.Join(dir, "idx.sqf"); b.tgt.Path != want {
		t.Errorf("target = %q, want %q", b.tgt.Path, want)
	}
}

// A definition names itself with #TARGET: the way a shell script does.
func TestExternalDefinitionTakesItsNameFromTarget(t *testing.T) {
	dir := t.TempDir()
	src := filepath.Join(dir, "tool.def")
	if err := os.WriteFile(src, []byte(
		"#TARGET: ubuntu24/tool\nBootstrap: docker\nFrom: ubuntu:24.04\n"), 0o644); err != nil {
		t.Fatal(err)
	}

	b, err := FromExternalSource(t.Context(), filepath.Join(dir, "idx"), src, true, dir, false)
	if err != nil {
		t.Fatalf("FromExternalSource: %v", err)
	}
	if got := b.spec.Image.Name; got != "ubuntu24/tool" {
		t.Errorf("name = %q, want the declared #TARGET:", got)
	}
}

// WithName decides the name and leaves the file where the path put it, over both
// a #TARGET: and the basename.
func TestWithNameOverridesTheNameNotTheFile(t *testing.T) {
	dir := t.TempDir()
	src := filepath.Join(dir, "build.sh")
	if err := os.WriteFile(src, []byte(
		"#!/usr/bin/env bash\n#TARGET: declared/1.0\necho hi\n"), 0o644); err != nil {
		t.Fatal(err)
	}

	b, err := FromExternalSource(t.Context(), filepath.Join(dir, "a--b"), src, false, dir, false, WithName("star=2.7"))
	if err != nil {
		t.Fatalf("FromExternalSource: %v", err)
	}
	if got := b.spec.Image.Name; got != "star/2.7" {
		t.Errorf("external name = %q, want the WithName value normalized", got)
	}
	if want := filepath.Join(dir, "a--b.sqf"); b.tgt.Path != want {
		t.Errorf("external target = %q, want %q", b.tgt.Path, want)
	}

	c, err := NewCondaObjectWithSource("plain@1", "", dir, false, WithName("star/2.7"))
	if err != nil {
		t.Fatalf("NewCondaObjectWithSource: %v", err)
	}
	if got := c.spec.Image.Name; got != "star/2.7" {
		t.Errorf("conda name = %q, want star/2.7", got)
	}
	if want := filepath.Join(dir, "plain@1.sqf"); c.tgt.Path != want {
		t.Errorf("conda target = %q, want %q", c.tgt.Path, want)
	}
}

// Without #TARGET: the name is the -p basename, which is the historical
// behaviour and stays.
func TestExternalSourceFallsBackToThePrefixBasename(t *testing.T) {
	dir := t.TempDir()
	src := filepath.Join(dir, "build.sh")
	if err := os.WriteFile(src, []byte("#!/usr/bin/env bash\necho hi\n"), 0o644); err != nil {
		t.Fatal(err)
	}

	b, err := FromExternalSource(t.Context(), filepath.Join(dir, "demo"), src, false, dir, false)
	if err != nil {
		t.Fatalf("FromExternalSource: %v", err)
	}
	if got := b.spec.Image.Name; got != "demo" {
		t.Errorf("name = %q, want the -p basename", got)
	}
}

// Declaring a dependency without naming yourself leaves the dependency's role —
// and so the artifact's equivalence — decided by where the file was written.
func TestExternalSourceRefusesDependenciesWithoutTarget(t *testing.T) {
	dir := t.TempDir()
	src := filepath.Join(dir, "build.sh")
	if err := os.WriteFile(src, []byte(
		"#!/usr/bin/env bash\n#TYPE: data\n#DEP: samtools/1.21\necho hi\n"), 0o644); err != nil {
		t.Fatal(err)
	}

	_, err := FromExternalSource(t.Context(), filepath.Join(dir, "idx"), src, false, dir, false)
	if err == nil {
		t.Fatal("a #DEP: without #TARGET: was accepted")
	}
	if !strings.Contains(err.Error(), "#TARGET:") {
		t.Errorf("refusal does not name the remedy: %v", err)
	}
}

// The data-only rule is not a catalog-recipe rule: an external build reaches it
// too, whatever it is named and wherever -p puts it.
func TestExternalSourceRefusesDependenciesOnAnApp(t *testing.T) {
	for _, tc := range []struct{ name, content string }{
		{"WithTarget", "#!/usr/bin/env bash\n#TYPE: app\n#TARGET: mytool/1.0\n#DEP: samtools/1.21\n"},
		{"WithoutTarget", "#!/usr/bin/env bash\n#TYPE: app\n#DEP: samtools/1.21\n"},
		{"TypeOmittedDefaultsToApp", "#!/usr/bin/env bash\n#TARGET: mytool/1.0\n#DEP: samtools/1.21\n"},
	} {
		t.Run(tc.name, func(t *testing.T) {
			dir := t.TempDir()
			src := filepath.Join(dir, "build.sh")
			if err := os.WriteFile(src, []byte(tc.content), 0o644); err != nil {
				t.Fatal(err)
			}
			_, err := FromExternalSource(t.Context(), filepath.Join(dir, "idx"), src, false, dir, false)
			if err == nil {
				t.Fatal("an app declaring #DEP: was accepted")
			}
			if !strings.Contains(err.Error(), "only data has build dependencies") {
				t.Errorf("refusal does not name the rule: %v", err)
			}
		})
	}
}

// A path names neither an artifact nor an identity, so nothing could re-resolve
// it elsewhere or regenerate a key from it.
func TestExternalSourceRefusesAPathDependency(t *testing.T) {
	dir := t.TempDir()
	src := filepath.Join(dir, "build.sh")
	if err := os.WriteFile(src, []byte(
		"#!/usr/bin/env bash\n#TYPE: data\n#TARGET: grch38/index\n#DEP: overlays/tool.sqf\n"), 0o644); err != nil {
		t.Fatal(err)
	}

	_, err := FromExternalSource(t.Context(), filepath.Join(dir, "idx"), src, false, dir, false)
	if err == nil {
		t.Fatal("a path #DEP: was accepted")
	}
	if !strings.Contains(err.Error(), "not an overlay path") {
		t.Errorf("refusal does not name the rule: %v", err)
	}
}
