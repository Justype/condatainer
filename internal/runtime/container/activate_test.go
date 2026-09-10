package container

import (
	"os"
	"os/exec"
	"path/filepath"
	"strings"
	"testing"

	"github.com/Justype/condatainer/catalog"
	"github.com/Justype/condatainer/internal/artifact/meta"
)

func TestActivationScriptEmptyWithNothingToActivate(t *testing.T) {
	if got := ActivationScript(nil, false, ActivationAll); got != "" {
		t.Errorf("ActivationScript(nil, false, ActivationAll) = %q, want empty", got)
	}
}

func TestActivationScriptNoneSourcesNothing(t *testing.T) {
	requireSquashfsTools(t)
	dir := t.TempDir()
	appSqf := filepath.Join(dir, "myapp.sqf")
	packRuntimeSqf(t, appSqf, meta.Runtime{
		SchemaVersion: meta.SchemaVersion,
		Name:          "myapp/1.0",
		Type:          catalog.TypeApp,
		Platform:      meta.NativePlatform(),
		Prefix:        "/cnt/myapp",
	})

	if got := ActivationScript([]string{appSqf}, true, ActivationNone); got != "" {
		t.Errorf("ActivationScript(..., ActivationNone) = %q, want empty", got)
	}
}

func TestActivationScriptEnvSourcesOnlyCntEnv(t *testing.T) {
	requireSquashfsTools(t)
	dir := t.TempDir()
	appSqf := filepath.Join(dir, "myapp.sqf")
	packRuntimeSqf(t, appSqf, meta.Runtime{
		SchemaVersion: meta.SchemaVersion,
		Name:          "myapp/1.0",
		Type:          catalog.TypeApp,
		Platform:      meta.NativePlatform(),
		Prefix:        "/cnt/myapp",
	})

	got := ActivationScript([]string{appSqf}, true, ActivationEnv)
	if strings.Contains(got, "'/cnt/myapp'/etc/conda/activate.d") {
		t.Errorf("ActivationEnv sourced an app block, want only /cnt_env:\n%s", got)
	}
	if !strings.Contains(got, "'/cnt_env'/etc/conda/activate.d") {
		t.Errorf("script missing /cnt_env block:\n%s", got)
	}
}

func TestActivationScriptIncludesAppPrefixAndEnv(t *testing.T) {
	requireSquashfsTools(t)
	dir := t.TempDir()
	appSqf := filepath.Join(dir, "myapp.sqf")
	packRuntimeSqf(t, appSqf, meta.Runtime{
		SchemaVersion: meta.SchemaVersion,
		Name:          "myapp/1.0",
		Type:          catalog.TypeApp,
		Platform:      meta.NativePlatform(),
		Prefix:        "/cnt/myapp",
	})
	envSqf := filepath.Join(dir, "env.sqf")
	packRuntimeSqf(t, envSqf, envRuntime())

	got := ActivationScript([]string{appSqf, envSqf}, true, ActivationAll)

	if !strings.Contains(got, "'/cnt/myapp'/etc/conda/activate.d") {
		t.Errorf("script missing app prefix block:\n%s", got)
	}
	if !strings.Contains(got, "'/cnt_env'/etc/conda/activate.d") {
		t.Errorf("script missing /cnt_env block for a mounted conda environment:\n%s", got)
	}
}

// A writable .img and its paired env.sqf snapshot both claim EnvPrefix, but
// Apptainer merges them into the one physical /cnt_env directory at mount
// time — the generated script must source that directory's activate.d only
// once, not once per overlay that claims the prefix.
func TestActivationScriptDedupesPairedImgAndSnapshot(t *testing.T) {
	requireSquashfsTools(t)
	dir := t.TempDir()
	envSqf := filepath.Join(dir, "env.sqf")
	packRuntimeSqf(t, envSqf, envRuntime())
	// imgContribution only reads an optional <path>.env sidecar, so a bare
	// path ending in .img exercises it without a real ext3 image.
	img := filepath.Join(dir, "env.img")

	got := ActivationScript([]string{envSqf, img}, true, ActivationAll)

	if n := strings.Count(got, "if [ -d '/cnt_env'/etc/conda/activate.d"); n != 1 {
		t.Errorf("want exactly one /cnt_env activate.d block, got %d:\n%s", n, got)
	}
}

func TestMMHelperScriptEmptyWithoutEnvMounted(t *testing.T) {
	if got := MMHelperScript(false); got != "" {
		t.Errorf("MMHelperScript(false) = %q, want empty", got)
	}
}

func TestMMHelperScriptDefinesAndExportsMM(t *testing.T) {
	got := MMHelperScript(true)
	if !strings.Contains(got, "export -f mm") {
		t.Errorf("MMHelperScript(true) = %q, want it to export the mm function", got)
	}
	if !strings.Contains(got, "reactivate --shell bash") {
		t.Errorf("MMHelperScript(true) = %q, want it to force bash reactivate syntax", got)
	}
}

// TestActivationScriptRunsWithScopedPrefix runs the generated script for real
// through bash against two fake overlay prefixes, each with its own
// activate.d script that reads CONDA_PREFIX the way libxml2's real one does
// (internal/libexec's own provisioned toolchain ships exactly this script) —
// proving each block sees its own prefix, not whichever ran last or first.
func TestActivationScriptRunsWithScopedPrefix(t *testing.T) {
	requireSquashfsTools(t)
	root := t.TempDir()

	writePrefixHook := func(name, marker string) string {
		prefix := filepath.Join(root, name)
		activateDir := filepath.Join(prefix, "etc", "conda", "activate.d")
		if err := os.MkdirAll(activateDir, 0o755); err != nil {
			t.Fatal(err)
		}
		script := "export CNT_TEST_" + marker + "=\"$CONDA_PREFIX\"\n"
		if err := os.WriteFile(filepath.Join(activateDir, "hook.sh"), []byte(script), 0o755); err != nil {
			t.Fatal(err)
		}
		return prefix
	}

	appAPrefix := writePrefixHook("appA", "A")
	appBPrefix := writePrefixHook("appB", "B")

	dir := t.TempDir()
	appASqf := filepath.Join(dir, "appA.sqf")
	appBSqf := filepath.Join(dir, "appB.sqf")
	packRuntimeSqf(t, appASqf, meta.Runtime{
		SchemaVersion: meta.SchemaVersion, Name: "appA/1.0", Type: catalog.TypeApp,
		Platform: meta.NativePlatform(), Prefix: appAPrefix,
	})
	packRuntimeSqf(t, appBSqf, meta.Runtime{
		SchemaVersion: meta.SchemaVersion, Name: "appB/1.0", Type: catalog.TypeApp,
		Platform: meta.NativePlatform(), Prefix: appBPrefix,
	})

	script := ActivationScript([]string{appASqf, appBSqf}, false, ActivationAll)
	script += "echo \"A=$CNT_TEST_A B=$CNT_TEST_B\"\n"

	out, err := exec.Command("bash", "-c", script).CombinedOutput()
	if err != nil {
		t.Fatalf("bash -c script: %v\n%s", err, out)
	}
	want := "A=" + appAPrefix + " B=" + appBPrefix + "\n"
	if string(out) != want {
		t.Errorf("output = %q, want %q", out, want)
	}
}
