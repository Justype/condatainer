package conda

import (
	"os"
	"path/filepath"
	"slices"
	"strings"
	"testing"
)

func TestResolveEnvironmentReadAndWrite(t *testing.T) {
	root := t.TempDir()
	if err := os.MkdirAll(filepath.Join(root, "conda-meta"), 0o755); err != nil {
		t.Fatal(err)
	}
	if err := os.WriteFile(filepath.Join(root, "conda-meta", "history"), nil, 0o644); err != nil {
		t.Fatal(err)
	}
	t.Setenv("IN_CONDATAINER", "1")
	t.Setenv("CNT_CONDA_ROOT", root)
	t.Setenv("CNT_CONDA_WRITABLE", "0")

	env, err := ResolveEnvironment(false)
	if err != nil {
		t.Fatal(err)
	}
	if env.Writable {
		t.Fatal("read-only environment reported writable")
	}
	if _, err := ResolveEnvironment(true); err == nil || !strings.Contains(err.Error(), "read-only") {
		t.Fatalf("write resolution error = %v", err)
	}

	t.Setenv("CNT_CONDA_WRITABLE", "1")
	env, err = ResolveEnvironment(true)
	if err != nil {
		t.Fatal(err)
	}
	if !env.Writable {
		t.Fatal("writable environment reported read-only")
	}
}

func TestResolveEnvironmentRequiresContainer(t *testing.T) {
	t.Setenv("IN_CONDATAINER", "")
	t.Setenv("CNT_CONDA_ROOT", t.TempDir())
	if _, err := ResolveEnvironment(false); err == nil || !strings.Contains(err.Error(), "inside CondaTainer") {
		t.Fatalf("error = %v", err)
	}
}

func TestResolveEnvironmentReportsUninitialized(t *testing.T) {
	root := filepath.Join(t.TempDir(), "env")
	t.Setenv("IN_CONDATAINER", "1")
	t.Setenv("CNT_CONDA_ROOT", root)
	t.Setenv("CNT_CONDA_WRITABLE", "1")

	if _, err := ResolveEnvironment(false); err == nil || !strings.Contains(err.Error(), "not initialized") {
		t.Fatalf("read resolution error = %v", err)
	}
	if _, err := os.Stat(root); !os.IsNotExist(err) {
		t.Fatalf("read resolution created environment root: %v", err)
	}
}

func TestResolveEnvironmentReportsExistingButUninitializedRoot(t *testing.T) {
	root := t.TempDir()
	t.Setenv("IN_CONDATAINER", "1")
	t.Setenv("CNT_CONDA_ROOT", root)
	t.Setenv("CNT_CONDA_WRITABLE", "1")

	if _, err := ResolveEnvironment(false); err == nil || !strings.Contains(err.Error(), "not initialized") {
		t.Fatalf("read resolution error = %v", err)
	}
}

func TestResolveInstallEnvironmentCreatesRoot(t *testing.T) {
	root := filepath.Join(t.TempDir(), "env")
	t.Setenv("IN_CONDATAINER", "1")
	t.Setenv("CNT_CONDA_ROOT", root)
	t.Setenv("CNT_CONDA_WRITABLE", "1")
	t.Setenv("CNT_CHANNELS", "unrelated")
	t.Setenv("CNT_CONDA_CHANNELS", "internal|conda-forge|internal")

	env, err := ResolveInstallEnvironment()
	if err != nil {
		t.Fatal(err)
	}
	if env.Root != root {
		t.Fatalf("root = %q, want %q", env.Root, root)
	}
	if got, want := env.DefaultChannels, []string{"internal", "conda-forge"}; !slices.Equal(got, want) {
		t.Fatalf("default channels = %v, want %v", got, want)
	}
	if info, err := os.Stat(root); err != nil || !info.IsDir() {
		t.Fatalf("environment root was not created: info=%v err=%v", info, err)
	}
}

func TestResolveInstallEnvironmentRequiresWritableOverlay(t *testing.T) {
	root := filepath.Join(t.TempDir(), "env")
	t.Setenv("IN_CONDATAINER", "1")
	t.Setenv("CNT_CONDA_ROOT", root)
	t.Setenv("CNT_CONDA_WRITABLE", "0")

	if _, err := ResolveInstallEnvironment(); err == nil || !strings.Contains(err.Error(), "read-only") {
		t.Fatalf("install resolution error = %v", err)
	}
	if _, err := os.Stat(root); !os.IsNotExist(err) {
		t.Fatalf("read-only resolution created environment root: %v", err)
	}
}

func TestPrefixState(t *testing.T) {
	root := t.TempDir()
	env := &Environment{Root: root}
	state, err := env.State()
	if err != nil || state != PrefixFresh {
		t.Fatalf("fresh state = %v, err = %v", state, err)
	}
	if err := os.Mkdir(filepath.Join(root, "conda-meta"), 0o755); err != nil {
		t.Fatal(err)
	}
	state, err = env.State()
	if err != nil || state != PrefixPartial {
		t.Fatalf("partial state = %v, err = %v", state, err)
	}
	if err := os.WriteFile(filepath.Join(root, "conda-meta", "history"), nil, 0o644); err != nil {
		t.Fatal(err)
	}
	state, err = env.State()
	if err != nil || state != PrefixInitialized {
		t.Fatalf("initialized state = %v, err = %v", state, err)
	}
}
