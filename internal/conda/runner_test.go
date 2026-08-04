package conda

import (
	"context"
	"os"
	"path/filepath"
	"reflect"
	"strings"
	"testing"
)

func TestInstallCreatesFreshPrefixWithNoRCAndExplicitChannels(t *testing.T) {
	env, argsPath := fakeEnvironment(t)
	if err := env.Install(context.Background(), []string{"-y", "-c", "pytorch", "python"}, IO{}); err != nil {
		t.Fatal(err)
	}
	got := readArgumentLines(t, argsPath)
	want := []string{"--no-rc", "install", "-c", "pytorch", "-c", "conda-forge", "-c", "bioconda", "-y", "-c", "pytorch", "python"}
	if !reflect.DeepEqual(got, want) {
		t.Fatalf("micromamba args = %v, want %v", got, want)
	}
	channels, err := ReadChannels(env.CondarcPath)
	if err != nil {
		t.Fatal(err)
	}
	if want := []string{"pytorch", "conda-forge", "bioconda"}; !reflect.DeepEqual(channels, want) {
		t.Fatalf("channels = %v, want %v", channels, want)
	}
}

func TestInstallCreatesEmptyPrefixWithoutInjectingPackage(t *testing.T) {
	env, argsPath := fakeEnvironment(t)
	env.DefaultChannels = []string{"internal", "conda-forge"}
	if err := env.Install(context.Background(), nil, IO{}); err != nil {
		t.Fatal(err)
	}
	if got, want := readArgumentLines(t, argsPath), []string{"--no-rc", "install", "-c", "internal", "-c", "conda-forge"}; !reflect.DeepEqual(got, want) {
		t.Fatalf("micromamba args = %v, want %v", got, want)
	}
	channels, err := ReadChannels(env.CondarcPath)
	if err != nil {
		t.Fatal(err)
	}
	if want := env.DefaultChannels; !reflect.DeepEqual(channels, want) {
		t.Fatalf("channels = %v, want %v", channels, want)
	}
}

func TestProjectChannelsFromEnvironmentFileWithoutChangingArgs(t *testing.T) {
	env, argsPath := fakeEnvironment(t)
	file := filepath.Join(t.TempDir(), "environment.yml")
	if err := os.WriteFile(file, []byte("channels:\n  - bioconda\n  - conda-forge\ndependencies:\n  - samtools\n"), 0o644); err != nil {
		t.Fatal(err)
	}
	args := []string{"-f", file, "--future-option", "value"}
	if err := env.Install(context.Background(), args, IO{}); err != nil {
		t.Fatal(err)
	}
	wantArgs := []string{"--no-rc", "install", "-c", "bioconda", "-c", "conda-forge"}
	wantArgs = append(wantArgs, args...)
	if got := readArgumentLines(t, argsPath); !reflect.DeepEqual(got, wantArgs) {
		t.Fatalf("micromamba args = %v, want %v", got, wantArgs)
	}
	channels, err := ReadChannels(env.CondarcPath)
	if err != nil {
		t.Fatal(err)
	}
	if want := []string{"bioconda", "conda-forge"}; !reflect.DeepEqual(channels, want) {
		t.Fatalf("channels = %v, want %v", channels, want)
	}
}

func TestFreshInstallPreservesExistingProjectCondarc(t *testing.T) {
	env, argsPath := fakeEnvironment(t)
	projectChannels := []string{"local", "conda-forge"}
	if err := WriteChannels(env.CondarcPath, projectChannels); err != nil {
		t.Fatal(err)
	}
	args := []string{"python", "--future-option", "value"}
	if err := env.Install(context.Background(), args, IO{}); err != nil {
		t.Fatal(err)
	}
	wantArgs := append([]string{"--rc-file", env.CondarcPath, "install"}, args...)
	if got := readArgumentLines(t, argsPath); !reflect.DeepEqual(got, wantArgs) {
		t.Fatalf("micromamba args = %v, want %v", got, wantArgs)
	}
	channels, err := ReadChannels(env.CondarcPath)
	if err != nil {
		t.Fatal(err)
	}
	if !reflect.DeepEqual(channels, projectChannels) {
		t.Fatalf("channels = %v, want %v", channels, projectChannels)
	}
}

func TestInstallUsesSavedConfigForInitializedPrefix(t *testing.T) {
	env, argsPath := fakeEnvironment(t)
	meta := filepath.Join(env.Root, "conda-meta")
	if err := os.MkdirAll(meta, 0o755); err != nil {
		t.Fatal(err)
	}
	if err := os.WriteFile(filepath.Join(meta, "history"), nil, 0o644); err != nil {
		t.Fatal(err)
	}
	if err := WriteChannels(env.CondarcPath, []string{"conda-forge"}); err != nil {
		t.Fatal(err)
	}
	if err := env.Install(context.Background(), []string{"numpy"}, IO{}); err != nil {
		t.Fatal(err)
	}
	got := readArgumentLines(t, argsPath)
	want := []string{"--rc-file", env.CondarcPath, "install", "numpy"}
	if !reflect.DeepEqual(got, want) {
		t.Fatalf("micromamba args = %v, want %v", got, want)
	}
}

func TestRunForwardsUnknownArgumentsUnchanged(t *testing.T) {
	env, argsPath := fakeEnvironment(t)
	args := []string{"update", "--future-option", "value", "numpy"}
	if err := env.Run(context.Background(), IO{}, args...); err != nil {
		t.Fatal(err)
	}
	want := append([]string{"--no-rc"}, args...)
	if got := readArgumentLines(t, argsPath); !reflect.DeepEqual(got, want) {
		t.Fatalf("micromamba args = %v, want %v", got, want)
	}
}

func TestCommandEnvSelectsMountedEnvironmentWithoutArguments(t *testing.T) {
	env, _ := fakeEnvironment(t)
	t.Setenv("CONDA_PREFIX", "/host/env")
	t.Setenv("CONDARC", "/host/.condarc")
	t.Setenv("MAMBA_NO_RC", "true")
	t.Setenv("MAMBA_TARGET_PREFIX", "/host/target")
	if err := WriteChannels(env.CondarcPath, []string{"conda-forge"}); err != nil {
		t.Fatal(err)
	}

	values := map[string]string{}
	for _, entry := range env.commandEnv() {
		name, value, found := strings.Cut(entry, "=")
		if found {
			values[name] = value
		}
	}
	for _, name := range []string{"CONDA_PREFIX", "MAMBA_ROOT_PREFIX"} {
		if values[name] != env.Root {
			t.Errorf("%s = %q, want %q", name, values[name], env.Root)
		}
	}
	if _, found := values["MAMBA_TARGET_PREFIX"]; found {
		t.Error("MAMBA_TARGET_PREFIX should be unset")
	}
	for _, name := range []string{"CONDARC", "MAMBA_NO_RC"} {
		if _, found := values[name]; found {
			t.Errorf("ambient %s should not be inherited", name)
		}
	}
}

func TestCommandEnvDisablesAmbientRcWhenOverlayHasNone(t *testing.T) {
	env, _ := fakeEnvironment(t)
	t.Setenv("CONDARC", "/host/.condarc")

	values := map[string]string{}
	for _, entry := range env.commandEnv() {
		name, value, found := strings.Cut(entry, "=")
		if found {
			values[name] = value
		}
	}
	if _, found := values["CONDARC"]; found {
		t.Error("ambient CONDARC should not be inherited")
	}
	if _, found := values["MAMBA_NO_RC"]; found {
		t.Error("ambient MAMBA_NO_RC should not be inherited")
	}
}

func fakeEnvironment(t *testing.T) (*Environment, string) {
	t.Helper()
	root := t.TempDir()
	argsPath := filepath.Join(t.TempDir(), "args")
	binary := filepath.Join(t.TempDir(), "micromamba")
	script := "#!/bin/sh\nprintf '%s\\n' \"$@\" > \"$MAMBA_TEST_ARGS\"\n"
	if err := os.WriteFile(binary, []byte(script), 0o755); err != nil {
		t.Fatal(err)
	}
	t.Setenv("MAMBA_TEST_ARGS", argsPath)
	return &Environment{
		Root:            root,
		CondarcPath:     filepath.Join(root, ".condarc"),
		PinnedPath:      filepath.Join(root, "conda-meta", "pinned"),
		Micromamba:      binary,
		DefaultChannels: []string{"conda-forge", "bioconda"},
		Writable:        true,
	}, argsPath
}

func readArgumentLines(t *testing.T, path string) []string {
	t.Helper()
	data, err := os.ReadFile(path)
	if err != nil {
		t.Fatal(err)
	}
	return strings.Split(strings.TrimSpace(string(data)), "\n")
}
