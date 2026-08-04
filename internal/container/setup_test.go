package container

import (
	"slices"
	"testing"

	"github.com/Justype/condatainer/internal/config"
)

func TestBuildEnvironmentCondaMarkers(t *testing.T) {
	previousChannels := config.Global.Build.Channels
	config.Global.Build.Channels = []string{"internal", "conda-forge"}
	t.Cleanup(func() { config.Global.Build.Channels = previousChannels })

	readOnly, _, _ := buildEnvironment(nil, "env.img", SetupConfig{})
	if !slices.Contains(readOnly, "CNT_CONDA_ROOT=/cnt_env") {
		t.Fatal("read-only environment is missing CNT_CONDA_ROOT")
	}
	if !slices.Contains(readOnly, "CNT_CONDA_WRITABLE=0") {
		t.Fatal("read-only environment is missing its writable indicator")
	}
	if !slices.Contains(readOnly, "CNT_CONDA_CHANNELS=internal|conda-forge") {
		t.Fatal("environment is missing the configured channels")
	}

	writable, _, _ := buildEnvironment(nil, "env.img", SetupConfig{
		WritableImg: true,
		EnvSettings: []string{
			"CNT_CONDA_ROOT=/wrong", "CNT_CONDA_WRITABLE=0",
			"CONDA_PREFIX=/wrong", "MAMBA_ROOT_PREFIX=/wrong",
		},
	})
	for _, want := range []string{
		"CNT_CONDA_ROOT=/cnt_env", "CNT_CONDA_WRITABLE=1",
		"CONDA_PREFIX=/cnt_env", "MAMBA_ROOT_PREFIX=/cnt_env",
	} {
		if !slices.Contains(writable, want) {
			t.Fatalf("writable environment is missing %s", want)
		}
	}
	wantTail := []string{
		"CNT_CONDA_ROOT=/cnt_env", "CONDA_PREFIX=/cnt_env",
		"MAMBA_ROOT_PREFIX=/cnt_env", "CNT_CONDA_WRITABLE=1",
	}
	if got := writable[len(writable)-4:]; !slices.Equal(got, wantTail) {
		t.Fatalf("runtime markers must override user settings; tail = %v", got)
	}
}
