package container

import (
	"slices"
	"strings"
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

// prefixes answers with a fixed prefix per path, standing in for the runtime
// metadata a real overlay carries.
func prefixes(byPath map[string]string) func(string) string {
	return func(path string) string { return byPath[path] }
}

// Two builds of one name record one prefix. Overlays are disjoint subtrees, so
// the later mount takes the whole subtree and the earlier contributes nothing.
func TestDistinctPrefixesRefusesTwoImagesClaimingOneSubtree(t *testing.T) {
	overlays := []string{"/images/tool--1.0.sqf", "/project/overlays/tool.sqf"}
	err := distinctPrefixes(overlays, prefixes(map[string]string{
		"/images/tool--1.0.sqf":      "/cnt/tool/1.0",
		"/project/overlays/tool.sqf": "/cnt/tool/1.0",
	}))
	if err == nil {
		t.Fatal("two overlays claiming one prefix were accepted")
	}
	for _, want := range []string{"/images/tool--1.0.sqf", "/project/overlays/tool.sqf", "/cnt/tool/1.0"} {
		if !strings.Contains(err.Error(), want) {
			t.Errorf("error %q does not name %q", err, want)
		}
	}
}

// Different names are different subtrees, which is the ordinary case.
func TestDistinctPrefixesAllowsDifferentNames(t *testing.T) {
	overlays := []string{"/images/samtools--1.22.sqf", "/images/bcftools--1.20.sqf"}
	err := distinctPrefixes(overlays, prefixes(map[string]string{
		"/images/samtools--1.22.sqf": "/cnt/samtools/1.22",
		"/images/bcftools--1.20.sqf": "/cnt/bcftools/1.20",
	}))
	if err != nil {
		t.Fatalf("distinct names were refused: %v", err)
	}
}

// A base, an OS image, and anything without readable metadata record no prefix
// and claim no subtree, so any number of them mount together.
func TestDistinctPrefixesIgnoresImagesWithNoPrefix(t *testing.T) {
	overlays := []string{"/images/os.sqf", "/images/base.sqf", "/images/unreadable.sqf"}
	err := distinctPrefixes(overlays, prefixes(nil))
	if err != nil {
		t.Fatalf("prefix-less images were refused: %v", err)
	}
}

// The same file named twice is redundant, not a collision.
func TestDistinctPrefixesAllowsOneFileNamedTwice(t *testing.T) {
	overlays := []string{"/images/tool--1.0.sqf", "/images/tool--1.0.sqf:ro"}
	err := distinctPrefixes(overlays, prefixes(map[string]string{
		"/images/tool--1.0.sqf": "/cnt/tool/1.0",
	}))
	if err != nil {
		t.Fatalf("one file named twice was refused: %v", err)
	}
}
