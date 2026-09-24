package container

import (
	"os"
	"path/filepath"
	"slices"
	"strings"
	"testing"

	"github.com/Justype/condatainer/catalog"
	"github.com/Justype/condatainer/internal/artifact/meta"
	"github.com/Justype/condatainer/internal/config"
)

func TestBuildEnvironmentCondaMarkers(t *testing.T) {
	previousChannels := config.Global.Build.Channels
	config.Global.Build.Channels = []string{"internal", "conda-forge"}
	t.Cleanup(func() { config.Global.Build.Channels = previousChannels })

	readOnly, _, _ := buildEnvironment(nil, "env.img", true, SetupConfig{})
	if !slices.Contains(readOnly, "CNT_CONDA_ROOT=/cnt_env") {
		t.Fatal("read-only environment is missing CNT_CONDA_ROOT")
	}
	if !slices.Contains(readOnly, "CNT_CONDA_WRITABLE=0") {
		t.Fatal("read-only environment is missing its writable indicator")
	}
	if !slices.Contains(readOnly, "CNT_CONDA_CHANNELS=internal|conda-forge") {
		t.Fatal("environment is missing the configured channels")
	}

	writable, _, _ := buildEnvironment(nil, "env.img", true, SetupConfig{
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

// A bare env-typed .sqf with no paired .img (envMounted true, lastImg "")
// still gets the conda markers, but never a writable one — there is no .img
// to write to, regardless of cfg.WritableImg.
func TestBuildEnvironmentSqfOnlyEnvMountIsNeverWritable(t *testing.T) {
	env, _, _ := buildEnvironment(nil, "", true, SetupConfig{WritableImg: true})
	if !slices.Contains(env, "CNT_CONDA_ROOT=/cnt_env") {
		t.Fatal("sqf-only environment is missing CNT_CONDA_ROOT")
	}
	if !slices.Contains(env, "CNT_CONDA_WRITABLE=0") {
		t.Fatal("sqf-only environment must report CNT_CONDA_WRITABLE=0 even with WritableImg set")
	}

	notMounted, _, _ := buildEnvironment(nil, "", false, SetupConfig{})
	if slices.ContainsFunc(notMounted, func(v string) bool { return strings.HasPrefix(v, "CNT_CONDA_ROOT=") }) {
		t.Fatal("no conda environment mounted must not set CNT_CONDA_ROOT")
	}
}

func TestNestedRootEnv(t *testing.T) {
	tests := []struct{ name, inherited, root, want string }{
		{"detected root is handed over", "", "/opt/cnt", "CNT_ROOT=/opt/cnt"},
		{"an existing CNT_ROOT is left alone", "/x", "/x", ""},
		{"no root, nothing to set", "", "", ""},
	}
	for _, tt := range tests {
		if got := nestedRootEnv(tt.inherited, tt.root); got != tt.want {
			t.Errorf("%s: got %q, want %q", tt.name, got, tt.want)
		}
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

// An .img has no identity of its own and is expected to sit on top of an
// env-typed .sqf, so it never participates in the collision check even though
// it contributes the same EnvPrefix for environment-variable purposes.
func TestDistinctPrefixesIgnoresImg(t *testing.T) {
	overlays := []string{"/proj/env.sqf", "/proj/env.img"}
	err := distinctPrefixes(overlays, prefixes(map[string]string{
		"/proj/env.sqf": meta.EnvPrefix,
		"/proj/env.img": meta.EnvPrefix,
	}))
	if err != nil {
		t.Fatalf("env.sqf + env.img was refused: %v", err)
	}
}

// Two environment snapshots together is still refused: there is no way to
// tell which one is meant.
func TestDistinctPrefixesRefusesTwoEnvSnapshots(t *testing.T) {
	overlays := []string{"/proj/env.sqf", "/proj/env-alice.sqf"}
	err := distinctPrefixes(overlays, prefixes(map[string]string{
		"/proj/env.sqf":       meta.EnvPrefix,
		"/proj/env-alice.sqf": meta.EnvPrefix,
	}))
	if err == nil {
		t.Fatal("two environment snapshots were accepted")
	}
}

// orderOverlays puts a writable .img last and, immediately beneath it, the one
// env-typed .sqf present — regardless of where either appeared originally.
func TestOrderOverlaysPlacesEnvSqfBeneathImg(t *testing.T) {
	requireSquashfsTools(t)
	dir := t.TempDir()
	envSqf := filepath.Join(dir, "env.sqf")
	packRuntimeSqf(t, envSqf, envRuntime())
	other := filepath.Join(dir, "samtools.sqf")
	if err := os.WriteFile(other, []byte("not read as metadata in this test"), 0o644); err != nil {
		t.Fatalf("write: %v", err)
	}
	img := filepath.Join(dir, "env.img")

	got := orderOverlays([]string{envSqf, other, img})
	want := []string{other, envSqf, img}
	if !slices.Equal(got, want) {
		t.Fatalf("order = %v, want %v", got, want)
	}
}

// With no .img, orderOverlays leaves the list exactly as given.
func TestOrderOverlaysNoImg(t *testing.T) {
	overlays := []string{"/images/a.sqf", "/images/b.sqf"}
	got := orderOverlays(overlays)
	if !slices.Equal(got, overlays) {
		t.Fatalf("order = %v, want unchanged %v", got, overlays)
	}
}

// os overlays sort ahead of everything else, each group keeping its own
// relative order regardless of how they were interleaved in the request.
func TestOrderOverlaysOsFirst(t *testing.T) {
	requireSquashfsTools(t)
	dir := t.TempDir()

	osA := filepath.Join(dir, "os-a.sqf")
	packRuntimeSqf(t, osA, meta.Runtime{
		SchemaVersion: meta.SchemaVersion, Name: "ubuntu24", Type: catalog.TypeOS, Platform: meta.NativePlatform(),
	})
	osB := filepath.Join(dir, "os-b.sqf")
	packRuntimeSqf(t, osB, meta.Runtime{
		SchemaVersion: meta.SchemaVersion, Name: "cuda-libs", Type: catalog.TypeOS, Platform: meta.NativePlatform(),
	})
	appA := filepath.Join(dir, "samtools--1.22.sqf")
	packRuntimeSqf(t, appA, meta.Runtime{
		SchemaVersion: meta.SchemaVersion, Name: "samtools/1.22", Type: catalog.TypeApp,
		Platform: meta.NativePlatform(), Prefix: "/cnt/samtools/1.22",
	})
	appB := filepath.Join(dir, "bcftools--1.20.sqf")
	packRuntimeSqf(t, appB, meta.Runtime{
		SchemaVersion: meta.SchemaVersion, Name: "bcftools/1.20", Type: catalog.TypeApp,
		Platform: meta.NativePlatform(), Prefix: "/cnt/bcftools/1.20",
	})
	envSqf := filepath.Join(dir, "env.sqf")
	packRuntimeSqf(t, envSqf, envRuntime())
	img := filepath.Join(dir, "env.img")

	got := orderOverlays([]string{appA, osA, envSqf, img, appB, osB})
	want := []string{osA, osB, appA, appB, envSqf, img}
	if !slices.Equal(got, want) {
		t.Fatalf("order = %v, want %v", got, want)
	}
}

// autoloadSnapshot finds a personal snapshot beside a writable .img and adds
// it to the overlay list, with a diagnostic naming what was autoloaded.
func TestAutoloadSnapshotFindsPersonalLine(t *testing.T) {
	requireSquashfsTools(t)
	t.Setenv("USER", "alice")
	dir := t.TempDir()
	envSqf := filepath.Join(dir, "env-alice.sqf")
	packRuntimeSqf(t, envSqf, envRuntime())
	img := filepath.Join(dir, "env-alice.img")

	overlays, diagnostics := autoloadSnapshot([]string{img})
	if !slices.Contains(overlays, envSqf) {
		t.Fatalf("overlays = %v, want %s autoloaded", overlays, envSqf)
	}
	if len(diagnostics) != 1 {
		t.Fatalf("diagnostics = %v, want exactly one", diagnostics)
	}
}

// An overlay the caller already listed explicitly is not autoloaded a second
// time.
func TestAutoloadSnapshotSkipsWhenAlreadyListed(t *testing.T) {
	requireSquashfsTools(t)
	t.Setenv("USER", "alice")
	dir := t.TempDir()
	envSqf := filepath.Join(dir, "env-alice.sqf")
	packRuntimeSqf(t, envSqf, envRuntime())
	img := filepath.Join(dir, "env-alice.img")

	overlays, diagnostics := autoloadSnapshot([]string{envSqf, img})
	if len(overlays) != 2 {
		t.Fatalf("overlays = %v, want no duplicate append", overlays)
	}
	if len(diagnostics) != 0 {
		t.Fatalf("diagnostics = %v, want none for an already-listed overlay", diagnostics)
	}
}

// A blocked slot (something occupying the derived path that isn't a
// snapshot) is treated as nothing found: the .img mounts alone rather than
// autoload guessing.
func TestAutoloadSnapshotSkipsWhenBlocked(t *testing.T) {
	requireSquashfsTools(t)
	t.Setenv("USER", "alice")
	dir := t.TempDir()
	if err := os.WriteFile(filepath.Join(dir, "env-alice.sqf"), []byte("junk"), 0o644); err != nil {
		t.Fatalf("write: %v", err)
	}
	img := filepath.Join(dir, "env-alice.img")

	overlays, diagnostics := autoloadSnapshot([]string{img})
	if !slices.Equal(overlays, []string{img}) {
		t.Fatalf("overlays = %v, want unchanged", overlays)
	}
	if len(diagnostics) != 0 {
		t.Fatalf("diagnostics = %v, want none when blocked", diagnostics)
	}
}

// eligibility is selectRootWith's lookup over a fixed set, so root selection
// can be exercised without a real image to read metadata out of.
func eligibility(rootPaths ...string) func(string) bool {
	set := map[string]bool{}
	for _, p := range rootPaths {
		set[p] = true
	}
	return func(path string) bool { return set[path] }
}

// The first root-eligible overlay is pulled out; every other entry, root-
// eligible or not, keeps its place in declaration order.
func TestSelectRootPullsFirstEligibleOverlay(t *testing.T) {
	overlays := []string{"/images/samtools--1.22.sqf", "/images/ubuntu24--base.sqf", "/images/other-os.sqf"}
	root, rest := selectRootWith(overlays, eligibility("/images/ubuntu24--base.sqf", "/images/other-os.sqf"))
	if root != "/images/ubuntu24--base.sqf" {
		t.Errorf("root = %q, want the first eligible entry", root)
	}
	want := []string{"/images/samtools--1.22.sqf", "/images/other-os.sqf"}
	if !slices.Equal(rest, want) {
		t.Errorf("rest = %v, want %v", rest, want)
	}
}

// With no root-eligible entry, selectRoot reports none and leaves the list
// untouched — the caller falls back to the configured default.
func TestSelectRootNoneEligible(t *testing.T) {
	overlays := []string{"/images/samtools--1.22.sqf", "/images/bcftools--1.20.sqf"}
	root, rest := selectRootWith(overlays, eligibility())
	if root != "" {
		t.Errorf("root = %q, want none", root)
	}
	if !slices.Equal(rest, overlays) {
		t.Errorf("rest = %v, want unchanged %v", rest, overlays)
	}
}

// A :ro/:rw suffix on the overlay path must not defeat the eligibility
// lookup, which is keyed by the clean path.
func TestSelectRootStripsMountSuffix(t *testing.T) {
	root, rest := selectRootWith([]string{"/images/ubuntu24--base.sqf:ro"}, eligibility("/images/ubuntu24--base.sqf"))
	if root != "/images/ubuntu24--base.sqf:ro" {
		t.Errorf("root = %q, want the suffixed entry returned as given", root)
	}
	if len(rest) != 0 {
		t.Errorf("rest = %v, want empty", rest)
	}
}

// A .sif wins root unconditionally when present, regardless of where it
// falls among the requested overlays — its only valid use is root, unlike an
// os-typed overlay, which has a legitimate non-root use too.
func TestSelectRootSifWinsRegardlessOfPosition(t *testing.T) {
	requireSquashfsTools(t)
	dir := t.TempDir()

	osSqf := filepath.Join(dir, "cuda-libs.sqf")
	packRuntimeSqf(t, osSqf, meta.Runtime{
		SchemaVersion: meta.SchemaVersion, Name: "cuda-libs", Type: catalog.TypeOS, Platform: meta.NativePlatform(),
	})
	sif := filepath.Join(dir, "foreign.sif")
	if err := os.WriteFile(sif, []byte("not a real sif, just needs the extension"), 0o644); err != nil {
		t.Fatal(err)
	}

	for _, overlays := range [][]string{{sif, osSqf}, {osSqf, sif}} {
		root, rest := selectRoot(overlays)
		if root != sif {
			t.Errorf("selectRoot(%v) root = %q, want the sif %q", overlays, root, sif)
		}
		if !slices.Equal(rest, []string{osSqf}) {
			t.Errorf("selectRoot(%v) rest = %v, want %v", overlays, rest, []string{osSqf})
		}
	}
}

// Two .sif overlays cannot both be the root, so the request is refused
// outright rather than silently mounting the loser as a no-op --overlay.
func TestEnsureAtMostOneSifRefusesTwo(t *testing.T) {
	err := ensureAtMostOneSif([]string{"/images/a.sif", "/images/b.sif:ro"})
	if err == nil {
		t.Fatal("two .sif overlays were accepted")
	}
	for _, want := range []string{"/images/a.sif", "/images/b.sif"} {
		if !strings.Contains(err.Error(), want) {
			t.Errorf("error %q does not name %q", err, want)
		}
	}
}

// One .sif, or none at all, is the ordinary case.
func TestEnsureAtMostOneSifAllowsOne(t *testing.T) {
	if err := ensureAtMostOneSif([]string{"/images/samtools--1.22.sqf", "/images/a.sif"}); err != nil {
		t.Fatalf("one .sif was refused: %v", err)
	}
	if err := ensureAtMostOneSif([]string{"/images/samtools--1.22.sqf"}); err != nil {
		t.Fatalf("a list with no .sif was refused: %v", err)
	}
}

// isRootEligible reads real metadata: an os and a base overlay are both
// eligible, an app is not, and HasRequestedRoot answers the same question
// over a raw (unresolved) overlay list the way callers actually have it.
func TestIsRootEligibleReadsRealMetadata(t *testing.T) {
	requireSquashfsTools(t)
	dir := t.TempDir()

	osSqf := filepath.Join(dir, "os.sqf")
	packRuntimeSqf(t, osSqf, meta.Runtime{
		SchemaVersion: meta.SchemaVersion, Name: "ubuntu24", Type: catalog.TypeOS, Platform: meta.NativePlatform(),
	})
	appSqf := filepath.Join(dir, "samtools--1.22.sqf")
	packRuntimeSqf(t, appSqf, meta.Runtime{
		SchemaVersion: meta.SchemaVersion, Name: "samtools/1.22", Type: catalog.TypeApp,
		Platform: meta.NativePlatform(), Prefix: "/cnt/samtools/1.22",
	})

	if !isRootEligible(osSqf) {
		t.Error("an os-typed image is not root-eligible")
	}
	if isRootEligible(appSqf) {
		t.Error("an app-typed image is root-eligible")
	}

	if !HasRequestedRoot([]string{appSqf, osSqf}) {
		t.Error("HasRequestedRoot missed the os overlay in the list")
	}
	if HasRequestedRoot([]string{appSqf}) {
		t.Error("HasRequestedRoot found a root among app-only overlays")
	}
}

// A plain .sif carries no condatainer metadata at all, but it is still
// Apptainer's own native root format and must be root-eligible on that basis
// alone — this is the only way left to run one directly, now that -b/
// --base-image is gone.
func TestIsRootEligibleAcceptsAForeignSif(t *testing.T) {
	sif := filepath.Join(t.TempDir(), "foreign.sif")
	if err := os.WriteFile(sif, []byte("not a real sif, just needs the extension"), 0o644); err != nil {
		t.Fatal(err)
	}
	if !isRootEligible(sif) {
		t.Error("a .sif with no condatainer metadata is not root-eligible")
	}
}

// libexec/ is bound for nested running with no conda environment mounted, and
// only when asked.
func TestSetupBindsLibexecOnlyForNestedRunning(t *testing.T) {
	// A bound temp dir would cover a libexec made under t.TempDir() and the
	// bind would be deduplicated away.
	for _, env := range append(tmpDirEnvVars, "SLURM_TMPDIR", "PBS_TMPDIR", "LSF_TMPDIR") {
		t.Setenv(env, "")
	}
	dir := filepath.Join(t.TempDir(), "libexec")
	if err := os.MkdirAll(filepath.Join(dir, "bin"), 0o755); err != nil {
		t.Fatal(err)
	}
	if err := os.WriteFile(filepath.Join(dir, "bin", "micromamba"), []byte("#!/bin/sh\n"), 0o755); err != nil {
		t.Fatal(err)
	}
	prev := config.GlobalDataPaths
	config.GlobalDataPaths.LibexecDirs = []string{dir}
	t.Cleanup(func() { config.GlobalDataPaths = prev })
	real, err := filepath.EvalSymlinks(dir)
	if err != nil {
		t.Fatal(err)
	}

	bound := func(cfg SetupConfig) bool {
		res, err := Setup(cfg)
		if err != nil {
			t.Fatalf("Setup: %v", err)
		}
		for _, p := range res.BindPaths {
			if strings.Contains(p, real) {
				return true
			}
		}
		return false
	}
	if bound(SetupConfig{}) {
		t.Error("libexec bound with no environment and no nested running")
	}
	if !bound(SetupConfig{BindLibexec: true}) {
		t.Error("libexec not bound for nested running")
	}
}
