package restore

import (
	"context"
	"errors"
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/Justype/condatainer/internal/artifact/capsule"
	"github.com/Justype/condatainer/internal/artifact/compare"
	"github.com/Justype/condatainer/internal/artifact/meta"
	"github.com/Justype/condatainer/internal/build"
	"github.com/Justype/condatainer/internal/conda"
	"github.com/Justype/condatainer/internal/project/lock"
)

// A project whose artifacts are all present acquires nothing and reports each
// one as adopted. This is the ordinary re-run, and it must touch no build.
func TestRunAdoptsWithoutBuilding(t *testing.T) {
	root := projectRoot(t)
	app := artifact(t, "star/2.7.11b", "echo star\n")
	appPath := vendor(t, root, app, "echo star\n")

	l := lock.New()
	l.Selections["star/2.7.11b"] = lock.Selection{Artifact: appPath}

	report, err := Run(context.Background(), root, l, Options{lookup: installed("star/2.7.11b")})
	if err != nil {
		t.Fatal(err)
	}
	if !report.Complete() {
		t.Fatalf("report = %#v", report)
	}
	if len(report.Results) != 1 {
		t.Fatalf("results = %#v, want one", report.Results)
	}
	result := report.Results[0]
	if result.Outcome != OutcomeAdopted || result.Verdict != compare.Exact {
		t.Errorf("result = %#v, want an exact adoption", result)
	}
	if result.Path == "" {
		t.Error("an adopted result has no path")
	}
}

// An equivalent adoption is reported as a substitution, with both identities.
func TestRunReportsAnEquivalentAdoption(t *testing.T) {
	root := projectRoot(t)
	app := artifact(t, "star/2.7.11b", "echo star\n")
	appPath := vendor(t, root, app, "echo star\n")

	l := lock.New()
	l.Selections["star/2.7.11b"] = lock.Selection{Artifact: appPath}

	report, err := Run(context.Background(), root, l, Options{lookup: substitute("star/2.7.11b")})
	if err != nil {
		t.Fatal(err)
	}
	result := report.Results[0]
	if result.Verdict != compare.Equivalent {
		t.Errorf("verdict = %q, want equivalent", result.Verdict)
	}
	if result.Identity == "" || result.Found == "" || result.Identity == result.Found {
		t.Errorf("result = %#v, want both identities reported", result)
	}
}

// Planning failures stop a restore before it acquires anything.
func TestRunStopsOnAPlanningProblem(t *testing.T) {
	root := projectRoot(t)
	l := lock.New()
	l.Selections["star/2.7.11b"] = lock.Selection{Artifact: "provenance/missing--1.0@abc123456789"}

	report, err := Run(context.Background(), root, l, Options{lookup: nothingInstalled})
	if !errors.Is(err, ErrIncomplete) {
		t.Fatalf("err = %v, want ErrIncomplete", err)
	}
	if len(report.Problems) == 0 {
		t.Error("no problem was reported")
	}
	if len(report.Results) != 0 || len(report.Failures) != 0 {
		t.Errorf("a rejected plan produced work: %#v", report)
	}
}

// A recorded remote is not silently rebuilt: the fetch is an acquisition this
// build cannot perform, and --no-prebuilt is how a caller asks to build from source.
func TestRunFailsAFetchRatherThanRebuildingSilently(t *testing.T) {
	root := projectRoot(t)
	app := artifact(t, "star/2.7.11b", "echo star\n")
	appPath := vendor(t, root, app, "echo star\n")

	l := lock.New()
	l.Selections["star/2.7.11b"] = lock.Selection{Artifact: appPath}
	digest := "sha256:" + strings.Repeat("a", 64)
	if err := l.AddRemote(appPath, lock.Remote{Repository: "ghcr.io/x/star", ManifestDigest: digest}); err != nil {
		t.Fatal(err)
	}

	report, err := Run(context.Background(), root, l, Options{lookup: nothingInstalled})
	if !errors.Is(err, ErrIncomplete) {
		t.Fatalf("err = %v, want ErrIncomplete", err)
	}
	if len(report.Failures) != 1 {
		t.Fatalf("failures = %#v, want one", report.Failures)
	}
	if !strings.Contains(report.Failures[0].Reason, "--no-prebuilt") {
		t.Errorf("failure does not name the remedy: %q", report.Failures[0].Reason)
	}
}

// A cancelled context stops the loop rather than running the next build.
func TestRunHonoursCancellation(t *testing.T) {
	root := projectRoot(t)
	app := artifact(t, "star/2.7.11b", "echo star\n")
	appPath := vendor(t, root, app, "echo star\n")

	l := lock.New()
	l.Selections["star/2.7.11b"] = lock.Selection{Artifact: appPath}

	ctx, cancel := context.WithCancel(context.Background())
	cancel()
	if _, err := Run(ctx, root, l, Options{lookup: installed("star/2.7.11b")}); !errors.Is(err, context.Canceled) {
		t.Fatalf("err = %v, want context.Canceled", err)
	}
}

// Edge order is the contract: a rebuild re-derives its edges from what it
// mounts, so the paths must arrive in the manifest's own order.
func TestDependencyPathsFollowManifestOrder(t *testing.T) {
	first := artifact(t, "aaa/1.0", "echo a\n")
	second := artifact(t, "zzz/1.0", "echo z\n")
	entry := &lock.Entry{Manifest: meta.Manifest{
		Name:         "index/1.0",
		Dependencies: []meta.Dependency{edge(second), edge(first)},
	}}
	available := map[string]string{
		lock.EntryPath(entryName(second)): "/images/zzz.sqf",
		lock.EntryPath(entryName(first)):  "/images/aaa.sqf",
	}

	deps, err := dependencyPaths(entry, available)
	if err != nil {
		t.Fatal(err)
	}
	if len(deps) != 2 || deps[0].Name != "zzz/1.0" || deps[1].Name != "aaa/1.0" {
		t.Fatalf("deps = %#v, want the manifest's order, not sorted", deps)
	}
	if deps[0].Path != "/images/zzz.sqf" || deps[1].Path != "/images/aaa.sqf" {
		t.Errorf("paths did not follow their edges: %#v", deps)
	}
}

// Nothing is resolved by name: a dependency that was not restored first is an
// error, never a lookup.
func TestDependencyPathsRefuseAnUnrestoredDependency(t *testing.T) {
	dep := artifact(t, "zlib/1.3", "echo zlib\n")
	entry := &lock.Entry{Manifest: meta.Manifest{
		Name: "index/1.0", Dependencies: []meta.Dependency{edge(dep)},
	}}
	if _, err := dependencyPaths(entry, map[string]string{}); err == nil {
		t.Fatal("a missing dependency path was accepted")
	}
}

func TestPlaceAtPublishesByRename(t *testing.T) {
	dir := t.TempDir()
	output := filepath.Join(dir, "built.sqf")
	if err := os.WriteFile(output, []byte("payload"), 0o644); err != nil {
		t.Fatal(err)
	}
	destination := filepath.Join(dir, "overlays", "tool.sqf")

	got, err := placeAt(destination, output)
	if err != nil {
		t.Fatal(err)
	}
	if got != destination {
		t.Fatalf("placeAt = %q, want %q", got, destination)
	}
	if data, err := os.ReadFile(destination); err != nil || string(data) != "payload" {
		t.Errorf("destination = %q, %v", data, err)
	}
	if _, err := os.Stat(output); !os.IsNotExist(err) {
		t.Error("the staged file survived the rename")
	}
}

// A symlink points somewhere the project does not describe, so replacing what
// it targets would write outside the checkout.
func TestPlaceAtRefusesASymlink(t *testing.T) {
	dir := t.TempDir()
	output := filepath.Join(dir, "built.sqf")
	if err := os.WriteFile(output, []byte("payload"), 0o644); err != nil {
		t.Fatal(err)
	}
	outside := filepath.Join(t.TempDir(), "real.sqf")
	if err := os.WriteFile(outside, []byte("someone else's"), 0o644); err != nil {
		t.Fatal(err)
	}
	destination := filepath.Join(dir, "tool.sqf")
	if err := os.Symlink(outside, destination); err != nil {
		t.Skipf("symlinks unavailable: %v", err)
	}

	if _, err := placeAt(destination, output); err == nil {
		t.Fatal("a symlinked destination was overwritten")
	}
	if data, _ := os.ReadFile(outside); string(data) != "someone else's" {
		t.Error("the symlink target was modified")
	}
}

func TestPlaceAtRefusesANonRegularDestination(t *testing.T) {
	dir := t.TempDir()
	output := filepath.Join(dir, "built.sqf")
	if err := os.WriteFile(output, []byte("payload"), 0o644); err != nil {
		t.Fatal(err)
	}
	destination := filepath.Join(dir, "tool.sqf")
	if err := os.MkdirAll(destination, 0o775); err != nil {
		t.Fatal(err)
	}
	if _, err := placeAt(destination, output); err == nil {
		t.Fatal("a directory was accepted as a destination")
	}
}

// #INPUT: answers make an interactive rebuild plannable; without them the plan
// refuses, because there is no terminal inside a restore.
func TestAnswersMakeAnInteractiveRebuildPlannable(t *testing.T) {
	root := projectRoot(t)
	recipe := "#INPUT: which release?\nread -r a\n"
	app := artifact(t, "star/2.7.11b", recipe)
	app.Source.RequiresInput = true
	appPath := vendor(t, root, app, recipe)

	l := lock.New()
	l.Selections["star/2.7.11b"] = lock.Selection{Artifact: appPath}

	if plan := Compute(root, l, Options{lookup: nothingInstalled}); plan.Complete() {
		t.Fatal("an interactive rebuild was planned with no answers")
	}
	plan := Compute(root, l, Options{lookup: nothingInstalled,
		Answers: map[string][]string{appPath: {"49"}}})
	if !plan.Complete() {
		t.Fatalf("problems: %v", plan.Problems)
	}
}

// The locked build entry point is reached with the vendored bytes and the
// supplied paths, so a restore cannot quietly consult the catalog.
func TestRebuildUsesTheVendoredSources(t *testing.T) {
	root := projectRoot(t)
	recipe := "#TYPE: app\necho star\n"
	app := artifact(t, "star/2.7.11b", recipe)
	appPath := vendor(t, root, app, recipe)

	verified, problems := lock.Verify(root, mustLock(appPath))
	if len(problems) > 0 {
		t.Fatalf("problems: %v", problems)
	}
	sources, err := lock.Sources(root, verified.Entries[appPath])
	if err != nil {
		t.Fatal(err)
	}
	if got := string(sources[meta.RecipeFileName]); got != recipe {
		t.Fatalf("vendored recipe = %q, want it verbatim", got)
	}

	object, err := build.NewLockedObject(context.Background(), build.LockedSpec{
		Manifest: verified.Entries[appPath].Manifest,
		Sources:  sources,
		Output:   filepath.Join(t.TempDir(), "out.sqf"),
		TmpRoot:  t.TempDir(),
	})
	if err != nil {
		t.Fatal(err)
	}
	if object.NameVersion() != "star/2.7.11b" {
		t.Errorf("name = %q", object.NameVersion())
	}
}

func mustLock(artifactPath string) *lock.Lock {
	l := lock.New()
	l.Selections["star/2.7.11b"] = lock.Selection{Artifact: artifactPath}
	return l
}

func entryName(manifest meta.Manifest) string {
	return capsule.EntryName(manifest.Name, manifest.Keys.Identity.Digest())
}

// A closure node no selection names is scaffolding: it exists so its dependent
// can be built and is discarded with the restore.
func TestIsTransientFollowsWhatAskedForTheArtifact(t *testing.T) {
	for _, tc := range []struct {
		name string
		step Step
		opts Options
		want bool
	}{
		{"BuildInput", Step{}, Options{}, true},
		{"SelectedByName", Step{Direct: true}, Options{}, false},
		{"SelectedAtAPath", Step{Direct: true, Destination: "overlays/x.sqf"}, Options{}, false},
		{"KeptOnRequest", Step{}, Options{KeepBuildDeps: true}, false},
	} {
		t.Run(tc.name, func(t *testing.T) {
			if got := isTransient(tc.step, tc.opts); got != tc.want {
				t.Fatalf("isTransient = %v, want %v", got, tc.want)
			}
		})
	}
}

// Blockage names the artifact that actually failed and propagates down a chain,
// so one broken dependency reads as one cause rather than one per artifact.
func TestBlockedByNamesTheRootCause(t *testing.T) {
	stopped := map[string]string{"a": "a"}
	step := Step{Artifact: "b", Action: ActionBuild, DependsOn: []string{"a"}}
	if got := blockedBy(step, stopped); got != "a" {
		t.Fatalf("blockedBy = %q, want the failed artifact", got)
	}

	// b is now blocked by a; c depends on b and must still report a.
	stopped["b"] = "a"
	deeper := Step{Artifact: "c", Action: ActionBuild, DependsOn: []string{"b"}}
	if got := blockedBy(deeper, stopped); got != "a" {
		t.Fatalf("blockedBy = %q, want the root cause rather than the nearest edge", got)
	}
}

// An adopted or fetched artifact opens none of its inputs, so a dependency that
// could not be produced does not stop it.
func TestBlockedByIgnoresAnArtifactThatIsNotBuilt(t *testing.T) {
	stopped := map[string]string{"a": "a"}
	for _, action := range []Action{ActionAdopt, ActionFetch} {
		step := Step{Artifact: "b", Action: action, DependsOn: []string{"a"}}
		if got := blockedBy(step, stopped); got != "" {
			t.Errorf("blockedBy(%s) = %q, want no blockage", action, got)
		}
	}
}

// The transient root is created once and removed whole, so an input two
// dependents need is produced in one place and nothing survives the restore.
func TestTransientRootIsCreatedOnceAndRemoved(t *testing.T) {
	t.Setenv("XDG_DATA_HOME", t.TempDir())
	inputs := &transientRoot{}

	first, err := inputs.dir()
	if err != nil {
		t.Fatal(err)
	}
	second, err := inputs.dir()
	if err != nil {
		t.Fatal(err)
	}
	if first != second {
		t.Fatalf("dir() = %q then %q, want one directory for the whole restore", first, second)
	}
	if _, err := os.Stat(first); err != nil {
		t.Fatalf("the transient root was not created: %v", err)
	}

	inputs.remove()
	if _, err := os.Stat(first); !os.IsNotExist(err) {
		t.Fatalf("the transient root survived the restore: %v", err)
	}
}

// remove is safe on a root nothing ever needed, which is the common case: most
// restores adopt everything and never produce a build input.
func TestTransientRootRemoveWithoutUse(t *testing.T) {
	(&transientRoot{}).remove()
}

// A dependency that cannot be acquired stops its dependent, and the dependent is
// reported blocked rather than failed: one cause, named once.
func TestRunBlocksADependentRatherThanFailingIt(t *testing.T) {
	root := projectRoot(t)
	app := artifact(t, "star/2.7.11b", "echo star\n")
	appPath := vendor(t, root, app, "echo star\n")
	index := artifact(t, "grch38/star-index", "echo index\n", edge(app))
	indexPath := vendor(t, root, index, "echo index\n")

	l := lock.New()
	l.Selections["grch38/star-index"] = lock.Selection{Artifact: indexPath}
	// The dependency is a fetch this build cannot perform, which is the one
	// acquisition failure reachable without running a real build.
	digest := "sha256:" + strings.Repeat("a", 64)
	if err := l.AddRemote(appPath, lock.Remote{Repository: "ghcr.io/x/star", ManifestDigest: digest}); err != nil {
		t.Fatal(err)
	}

	report, err := Run(context.Background(), root, l, Options{lookup: nothingInstalled})
	if !errors.Is(err, ErrIncomplete) {
		t.Fatalf("err = %v, want ErrIncomplete", err)
	}
	if len(report.Failures) != 1 || report.Failures[0].Name != "star/2.7.11b" {
		t.Fatalf("failures = %#v, want only the dependency", report.Failures)
	}
	if len(report.Blocked) != 1 {
		t.Fatalf("blocked = %#v, want the dependent", report.Blocked)
	}
	if got := report.Blocked[0]; got.Name != "grch38/star-index" || got.Cause != appPath {
		t.Errorf("blocked = %#v, want grch38/star-index blocked by %s", got, appPath)
	}
	if report.Complete() {
		t.Error("a report with a blocked artifact reports itself complete")
	}
}

// explicit.txt reproduces the identity and is always tried first. Its URLs rot,
// so under the default mode environment.yml is a second chance — it is the
// equivalence preimage, so a result matching it satisfies the mode.
func TestCondaSourcesTierByMatchMode(t *testing.T) {
	condaEntry := func(files ...string) *lock.Entry {
		return &lock.Entry{Manifest: meta.Manifest{Name: "cutadapt/5.0", BuildType: "conda"}}
	}
	both := map[string][]byte{
		conda.ExplicitFileName:    []byte("@EXPLICIT\n"),
		conda.EnvironmentFileName: []byte("dependencies: []\n"),
	}
	explicitOnly := map[string][]byte{conda.ExplicitFileName: []byte("@EXPLICIT\n")}

	for _, tc := range []struct {
		name    string
		sources map[string][]byte
		match   Match
		want    []string
	}{
		{"EquivalentTriesBoth", both, MatchEquivalent,
			[]string{conda.ExplicitFileName, conda.EnvironmentFileName}},
		{"IdentityReplaysExplicitOnly", both, MatchIdentity,
			[]string{conda.ExplicitFileName}},
		{"NothingToFallBackOn", explicitOnly, MatchEquivalent,
			[]string{conda.ExplicitFileName}},
	} {
		t.Run(tc.name, func(t *testing.T) {
			got := condaSources(condaEntry(), tc.sources, tc.match)
			if len(got) != len(tc.want) {
				t.Fatalf("sources = %v, want %v", got, tc.want)
			}
			for i := range got {
				if got[i] != tc.want[i] {
					t.Fatalf("sources = %v, want %v", got, tc.want)
				}
			}
		})
	}
}

// A script or definition rebuild has one source and no tiering: its recipe is
// vendored verbatim and reproduces the identity outright.
func TestCondaSourcesLeaveOtherBuildTypesAlone(t *testing.T) {
	entry := &lock.Entry{Manifest: meta.Manifest{Name: "star/2.7.11b", BuildType: "script"}}
	got := condaSources(entry, map[string][]byte{conda.EnvironmentFileName: nil}, MatchEquivalent)
	if len(got) != 1 || got[0] != "" {
		t.Fatalf("sources = %v, want one unnamed source", got)
	}
}

// The transient tree is what someone reads when a build fails, so each staging
// directory is named for the artifact in it rather than randomly.
func TestTransientStagingIsNamedForTheArtifact(t *testing.T) {
	root := projectRoot(t)
	dep := artifact(t, "samtools/1.21", "echo samtools\n")
	depPath := vendor(t, root, dep, "echo samtools\n")

	entry, problems := func() (*lock.Entry, []lock.Problem) {
		l := lock.New()
		l.Selections["samtools/1.21"] = lock.Selection{Artifact: depPath}
		verified, problems := lock.Verify(root, l)
		return verified.Entries[depPath], problems
	}()
	if len(problems) > 0 {
		t.Fatalf("problems = %v", problems)
	}

	want := capsule.EntryName(entry.Manifest.Name, entry.Identity.Digest())
	if !strings.Contains(want, "samtools--1.21@") {
		t.Fatalf("entry name = %q, want the encoded name and a short digest", want)
	}
}
