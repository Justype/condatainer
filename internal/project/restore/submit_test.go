package restore

import (
	"strings"
	"testing"

	"github.com/Justype/condatainer/catalog"
	"github.com/Justype/condatainer/internal/config"
)

// buildStep is one plan entry, named so a partition case reads as a graph.
func buildStep(artifact string, direct bool, dependsOn ...string) Step {
	return Step{Artifact: artifact, Name: artifact, Action: ActionBuild,
		Direct: direct, DependsOn: dependsOn}
}

// declares marks the steps whose recipes carry scheduler directives.
func declares(names ...string) map[string]bool {
	out := map[string]bool{}
	for _, name := range names {
		out[artifactKey(name)] = true
	}
	return out
}

// transientSet marks the steps that are build dependencies of this restore.
func transientSet(steps []Step, opts Options) map[string]bool {
	out := map[string]bool{}
	for _, step := range steps {
		out[stepKey(step)] = isTransient(step, opts)
	}
	return out
}

func TestPartitionSubmitsOnlyStepsWithDirectives(t *testing.T) {
	steps := []Step{buildStep("light", true), buildStep("heavy", true)}
	submit, deferred := partition(steps, declares("heavy"), transientSet(steps, Options{}))

	if submit[artifactKey("light")] {
		t.Error("a step with no directives was submitted")
	}
	if !submit[artifactKey("heavy")] {
		t.Error("a step declaring directives was not submitted")
	}
	if len(deferred) != 0 {
		t.Errorf("nothing is a build dependency here, deferred = %v", deferred)
	}
}

func TestPartitionSubmitsWhateverWaitsOnASubmittedStep(t *testing.T) {
	// tool declares directives; report pins nothing of its own but cannot run
	// until tool exists.
	steps := []Step{buildStep("tool", true), buildStep("report", true, "tool")}
	submit, _ := partition(steps, declares("tool"), transientSet(steps, Options{}))

	if !submit[artifactKey("report")] {
		t.Error("a step waiting on a submitted dependency was left to run locally")
	}
}

func TestPartitionCarriesABuildDependencysDirectivesUp(t *testing.T) {
	// index is a build dependency: nothing pins it, so it is produced inside
	// whichever job builds genome. Its directives are what that job must use.
	steps := []Step{buildStep("index", false), buildStep("genome", true, "index")}
	submit, deferred := partition(steps, declares("index"), transientSet(steps, Options{}))

	if !submit[artifactKey("genome")] {
		t.Error("a heavy build dependency did not put its dependent in a job")
	}
	if submit[artifactKey("index")] {
		t.Error("a build dependency was given its own job")
	}
	if !deferred[artifactKey("index")] {
		t.Error("a build dependency of a submitted step was not left to that job")
	}
}

func TestPartitionCarriesDirectivesThroughAChainOfBuildDependencies(t *testing.T) {
	// Two levels of build dependency. The directives are at the bottom, and only
	// reach the top by being carried through the middle.
	steps := []Step{
		buildStep("bottom", false),
		buildStep("middle", false, "bottom"),
		buildStep("top", true, "middle"),
	}
	submit, deferred := partition(steps, declares("bottom"), transientSet(steps, Options{}))

	if !submit[artifactKey("top")] {
		t.Error("directives two levels down did not reach the pin")
	}
	for _, inline := range []string{"bottom", "middle"} {
		if !deferred[artifactKey(inline)] {
			t.Errorf("%s was not left to the job that needs it", inline)
		}
	}
}

func TestPartitionKeepsABuildDependencyALocalDependentStillNeeds(t *testing.T) {
	// shared feeds two dependents: one goes to a job, one stays here. The local
	// one needs it now, so it is produced now and the job rebuilds its own copy.
	steps := []Step{
		buildStep("shared", false),
		buildStep("queued", true, "shared"),
		buildStep("here", true, "shared"),
	}
	submit, deferred := partition(steps, declares("queued"), transientSet(steps, Options{}))

	if !submit[artifactKey("queued")] || submit[artifactKey("here")] {
		t.Fatalf("unexpected submission set: %v", submit)
	}
	if deferred[artifactKey("shared")] {
		t.Error("a build dependency a local step still opens was skipped")
	}
}

func TestPartitionGivesAKeptBuildDependencyItsOwnJob(t *testing.T) {
	// --keep-build-deps installs it, so it is an ordinary step: its own job, its
	// own directives, and an afterok edge into its dependent.
	steps := []Step{buildStep("index", false), buildStep("genome", true, "index")}
	opts := Options{KeepBuildDeps: true}
	submit, deferred := partition(steps, declares("index"), transientSet(steps, opts))

	if !submit[artifactKey("index")] {
		t.Error("a kept build dependency did not get its own job")
	}
	if !submit[artifactKey("genome")] {
		t.Error("the dependent of a submitted step was not submitted")
	}
	if len(deferred) != 0 {
		t.Errorf("nothing is discarded under --keep-build-deps, deferred = %v", deferred)
	}
}

func TestPartitionTreatsTwoDestinationsAsTwoJobs(t *testing.T) {
	first := buildStep("tool", true)
	first.Destination = "overlays/a.sqf"
	second := buildStep("tool", true)
	second.Destination = "overlays/b.sqf"
	steps := []Step{first, second}

	needs := map[string]bool{stepKey(first): true, stepKey(second): true}
	submit, _ := partition(steps, needs, transientSet(steps, Options{}))

	if !submit[stepKey(first)] || !submit[stepKey(second)] {
		t.Fatalf("each declared destination is one output and one job: %v", submit)
	}
	if len(submit) != 2 {
		t.Errorf("two destinations collapsed into %d job(s)", len(submit))
	}
}

func TestWaitForNamesOnlySubmittedDependencies(t *testing.T) {
	j := &jobs{ids: map[string]string{"tool": "4001"}}
	step := buildStep("report", true, "installed", "tool")

	ids := j.waitFor(step)
	if len(ids) != 1 || ids[0] != "4001" {
		t.Fatalf("waitFor = %v, want just the submitted dependency's job", ids)
	}
}

func TestRestoreCommandCarriesThePolicyThatDecidedThePlan(t *testing.T) {
	step := Step{Artifact: "provenance/samtools-1.23.1-sha256-abc"}
	opts := Options{Match: MatchIdentity, SkipPrebuilt: true, KeepBuildDeps: true}

	cmd := restoreCommand("/projects/rna", step, opts)
	for _, want := range []string{
		"--project '/projects/rna'",
		"--only 'provenance/samtools-1.23.1-sha256-abc'",
		"--match identity", "--no-prebuilt", "--keep-build-deps",
	} {
		if !strings.Contains(cmd, want) {
			t.Errorf("command %q is missing %q", cmd, want)
		}
	}
}

func TestRestoreCommandOmitsUnsetPolicy(t *testing.T) {
	cmd := restoreCommand("/projects/rna", Step{Artifact: "provenance/x"}, Options{})
	for _, unwanted := range []string{"--match", "--no-prebuilt", "--keep-build-deps"} {
		if strings.Contains(cmd, unwanted) {
			t.Errorf("command %q carries %q, which was not asked for", cmd, unwanted)
		}
	}
}

func TestShellQuoteSurvivesAQuoteInAPath(t *testing.T) {
	got := shellQuote("/home/o'brien/rna project")
	if want := `'/home/o'\''brien/rna project'`; got != want {
		t.Fatalf("shellQuote = %s, want %s", got, want)
	}
}

// build.always_submit_data calls for a job for a data rebuild and for nothing
// else; a recipe's own directives still do, whatever its type.
func TestCallsForJobFollowsAlwaysSubmitData(t *testing.T) {
	prev := config.Global.Build.AlwaysSubmitData
	config.Global.Build.AlwaysSubmitData = true
	t.Cleanup(func() { config.Global.Build.AlwaysSubmitData = prev })

	if !callsForJob(catalog.TypeData, nil) {
		t.Error("a data rebuild does not call for a job under always_submit_data")
	}
	if callsForJob(catalog.TypeApp, nil) {
		t.Error("an app rebuild calls for a job under always_submit_data")
	}
	config.Global.Build.AlwaysSubmitData = false
	if callsForJob(catalog.TypeData, nil) {
		t.Error("a data rebuild calls for a job with the flag off")
	}
}
