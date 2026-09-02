package restore

import (
	"encoding/json"
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/Justype/condatainer/catalog"
	"github.com/Justype/condatainer/internal/artifact/capsule"
	"github.com/Justype/condatainer/internal/artifact/key"
	"github.com/Justype/condatainer/internal/artifact/meta"
	"github.com/Justype/condatainer/internal/project/lock"
	"github.com/Justype/condatainer/internal/store"
)

func projectRoot(t *testing.T) string {
	t.Helper()
	root := t.TempDir()
	if err := os.MkdirAll(filepath.Join(root, lock.DirName), 0o775); err != nil {
		t.Fatal(err)
	}
	return root
}

// artifact builds a manifest whose keys its recipe really derives, so the plan
// runs over a closure that verifies rather than a stub.
func artifact(t *testing.T, name, recipe string, deps ...meta.Dependency) meta.Manifest {
	t.Helper()
	typ := catalog.TypeApp
	if len(deps) > 0 {
		typ = catalog.TypeData
	}
	manifest := meta.Manifest{
		SchemaVersion: meta.SchemaVersion,
		Name:          name,
		Type:          typ,
		BuildType:     "script",
		Platform:      meta.Platform{OS: "linux", Arch: meta.NativeArch()},
		Source:        meta.Source{Files: []string{meta.RecipeFileName}},
		Dependencies:  deps,
	}
	derived, err := key.Generate(manifest, key.Sources{meta.RecipeFileName: []byte(recipe)})
	if err != nil {
		t.Fatalf("generate keys for %s: %v", name, err)
	}
	manifest.Keys = derived.Keys()
	return manifest
}

func vendor(t *testing.T, root string, manifest meta.Manifest, recipe string) string {
	t.Helper()
	data, err := json.MarshalIndent(manifest, "", "  ")
	if err != nil {
		t.Fatal(err)
	}
	relative, err := lock.StageEntry(root, capsule.EntryName(manifest.Name, manifest.Keys.Identity.Digest()),
		map[string][]byte{meta.FileName: data, meta.RecipeFileName: []byte(recipe)})
	if err != nil {
		t.Fatal(err)
	}
	return relative
}

// A frozen environment records one key under snapshot-env-v1 and vendors no
// source: it was captured, not built.
func frozenEnv(t *testing.T, root string) string {
	t.Helper()
	ref := meta.KeyRef{Scheme: string(key.SnapshotEnvV1), SHA256: strings.Repeat("c", 64)}
	manifest := meta.Manifest{
		SchemaVersion: meta.SchemaVersion,
		Name:          meta.EnvName,
		Type:          catalog.TypeEnv,
		BuildType:     meta.BuildTypeSnapshot,
		Platform:      meta.Platform{OS: "linux", Arch: meta.NativeArch()},
		Keys:          meta.Keys{Identity: ref, Equiv: ref},
	}
	data, err := json.MarshalIndent(manifest, "", "  ")
	if err != nil {
		t.Fatal(err)
	}
	relative, err := lock.StageEntry(root, capsule.EntryName(manifest.Name, ref.Digest()),
		map[string][]byte{meta.FileName: data})
	if err != nil {
		t.Fatal(err)
	}
	return relative
}

// A frozen environment is refused while planning rather than planned as a
// rebuild: nothing it could be rebuilt from exists, so a registry copy is the
// only thing that produces it.
func TestComputeRefusesToRebuildAFrozenEnvironment(t *testing.T) {
	root := projectRoot(t)
	artifact := frozenEnv(t, root)

	l := lock.New()
	l.Pins[lock.PathPrefix+"overlays/env.sqf"] = lock.PinEntry{Artifact: artifact}

	plan := Compute(root, l, Options{lookup: nothingInstalled})
	if plan.Complete() {
		t.Fatal("a rebuild was planned for an artifact with no sources")
	}
	if len(plan.Steps) != 1 || plan.Steps[0].Action != ActionUnavailable {
		t.Fatalf("steps = %#v, want one unavailable step", plan.Steps)
	}
	if !strings.Contains(strings.Join(plan.Problems, " "), "project push") {
		t.Errorf("the refusal does not name the remedy: %v", plan.Problems)
	}

	// A recorded remote is the one thing that makes it obtainable — including
	// under --no-prebuilt, which chooses building over fetching and so has
	// nothing to say about something that cannot be built.
	l.Remotes = map[string][]lock.Remote{artifact: {{Repository: "ghcr.io/x/y", ManifestDigest: "sha256:" + strings.Repeat("d", 64)}}}
	for _, skip := range []bool{false, true} {
		plan := Compute(root, l, Options{lookup: nothingInstalled, SkipPrebuilt: skip})
		if !plan.Complete() || plan.Steps[0].Action != ActionFetch {
			t.Fatalf("SkipPrebuilt=%v: steps = %#v, problems = %v", skip, plan.Steps, plan.Problems)
		}
	}
}

func edge(manifest meta.Manifest) meta.Dependency {
	return meta.Dependency{Name: manifest.Name, Type: manifest.Type,
		Identity: manifest.Keys.Identity, Equiv: manifest.Keys.Equiv, Role: meta.RoleApp}
}

// nothingInstalled is a host with no images at all.
func nothingInstalled(string, meta.Keys, Match, []string) (store.Candidate, bool) {
	return store.Candidate{}, false
}

// installed reports an exact hit for the named artifacts.
func installed(names ...string) lookupFunc {
	have := make(map[string]bool, len(names))
	for _, name := range names {
		have[name] = true
	}
	return func(name string, keys meta.Keys, _ Match, _ []string) (store.Candidate, bool) {
		if !have[name] {
			return store.Candidate{}, false
		}
		return store.Candidate{Name: name, Path: "/images/" + name + ".sqf",
			Layout: store.LayoutFlat, Identity: keys.Identity, Equiv: keys.Equiv}, true
	}
}

// substitute reports a hit for the named artifact at a different identity but
// the locked equivalence — what a Conda re-solve onto the same versions looks
// like. Under MatchIdentity it is not a hit at all.
func substitute(name string) lookupFunc {
	return func(got string, keys meta.Keys, match Match, _ []string) (store.Candidate, bool) {
		if got != name || match.Normalize() == MatchIdentity {
			return store.Candidate{}, false
		}
		other := meta.KeyRef{Scheme: keys.Identity.Scheme, SHA256: strings.Repeat("e", 64)}
		return store.Candidate{Name: got, Path: "/images/" + got + ".sqf",
			Layout: store.LayoutFlat, Identity: other, Equiv: keys.Equiv}, true
	}
}

func steps(plan *Plan) map[string]Step {
	out := make(map[string]Step, len(plan.Steps))
	for _, step := range plan.Steps {
		out[step.Name] = step
	}
	return out
}

func TestComputeAdoptsWhatIsAlreadyPresent(t *testing.T) {
	root := projectRoot(t)
	app := artifact(t, "star/2.7.11b", "echo star\n")
	appPath := vendor(t, root, app, "echo star\n")

	l := lock.New()
	l.Pins["star/2.7.11b"] = lock.PinEntry{Artifact: appPath}

	plan := Compute(root, l, Options{lookup: installed("star/2.7.11b")})
	if !plan.Complete() {
		t.Fatalf("problems: %v", plan.Problems)
	}
	if plan.Work() != 0 {
		t.Errorf("a satisfied project reported %d step(s) of work", plan.Work())
	}
	step := steps(plan)["star/2.7.11b"]
	if step.Action != ActionAdopt || step.Path == "" || !step.Direct {
		t.Fatalf("step = %#v", step)
	}
	if len(step.Requests) != 1 || step.Requests[0] != "star/2.7.11b" {
		t.Errorf("requests = %v", step.Requests)
	}
}

// A recorded prebuilt is preferred to a build, and --no-prebuilt suppresses it.
func TestComputeChoosesFetchOrBuild(t *testing.T) {
	root := projectRoot(t)
	app := artifact(t, "star/2.7.11b", "echo star\n")
	appPath := vendor(t, root, app, "echo star\n")

	l := lock.New()
	l.Pins["star/2.7.11b"] = lock.PinEntry{Artifact: appPath}
	digest := "sha256:" + strings.Repeat("a", 64)
	if err := l.AddRemote(appPath, lock.Remote{Repository: "ghcr.io/x/star", ManifestDigest: digest}); err != nil {
		t.Fatal(err)
	}

	online := Compute(root, l, Options{lookup: nothingInstalled})
	if got := steps(online)["star/2.7.11b"].Action; got != ActionFetch {
		t.Errorf("action = %q, want fetch", got)
	}
	fromSource := Compute(root, l, Options{SkipPrebuilt: true, lookup: nothingInstalled})
	if got := steps(fromSource)["star/2.7.11b"].Action; got != ActionBuild {
		t.Errorf("--no-prebuilt action = %q, want build", got)
	}
	if !online.Complete() || !fromSource.Complete() {
		t.Errorf("problems: %v %v", online.Problems, fromSource.Problems)
	}
}

// Without an remote there is nowhere to fetch from, so a miss is a rebuild.
func TestComputeBuildsWhenNoOriginIsRecorded(t *testing.T) {
	root := projectRoot(t)
	app := artifact(t, "star/2.7.11b", "echo star\n")
	appPath := vendor(t, root, app, "echo star\n")

	l := lock.New()
	l.Pins["star/2.7.11b"] = lock.PinEntry{Artifact: appPath}

	if got := steps(Compute(root, l, Options{lookup: nothingInstalled}))["star/2.7.11b"].Action; got != ActionBuild {
		t.Fatalf("action = %q, want build", got)
	}
}

// Dependencies come first, and a closure-only node is planned but not direct.
func TestComputeOrdersDependenciesFirst(t *testing.T) {
	root := projectRoot(t)
	dep := artifact(t, "zlib/1.3", "echo zlib\n")
	depPath := vendor(t, root, dep, "echo zlib\n")
	data := artifact(t, "index/1.0", "echo index\n", edge(dep))
	dataPath := vendor(t, root, data, "echo index\n")

	l := lock.New()
	l.Pins["index/1.0"] = lock.PinEntry{Artifact: dataPath}

	plan := Compute(root, l, Options{lookup: nothingInstalled})
	if !plan.Complete() {
		t.Fatalf("problems: %v", plan.Problems)
	}
	if len(plan.Steps) != 2 {
		t.Fatalf("steps = %#v, want the closure planned", plan.Steps)
	}
	if plan.Steps[0].Name != "zlib/1.3" || plan.Steps[1].Name != "index/1.0" {
		t.Fatalf("order = %s then %s, want the dependency first", plan.Steps[0].Name, plan.Steps[1].Name)
	}
	if plan.Steps[0].Direct {
		t.Error("a closure-only node was reported as direct")
	}
	if got := plan.Steps[1].DependsOn; len(got) != 1 || got[0] != depPath {
		t.Errorf("depends_on = %v, want %q", got, depPath)
	}
}

// A present dependency is adopted even when its dependent must be rebuilt.
func TestComputeMixesAdoptionAndRebuild(t *testing.T) {
	root := projectRoot(t)
	dep := artifact(t, "zlib/1.3", "echo zlib\n")
	vendor(t, root, dep, "echo zlib\n")
	data := artifact(t, "index/1.0", "echo index\n", edge(dep))
	dataPath := vendor(t, root, data, "echo index\n")

	l := lock.New()
	l.Pins["index/1.0"] = lock.PinEntry{Artifact: dataPath}

	plan := Compute(root, l, Options{lookup: installed("zlib/1.3")})
	byName := steps(plan)
	if byName["zlib/1.3"].Action != ActionAdopt {
		t.Errorf("dependency action = %q, want adopt", byName["zlib/1.3"].Action)
	}
	if byName["index/1.0"].Action != ActionBuild {
		t.Errorf("root action = %q, want build", byName["index/1.0"].Action)
	}
	if plan.Work() != 1 {
		t.Errorf("work = %d, want only the rebuild", plan.Work())
	}
}

// A locked artifact is arch-bound unless its recipe declared noarch.
func TestComputeRejectsAForeignArchitecture(t *testing.T) {
	root := projectRoot(t)
	app := artifact(t, "star/2.7.11b", "echo star\n")
	appPath := vendor(t, root, app, "echo star\n")

	l := lock.New()
	l.Pins["star/2.7.11b"] = lock.PinEntry{Artifact: appPath}

	plan := Compute(root, l, Options{HostArch: "s390x", lookup: nothingInstalled})
	if plan.Complete() {
		t.Fatal("planned an acquisition for a foreign architecture")
	}
	if !strings.Contains(strings.Join(plan.Problems, "\n"), "this host is s390x") {
		t.Fatalf("problems = %v", plan.Problems)
	}
}

// noarch runs anywhere, so a foreign host is fine.
func TestComputeAllowsNoarchAnywhere(t *testing.T) {
	root := projectRoot(t)
	app := artifact(t, "refdata/1.0", "echo data\n")
	app.Platform.Arch = meta.ArchNone
	derived, err := key.Generate(app, key.Sources{meta.RecipeFileName: []byte("echo data\n")})
	if err != nil {
		t.Fatal(err)
	}
	app.Keys = derived.Keys()
	appPath := vendor(t, root, app, "echo data\n")

	l := lock.New()
	l.Pins["refdata/1.0"] = lock.PinEntry{Artifact: appPath}

	plan := Compute(root, l, Options{HostArch: "s390x", lookup: nothingInstalled})
	if !plan.Complete() {
		t.Fatalf("noarch was rejected: %v", plan.Problems)
	}
}

// An already-present artifact is adopted whatever its architecture: the file is
// here and verified, so there is nothing to materialize.
func TestComputeDoesNotCheckArchForAnAdoption(t *testing.T) {
	root := projectRoot(t)
	app := artifact(t, "star/2.7.11b", "echo star\n")
	appPath := vendor(t, root, app, "echo star\n")

	l := lock.New()
	l.Pins["star/2.7.11b"] = lock.PinEntry{Artifact: appPath}

	plan := Compute(root, l, Options{HostArch: "s390x", lookup: installed("star/2.7.11b")})
	if !plan.Complete() {
		t.Fatalf("an adoption was blocked on architecture: %v", plan.Problems)
	}
}

// A rebuild that would prompt cannot run unattended, and planning says so
// before anything starts.
func TestComputeReportsInteractiveRebuilds(t *testing.T) {
	root := projectRoot(t)
	app := artifact(t, "star/2.7.11b", "echo star\n")
	app.Source.RequiresInput = true
	derived, err := key.Generate(app, key.Sources{meta.RecipeFileName: []byte("echo star\n")})
	if err != nil {
		t.Fatal(err)
	}
	app.Keys = derived.Keys()
	appPath := vendor(t, root, app, "echo star\n")

	l := lock.New()
	l.Pins["star/2.7.11b"] = lock.PinEntry{Artifact: appPath}

	if plan := Compute(root, l, Options{lookup: nothingInstalled}); plan.Complete() {
		t.Fatal("an interactive rebuild was planned without comment")
	}
	// Adopting it needs no rebuild, so the prompt never arises.
	if plan := Compute(root, l, Options{lookup: installed("star/2.7.11b")}); !plan.Complete() {
		t.Fatalf("an adopted artifact was blocked on #INPUT:: %v", plan.Problems)
	}
}

// An unsound lock stops planning: acquiring against an inconsistent
// specification would only fail later, further from the cause.
func TestComputeStopsOnAnInvalidLock(t *testing.T) {
	root := projectRoot(t)
	l := lock.New()
	l.Pins["ghost/1.0"] = lock.PinEntry{Artifact: lock.EntryPath("ghost--1.0@aaaaaaaaaaaa")}

	plan := Compute(root, l, Options{lookup: nothingInstalled})
	if plan.Complete() || len(plan.Steps) != 0 {
		t.Fatalf("plan = %#v", plan)
	}
}

// A diamond is planned once, with both parents depending on it.
func TestComputePlansADiamondOnce(t *testing.T) {
	root := projectRoot(t)
	shared := artifact(t, "zlib/1.3", "echo zlib\n")
	sharedPath := vendor(t, root, shared, "echo zlib\n")
	left := artifact(t, "index/1.0", "echo index\n", edge(shared))
	leftPath := vendor(t, root, left, "echo index\n")
	right := artifact(t, "other/1.0", "echo other\n", edge(shared))
	rightPath := vendor(t, root, right, "echo other\n")

	l := lock.New()
	l.Pins["index/1.0"] = lock.PinEntry{Artifact: leftPath}
	l.Pins["other/1.0"] = lock.PinEntry{Artifact: rightPath}

	plan := Compute(root, l, Options{lookup: nothingInstalled})
	if !plan.Complete() {
		t.Fatalf("problems: %v", plan.Problems)
	}
	if len(plan.Steps) != 3 {
		t.Fatalf("steps = %d, want the shared node planned once", len(plan.Steps))
	}
	if plan.Steps[0].Artifact != sharedPath {
		t.Errorf("order = %v, want the shared dependency first", plan.Steps[0].Artifact)
	}
}

// A name key must equal the artifact's own name, so two name keys can never
// select one artifact. Each direct step answers exactly one request.
func TestComputeGivesEachDirectStepOneRequest(t *testing.T) {
	root := projectRoot(t)
	app := artifact(t, "star/2.7.11b", "echo star\n")
	appPath := vendor(t, root, app, "echo star\n")
	other := artifact(t, "cutadapt/5.0", "echo cutadapt\n")
	otherPath := vendor(t, root, other, "echo cutadapt\n")

	l := lock.New()
	l.Pins["star/2.7.11b"] = lock.PinEntry{Artifact: appPath}
	l.Pins["cutadapt/5.0"] = lock.PinEntry{Artifact: otherPath}

	plan := Compute(root, l, Options{lookup: installed("star/2.7.11b", "cutadapt/5.0")})
	if !plan.Complete() {
		t.Fatalf("problems: %v", plan.Problems)
	}
	if len(plan.Steps) != 2 {
		t.Fatalf("steps = %d, want one per artifact", len(plan.Steps))
	}
	for _, step := range plan.Steps {
		if len(step.Requests) != 1 {
			t.Errorf("%s answers %v, want exactly one request", step.Name, step.Requests)
		}
	}
}

// Equivalence is the default: an artifact that can substitute for the locked one
// is adopted, and the plan says which identity it actually got.
func TestComputeAdoptsAnEquivalentSubstituteByDefault(t *testing.T) {
	root := projectRoot(t)
	app := artifact(t, "star/2.7.11b", "echo star\n")
	appPath := vendor(t, root, app, "echo star\n")

	l := lock.New()
	l.Pins["star/2.7.11b"] = lock.PinEntry{Artifact: appPath}

	plan := Compute(root, l, Options{lookup: substitute("star/2.7.11b")})
	if !plan.Complete() {
		t.Fatalf("problems: %v", plan.Problems)
	}
	if plan.Match != MatchEquivalent {
		t.Errorf("match = %q, want equivalent by default", plan.Match)
	}
	step := plan.Steps[0]
	if step.Action != ActionAdopt {
		t.Fatalf("action = %q, want adopt", step.Action)
	}
	if step.Identity != app.Keys.Identity.Digest() {
		t.Errorf("identity = %q, want the locked one", step.Identity)
	}
	if step.Found == "" || step.Found == step.Identity {
		t.Errorf("found = %q, want the substitute's own identity", step.Found)
	}
	if plan.Work() != 0 {
		t.Errorf("work = %d, want none: an equivalent copy is already here", plan.Work())
	}
}

// MatchIdentity refuses the substitute and plans to produce the exact build.
func TestComputeRefusesASubstituteUnderMatchIdentity(t *testing.T) {
	root := projectRoot(t)
	app := artifact(t, "star/2.7.11b", "echo star\n")
	appPath := vendor(t, root, app, "echo star\n")

	l := lock.New()
	l.Pins["star/2.7.11b"] = lock.PinEntry{Artifact: appPath}

	plan := Compute(root, l, Options{Match: MatchIdentity, lookup: substitute("star/2.7.11b")})
	if !plan.Complete() {
		t.Fatalf("problems: %v", plan.Problems)
	}
	if plan.Match != MatchIdentity {
		t.Errorf("match = %q", plan.Match)
	}
	step := plan.Steps[0]
	if step.Action != ActionBuild {
		t.Fatalf("action = %q, want build: the exact identity is absent", step.Action)
	}
	if step.Found != "" {
		t.Errorf("found = %q, want empty: nothing was adopted", step.Found)
	}
}

// An exact hit is never reported as a substitution, in either mode.
func TestComputeLeavesFoundEmptyForAnExactHit(t *testing.T) {
	root := projectRoot(t)
	app := artifact(t, "star/2.7.11b", "echo star\n")
	appPath := vendor(t, root, app, "echo star\n")

	l := lock.New()
	l.Pins["star/2.7.11b"] = lock.PinEntry{Artifact: appPath}

	for _, match := range []Match{MatchEquivalent, MatchIdentity} {
		plan := Compute(root, l, Options{Match: match, lookup: installed("star/2.7.11b")})
		step := plan.Steps[0]
		if step.Action != ActionAdopt || step.Found != "" {
			t.Errorf("%s: step = %#v, want a plain adoption", match, step)
		}
	}
}

// A `path:` pin is a project output, not something to find in the store.
// Resolving it against the images roots would report success while the declared
// path stays empty and the script still fails at mount time.
func TestComputePlansAPathSelectionAtItsProjectDestination(t *testing.T) {
	root := projectRoot(t)
	app := artifact(t, "tool/1.0", "echo tool\n")
	appPath := vendor(t, root, app, "echo tool\n")

	l := lock.New()
	l.Pins[lock.PathPrefix+"overlays/tool.sqf"] = lock.PinEntry{Artifact: appPath}

	// installed() would answer for this artifact if the store were consulted.
	plan := Compute(root, l, Options{lookup: installed("tool/1.0")})
	if !plan.Complete() {
		t.Fatalf("problems: %v", plan.Problems)
	}
	step := plan.Steps[0]
	if step.Destination != "overlays/tool.sqf" {
		t.Fatalf("destination = %q, want the declared project path", step.Destination)
	}
	if step.Action != ActionBuild {
		t.Errorf("action = %q, want build: nothing is at the project path, and a store copy is not a substitute", step.Action)
	}
	if step.Layout != "" {
		t.Errorf("layout = %q, want none: a project path is not a store namespace", step.Layout)
	}
}

// A project path already holding the right artifact is adopted in place.
func TestComputeAdoptsAnExistingProjectPath(t *testing.T) {
	root := projectRoot(t)
	app := artifact(t, "tool/1.0", "echo tool\n")
	appPath := vendor(t, root, app, "echo tool\n")

	l := lock.New()
	l.Pins[lock.PathPrefix+"overlays/tool.sqf"] = lock.PinEntry{Artifact: appPath}

	at := filepath.Join(root, "overlays", "tool.sqf")
	plan := Compute(root, l, Options{lookup: nothingInstalled, lookupAt: presentAt(at, app)})
	if !plan.Complete() {
		t.Fatalf("problems: %v", plan.Problems)
	}
	step := plan.Steps[0]
	if step.Action != ActionAdopt || step.Path != at {
		t.Fatalf("step = %#v, want adoption at %s", step, at)
	}
	if plan.Work() != 0 {
		t.Errorf("work = %d, want none", plan.Work())
	}
}

// A project path is the only place a restore can destroy a file, so one holding
// something the lock does not name is refused before anything is acquired.
func TestComputeRefusesToOverwriteAProjectPath(t *testing.T) {
	root := projectRoot(t)
	app := artifact(t, "tool/1.0", "echo tool\n")
	appPath := vendor(t, root, app, "echo tool\n")

	l := lock.New()
	l.Pins[lock.PathPrefix+"overlays/tool.sqf"] = lock.PinEntry{Artifact: appPath}

	at := filepath.Join(root, "overlays", "tool.sqf")
	if err := os.MkdirAll(filepath.Dir(at), 0o755); err != nil {
		t.Fatal(err)
	}
	if err := os.WriteFile(at, []byte("copied from somewhere else"), 0o644); err != nil {
		t.Fatal(err)
	}

	plan := Compute(root, l, Options{lookup: nothingInstalled, lookupAt: missingAt})
	if plan.Complete() {
		t.Fatal("the file at the project path would have been overwritten without a word")
	}
	if step := plan.Steps[0]; step.Replaces == "" {
		t.Error("replaces is empty, so a dry run cannot say what is about to be lost")
	}

	// --replace is the deliberate act, and it plans the overwrite.
	forced := Compute(root, l, Options{lookup: nothingInstalled, lookupAt: missingAt, Replace: true})
	if !forced.Complete() {
		t.Fatalf("problems with --replace: %v", forced.Problems)
	}
	if step := forced.Steps[0]; step.Action != ActionBuild || step.Replaces == "" {
		t.Fatalf("step = %#v, want a build that reports what it replaces", step)
	}
}

// Nothing is being replaced when the path already holds the locked artifact.
func TestComputeReportsNoReplacementWhenAdopting(t *testing.T) {
	root := projectRoot(t)
	app := artifact(t, "tool/1.0", "echo tool\n")
	appPath := vendor(t, root, app, "echo tool\n")

	l := lock.New()
	l.Pins[lock.PathPrefix+"overlays/tool.sqf"] = lock.PinEntry{Artifact: appPath}

	at := filepath.Join(root, "overlays", "tool.sqf")
	plan := Compute(root, l, Options{lookup: nothingInstalled, lookupAt: presentAt(at, app)})
	if step := plan.Steps[0]; step.Replaces != "" {
		t.Fatalf("replaces = %q, want none", step.Replaces)
	}
}

// Two paths pinning one artifact are two files to produce.
func TestComputePlansEachProjectDestinationSeparately(t *testing.T) {
	root := projectRoot(t)
	app := artifact(t, "tool/1.0", "echo tool\n")
	appPath := vendor(t, root, app, "echo tool\n")

	l := lock.New()
	l.Pins[lock.PathPrefix+"a/tool.sqf"] = lock.PinEntry{Artifact: appPath}
	l.Pins[lock.PathPrefix+"b/tool.sqf"] = lock.PinEntry{Artifact: appPath}

	plan := Compute(root, l, Options{lookup: nothingInstalled})
	if len(plan.Steps) != 2 {
		t.Fatalf("steps = %d, want one per destination", len(plan.Steps))
	}
	if plan.Steps[0].Destination != "a/tool.sqf" || plan.Steps[1].Destination != "b/tool.sqf" {
		t.Fatalf("destinations = %q, %q", plan.Steps[0].Destination, plan.Steps[1].Destination)
	}
}

// One artifact cannot be both store-addressed and pinned to a project path.
func TestComputeRefusesAnArtifactSelectedBothWays(t *testing.T) {
	root := projectRoot(t)
	app := artifact(t, "tool/1.0", "echo tool\n")
	appPath := vendor(t, root, app, "echo tool\n")

	l := lock.New()
	l.Pins["tool/1.0"] = lock.PinEntry{Artifact: appPath}
	l.Pins[lock.PathPrefix+"overlays/tool.sqf"] = lock.PinEntry{Artifact: appPath}

	plan := Compute(root, l, Options{lookup: nothingInstalled})
	if plan.Complete() {
		t.Fatal("a lock with two destinations for one artifact was accepted")
	}
	if !strings.Contains(strings.Join(plan.Problems, "\n"), "two destinations") {
		t.Errorf("problems = %v, want the conflict named", plan.Problems)
	}
}

// A closure dependency of a project-path root is an ordinary named artifact and
// still goes to the store — only the root itself is project-owned.
func TestComputeKeepsAPathRootsChildrenStoreBound(t *testing.T) {
	root := projectRoot(t)
	dep := artifact(t, "zlib/1.3", "echo zlib\n")
	depPath := vendor(t, root, dep, "echo zlib\n")
	data := artifact(t, "index/1.0", "echo index\n", edge(dep))
	dataPath := vendor(t, root, data, "echo index\n")

	l := lock.New()
	l.Pins[lock.PathPrefix+"overlays/index.sqf"] = lock.PinEntry{Artifact: dataPath}

	plan := Compute(root, l, Options{lookup: nothingInstalled})
	if !plan.Complete() {
		t.Fatalf("problems: %v", plan.Problems)
	}
	if len(plan.Steps) != 2 {
		t.Fatalf("steps = %d, want the dependency and the root", len(plan.Steps))
	}
	child, rootStep := plan.Steps[0], plan.Steps[1]
	if child.Artifact != depPath || child.Destination != "" {
		t.Errorf("child = %#v, want a store-bound dependency", child)
	}
	if rootStep.Destination != "overlays/index.sqf" {
		t.Errorf("root destination = %q", rootStep.Destination)
	}
}

// missingAt reports nothing usable at the path, whatever bytes are actually
// there.
func missingAt(string, string, meta.Keys, Match) (store.Candidate, bool) {
	return store.Candidate{}, false
}

// presentAt reports the artifact as already sitting at one exact path.
func presentAt(path string, manifest meta.Manifest) lookupAtFunc {
	return func(got, name string, keys meta.Keys, _ Match) (store.Candidate, bool) {
		if got != path || name != manifest.Name {
			return store.Candidate{}, false
		}
		return store.Candidate{Name: name, Path: got,
			Identity: keys.Identity, Equiv: keys.Equiv}, true
	}
}

// A closure entry is a build input. When its dependent is adopted, nothing ever
// opens it, so acquiring it would be pure waste.
func TestComputeSkipsAClosureInputNothingNeeds(t *testing.T) {
	root := projectRoot(t)
	dep := artifact(t, "samtools/1.21", "echo old\n")
	vendor(t, root, dep, "echo old\n")
	data := artifact(t, "index/1.0", "echo index\n", edge(dep))
	dataPath := vendor(t, root, data, "echo index\n")

	l := lock.New()
	l.Pins["index/1.0"] = lock.PinEntry{Artifact: dataPath}

	plan := Compute(root, l, Options{lookup: installed("index/1.0")})
	if !plan.Complete() {
		t.Fatalf("problems: %v", plan.Problems)
	}
	if len(plan.Steps) != 1 || plan.Steps[0].Name != "index/1.0" {
		t.Fatalf("steps = %#v, want only the adopted pin", plan.Steps)
	}
	if plan.Work() != 0 {
		t.Errorf("work = %d, want none", plan.Work())
	}
}

// The same input is planned once its dependent actually has to be built.
func TestComputeKeepsAClosureInputABuildNeeds(t *testing.T) {
	root := projectRoot(t)
	dep := artifact(t, "samtools/1.21", "echo old\n")
	depPath := vendor(t, root, dep, "echo old\n")
	data := artifact(t, "index/1.0", "echo index\n", edge(dep))
	dataPath := vendor(t, root, data, "echo index\n")

	l := lock.New()
	l.Pins["index/1.0"] = lock.PinEntry{Artifact: dataPath}

	plan := Compute(root, l, Options{lookup: nothingInstalled})
	if len(plan.Steps) != 2 {
		t.Fatalf("steps = %#v, want the input and its dependent", plan.Steps)
	}
	if plan.Steps[0].Artifact != depPath || plan.Steps[0].Direct {
		t.Errorf("first step = %#v, want the closure-only input", plan.Steps[0])
	}
}

// A pruned step's problems are not the plan's: an input nothing needs cannot
// refuse a restore for being built for another architecture.
func TestComputeIgnoresProblemsFromAPrunedInput(t *testing.T) {
	root := projectRoot(t)
	dep := artifact(t, "samtools/1.21", "echo old\n")
	dep.Platform.Arch = "s390x"
	vendor(t, root, dep, "echo old\n")
	data := artifact(t, "index/1.0", "echo index\n", edge(dep))
	dataPath := vendor(t, root, data, "echo index\n")

	l := lock.New()
	l.Pins["index/1.0"] = lock.PinEntry{Artifact: dataPath}

	plan := Compute(root, l, Options{lookup: installed("index/1.0")})
	if !plan.Complete() {
		t.Fatalf("a foreign-arch input nothing needs refused the plan: %v", plan.Problems)
	}
}

// A fetch downloads a finished artifact and opens none of its inputs, so it
// does not drag them in either.
func TestComputeDoesNotMaterializeInputsForAFetch(t *testing.T) {
	root := projectRoot(t)
	dep := artifact(t, "samtools/1.21", "echo old\n")
	vendor(t, root, dep, "echo old\n")
	data := artifact(t, "index/1.0", "echo index\n", edge(dep))
	dataPath := vendor(t, root, data, "echo index\n")

	l := lock.New()
	l.Pins["index/1.0"] = lock.PinEntry{Artifact: dataPath}
	digest := "sha256:" + strings.Repeat("b", 64)
	if err := l.AddRemote(dataPath, lock.Remote{Repository: "ghcr.io/x/index", ManifestDigest: digest}); err != nil {
		t.Fatal(err)
	}

	plan := Compute(root, l, Options{lookup: nothingInstalled})
	if len(plan.Steps) != 1 || plan.Steps[0].Action != ActionFetch {
		t.Fatalf("steps = %#v, want only the fetch", plan.Steps)
	}
}

// A build input is resolved by the role its edge records, not by its own keys:
// what may stand in is whatever leaves the dependent's equivalence unchanged.
func TestComputeResolvesABuildInputByItsEdgeRole(t *testing.T) {
	root := projectRoot(t)
	dep := artifact(t, "samtools/1.21", "echo samtools\n")
	vendor(t, root, dep, "echo samtools\n")
	data := artifact(t, "grch38/index", "echo index\n", history(dep))
	dataPath := vendor(t, root, data, "echo index\n")

	l := lock.New()
	l.Pins["grch38/index"] = lock.PinEntry{Artifact: dataPath}

	var asked meta.Dependency
	plan := Compute(root, l, Options{
		lookup: nothingInstalled,
		lookupInput: func(got meta.Dependency, _ Match, _ []string) (store.Candidate, bool) {
			asked = got
			// A different version entirely, which a history edge admits.
			return store.Candidate{Name: "samtools/1.25", Path: "/images/samtools--1.25.sqf",
				Layout: store.LayoutFlat}, true
		},
	})
	if asked.Role != meta.RoleHistory || asked.Name != "samtools/1.21" {
		t.Fatalf("edge = %#v, want the recorded history edge", asked)
	}
	byName := steps(plan)
	if got := byName["samtools/1.21"].Action; got != ActionAdopt {
		t.Errorf("input action = %q, want adopt — a newer version satisfies a history edge", got)
	}
	if got := byName["samtools/1.21"].Path; got != "/images/samtools--1.25.sqf" {
		t.Errorf("input path = %q, want the substitute", got)
	}
	if got := byName["grch38/index"].Action; got != ActionBuild {
		t.Errorf("dependent action = %q, want build", got)
	}
}

// A pin answers for its own keys however it is also depended on: someone
// asked for it by name, so the edge's latitude does not apply.
func TestComputeResolvesASelectionByItsOwnKeys(t *testing.T) {
	root := projectRoot(t)
	dep := artifact(t, "samtools/1.21", "echo samtools\n")
	depPath := vendor(t, root, dep, "echo samtools\n")
	data := artifact(t, "grch38/index", "echo index\n", history(dep))
	dataPath := vendor(t, root, data, "echo index\n")

	l := lock.New()
	l.Pins["grch38/index"] = lock.PinEntry{Artifact: dataPath}
	l.Pins["samtools/1.21"] = lock.PinEntry{Artifact: depPath}

	plan := Compute(root, l, Options{
		lookup: installed("samtools/1.21"),
		lookupInput: func(meta.Dependency, Match, []string) (store.Candidate, bool) {
			t.Error("a pinned artifact was resolved as a build input")
			return store.Candidate{}, false
		},
	})
	if got := steps(plan)["samtools/1.21"]; !got.Direct || got.Action != ActionAdopt {
		t.Fatalf("step = %#v, want a direct adoption", got)
	}
}

// One artifact depended on at two roles must satisfy the stricter claim, or the
// dependent that pinned its version would find its equivalence moved.
//
// App-versus-history is the case that arises: key.Role reads the dependent's
// name, so star/2.7.11b is an app to star/2.7.11b/index and history to
// grch38/star-gencode49. Data is never in tension with either — a data-typed
// dependency is always RoleData, which key generation enforces.
func TestInboundEdgesKeepsTheStrictestRole(t *testing.T) {
	root := projectRoot(t)
	shared := artifact(t, "star/2.7.11b", "echo star\n")
	vendor(t, root, shared, "echo star\n")
	loose := artifact(t, "grch38/star-gencode49", "echo loose\n", history(shared))
	vendor(t, root, loose, "echo loose\n")
	strict := artifact(t, "star/2.7.11b/index", "echo strict\n", edge(shared))
	vendor(t, root, strict, "echo strict\n")

	l := lock.New()
	l.Pins["grch38/star-gencode49"] = lock.PinEntry{Artifact: lock.EntryPath(entryName(loose))}
	l.Pins["star/2.7.11b/index"] = lock.PinEntry{Artifact: lock.EntryPath(entryName(strict))}

	verified, problems := lock.Verify(root, l)
	if len(problems) > 0 {
		t.Fatalf("problems = %v", problems)
	}
	got, ok := inboundEdges(verified)[lock.EntryPath(entryName(shared))]
	if !ok {
		t.Fatal("the shared dependency has no inbound edge")
	}
	if got.Role != meta.RoleApp {
		t.Fatalf("role = %q, want the stricter app claim", got.Role)
	}
}

// history is one dependency edge demoted to build history; edge() is the app
// case, and a data edge is only valid on a data-typed dependency.
func history(manifest meta.Manifest) meta.Dependency {
	dep := edge(manifest)
	dep.Role = meta.RoleHistory
	return dep
}

// --only is what a submitted job carries: the artifact it was sent for, plus
// what that artifact needs to be built.
func TestComputeRestrictsToOneArtifactAndItsClosure(t *testing.T) {
	root := projectRoot(t)
	dep := artifact(t, "zlib/1.3", "echo zlib\n")
	depPath := vendor(t, root, dep, "echo zlib\n")
	wanted := artifact(t, "index/1.0", "echo index\n", edge(dep))
	wantedPath := vendor(t, root, wanted, "echo index\n")
	other := artifact(t, "samtools/1.23.1", "echo samtools\n")
	otherPath := vendor(t, root, other, "echo samtools\n")

	l := lock.New()
	l.Pins["index/1.0"] = lock.PinEntry{Artifact: wantedPath}
	l.Pins["samtools/1.23.1"] = lock.PinEntry{Artifact: otherPath}

	plan := Compute(root, l, Options{lookup: nothingInstalled, Only: wantedPath})
	if !plan.Complete() {
		t.Fatalf("problems: %v", plan.Problems)
	}
	byName := steps(plan)
	if len(plan.Steps) != 2 {
		t.Fatalf("steps = %v, want only the named artifact and its dependency", byName)
	}
	if _, ok := byName["samtools/1.23.1"]; ok {
		t.Error("an unrelated pin survived --only")
	}
	if got := byName["zlib/1.3"]; got.Artifact != depPath {
		t.Errorf("the dependency was dropped: %#v", got)
	}
}

// A pin that is also a dependency stays a pin: narrowing what a
// restore covers says nothing about what its members are.
func TestComputeKeepsDirectnessUnderOnly(t *testing.T) {
	root := projectRoot(t)
	dep := artifact(t, "zlib/1.3", "echo zlib\n")
	depPath := vendor(t, root, dep, "echo zlib\n")
	wanted := artifact(t, "index/1.0", "echo index\n", edge(dep))
	wantedPath := vendor(t, root, wanted, "echo index\n")

	l := lock.New()
	l.Pins["index/1.0"] = lock.PinEntry{Artifact: wantedPath}
	l.Pins["zlib/1.3"] = lock.PinEntry{Artifact: depPath}

	plan := Compute(root, l, Options{lookup: nothingInstalled, Only: wantedPath})
	if got := steps(plan)["zlib/1.3"]; !got.Direct {
		t.Errorf("a selected dependency became a build dependency under --only: %#v", got)
	}
}

func TestComputeRefusesAnOnlyItDoesNotCover(t *testing.T) {
	root := projectRoot(t)
	wanted := artifact(t, "index/1.0", "echo index\n")
	wantedPath := vendor(t, root, wanted, "echo index\n")

	l := lock.New()
	l.Pins["index/1.0"] = lock.PinEntry{Artifact: wantedPath}

	plan := Compute(root, l, Options{lookup: nothingInstalled, Only: "provenance/absent"})
	if plan.Complete() {
		t.Fatal("--only naming nothing was accepted")
	}
	if len(plan.Steps) != 0 {
		t.Errorf("steps = %v, want none planned", plan.Steps)
	}
}

// A problem belonging to an artifact --only excluded is not this restore's.
func TestComputeDropsProblemsFromOutsideOnly(t *testing.T) {
	root := projectRoot(t)
	wanted := artifact(t, "index/1.0", "echo index\n")
	wantedPath := vendor(t, root, wanted, "echo index\n")
	interactive := artifact(t, "genome/1.0", "echo genome\n")
	interactive.Source.RequiresInput = true
	interactivePath := vendor(t, root, interactive, "echo genome\n")

	l := lock.New()
	l.Pins["index/1.0"] = lock.PinEntry{Artifact: wantedPath}
	l.Pins["genome/1.0"] = lock.PinEntry{Artifact: interactivePath}

	if plan := Compute(root, l, Options{lookup: nothingInstalled}); plan.Complete() {
		t.Fatal("the interactive rebuild was expected to be a problem")
	}
	plan := Compute(root, l, Options{lookup: nothingInstalled, Only: wantedPath})
	if !plan.Complete() {
		t.Errorf("problems from outside --only were reported: %v", plan.Problems)
	}
}

// stepsFor returns every step for one artifact path, since a name can be shared
// by two identities and a path by two destinations.
func stepsFor(plan *Plan, artifact string) []Step {
	var out []Step
	for _, step := range plan.Steps {
		if step.Artifact == artifact {
			out = append(out, step)
		}
	}
	return out
}

// Two identities of one name contend for the bare name. It belongs to the
// pin, not to whatever a rebuild happened to need — and dependencies are
// built first, so build order alone would decide it backwards.
func TestComputeGivesTheFlatNameToTheSelectionNotTheBuildDependency(t *testing.T) {
	root := projectRoot(t)
	older := artifact(t, "samtools/1.23.1", "echo older\n")
	olderPath := vendor(t, root, older, "echo older\n")
	genome := artifact(t, "grch38/genome", "echo genome\n", edge(older))
	genomePath := vendor(t, root, genome, "echo genome\n")
	selected := artifact(t, "samtools/1.23.1", "echo selected\n")
	selectedPath := vendor(t, root, selected, "echo selected\n")

	l := lock.New()
	l.Pins["grch38/genome"] = lock.PinEntry{Artifact: genomePath}
	l.Pins["samtools/1.23.1"] = lock.PinEntry{Artifact: selectedPath}

	plan := Compute(root, l, Options{lookup: nothingInstalled, KeepBuildDeps: true})
	if !plan.Complete() {
		t.Fatalf("problems: %v", plan.Problems)
	}
	dependency := stepsFor(plan, olderPath)
	if len(dependency) != 1 {
		t.Fatalf("steps for the build dependency = %v", dependency)
	}
	if !dependency[0].StoreOnly {
		t.Error("the build dependency was left free to take the name its pin answers to")
	}
	pin := stepsFor(plan, selectedPath)
	if len(pin) != 1 || pin[0].StoreOnly {
		t.Errorf("the pin yielded its own name: %#v", pin)
	}
}

// Nothing else claims the name, so the dependency is entitled to it.
func TestComputeLeavesAnUncontestedBuildDependencyFlat(t *testing.T) {
	root := projectRoot(t)
	dep := artifact(t, "zlib/1.3", "echo zlib\n")
	depPath := vendor(t, root, dep, "echo zlib\n")
	genome := artifact(t, "grch38/genome", "echo genome\n", edge(dep))
	genomePath := vendor(t, root, genome, "echo genome\n")

	l := lock.New()
	l.Pins["grch38/genome"] = lock.PinEntry{Artifact: genomePath}

	plan := Compute(root, l, Options{lookup: nothingInstalled, KeepBuildDeps: true})
	if got := stepsFor(plan, depPath); len(got) != 1 || got[0].StoreOnly {
		t.Errorf("an uncontested build dependency was pushed into the store: %#v", got)
	}
}

// A discarded build dependency is never installed, so it has no name to yield.
func TestComputeLeavesADiscardedBuildDependencyUnmarked(t *testing.T) {
	root := projectRoot(t)
	older := artifact(t, "samtools/1.23.1", "echo older\n")
	olderPath := vendor(t, root, older, "echo older\n")
	genome := artifact(t, "grch38/genome", "echo genome\n", edge(older))
	genomePath := vendor(t, root, genome, "echo genome\n")
	selected := artifact(t, "samtools/1.23.1", "echo selected\n")
	selectedPath := vendor(t, root, selected, "echo selected\n")

	l := lock.New()
	l.Pins["grch38/genome"] = lock.PinEntry{Artifact: genomePath}
	l.Pins["samtools/1.23.1"] = lock.PinEntry{Artifact: selectedPath}

	plan := Compute(root, l, Options{lookup: nothingInstalled})
	if got := stepsFor(plan, olderPath); len(got) != 1 || got[0].StoreOnly {
		t.Errorf("a transient step was marked store-only: %#v", got)
	}
}
