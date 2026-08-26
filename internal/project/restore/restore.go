package restore

import (
	"context"
	"errors"
	"fmt"
	"os"
	"path/filepath"
	"strings"

	ocispec "github.com/opencontainers/image-spec/specs-go/v1"

	"github.com/Justype/condatainer/internal/artifact/capsule"
	"github.com/Justype/condatainer/internal/artifact/compare"
	"github.com/Justype/condatainer/internal/artifact/meta"
	"github.com/Justype/condatainer/internal/build"
	"github.com/Justype/condatainer/internal/conda"
	"github.com/Justype/condatainer/internal/config"
	"github.com/Justype/condatainer/internal/image/producer"
	"github.com/Justype/condatainer/internal/logging"
	"github.com/Justype/condatainer/internal/project/lock"
	"github.com/Justype/condatainer/internal/registry"
	"github.com/Justype/condatainer/internal/store"
	"github.com/Justype/condatainer/internal/utils"
)

// ErrIncomplete reports that a restore did not make every pin available.
var ErrIncomplete = errors.New("restore is incomplete")

// ErrNotAcquirable reports an acquisition this build cannot perform.
var ErrNotAcquirable = errors.New("cannot acquire artifact")

// ErrJobsSubmitted reports that the restore handed work to the scheduler, so
// the project is not available yet. It is not a failure: re-running the restore
// after the jobs finish adopts what they published, and is also how a partly
// finished submission is resumed.
var ErrJobsSubmitted = errors.New("restore submitted scheduler jobs")

// Outcome is how one artifact was made available.
type Outcome string

const (
	// OutcomeAdopted found it already here.
	OutcomeAdopted Outcome = "adopted"
	// OutcomeBuilt rebuilt it from the vendored sources.
	OutcomeBuilt Outcome = "built"
	// OutcomeFetched downloaded it from a recorded remote.
	OutcomeFetched Outcome = "fetched"
)

// Result is one artifact that is now available, and how it got that way.
type Result struct {
	Artifact string `json:"artifact"`
	Name     string `json:"name"`
	// Identity is what the lock records; Found is what is actually there, set
	// only when an equivalent artifact stands in.
	Identity string          `json:"identity"`
	Found    string          `json:"found,omitempty"`
	Verdict  compare.Verdict `json:"verdict"`
	Outcome  Outcome         `json:"outcome"`
	// Path is the absolute image, and Destination the project-relative path it
	// was required at, empty for a store-addressed artifact.
	Path        string       `json:"path"`
	Destination string       `json:"destination,omitempty"`
	Layout      store.Layout `json:"layout,omitempty"`
	// Transient reports that this artifact was produced only to build something
	// else and no longer exists. --keep-build-deps installs it instead.
	Transient bool `json:"transient,omitempty"`
}

// Failure is one artifact that could not be made available.
type Failure struct {
	Artifact string `json:"artifact"`
	Name     string `json:"name"`
	Reason   string `json:"reason"`
	// Diffs name the inputs that moved, when a result was produced and rejected.
	Diffs []string `json:"diffs,omitempty"`
	// Rejected is where a rejected result was kept, when it was kept.
	Rejected string `json:"rejected,omitempty"`
}

// Blocked is one artifact that was never attempted because something beneath it
// failed. It is kept apart from Failure so one cause reads as one cause: a
// broken dependency produces a single failure and a list of what it stopped,
// rather than a failure per artifact with the real one buried among them.
type Blocked struct {
	Artifact string `json:"artifact"`
	Name     string `json:"name"`
	// Cause is the artifact that actually failed, which may lie further down
	// than this one's own edges.
	Cause string `json:"cause"`
}

// Report is everything one restore did.
type Report struct {
	Root     string    `json:"root"`
	Match    Match     `json:"match"`
	Results  []Result  `json:"results,omitempty"`
	Failures []Failure `json:"failures,omitempty"`
	Blocked  []Blocked `json:"blocked,omitempty"`
	// Submitted are the rebuilds the scheduler will run. They are not results:
	// nothing is available until those jobs finish.
	Submitted []Submitted `json:"submitted,omitempty"`
	// Problems are planning failures, which stop a restore before it acquires
	// anything.
	Problems []string `json:"problems,omitempty"`
}

// Complete reports whether every pin is now available. A submission is
// not completion: the artifact exists only once its job has run.
func (r *Report) Complete() bool {
	return len(r.Failures) == 0 && len(r.Blocked) == 0 && len(r.Problems) == 0 &&
		len(r.Submitted) == 0
}

// Run makes every locked pin available, dependency-first.
//
// Atomicity is per artifact, not across the project. A later failure leaves
// earlier results in place; they are exact results someone can use, not a
// partial transaction to roll back, and re-running adopts them.
func Run(ctx context.Context, root string, l *lock.Lock, opts Options) (*Report, error) {
	plan := Compute(root, l, opts)
	report := &Report{Root: plan.Root, Match: plan.Match, Problems: plan.Problems}
	if !plan.Complete() {
		return report, fmt.Errorf("%w: %d planning problem(s)", ErrIncomplete, len(plan.Problems))
	}

	verified, problems := lock.Verify(root, l)
	if len(problems) > 0 {
		// Compute already verified, so this cannot normally differ. It is read
		// again because the sources are needed, and a checkout that changed in
		// between must not be built from.
		for _, problem := range problems {
			report.Problems = append(report.Problems, problem.String())
		}
		return report, fmt.Errorf("%w: the checkout changed during planning", ErrIncomplete)
	}

	// One directory for every build dependency this restore has to produce, removed
	// when it returns — on success, on failure, and on cancellation alike.
	scratch := &transientRoot{}
	defer scratch.remove()

	// Which rebuilds the scheduler runs, decided before the first step so a
	// build dependency knows whether it is produced here or inside a job.
	queue := planJobs(ctx, root, verified, plan, opts)

	// Where each artifact ended up, so a dependent mounts what was just made
	// rather than resolving its name again.
	available := map[string]string{}
	// Which artifacts cannot be attempted, mapped to the artifact that failed
	// beneath them. Blockage propagates, so a chain reports one cause.
	stopped := map[string]string{}

	for _, step := range plan.Steps {
		if err := ctx.Err(); err != nil {
			return report, err
		}
		key := stepKey(step)
		// Produced inside a submitted job instead of here, together with the
		// dependent that needs it.
		if queue.deferred[key] {
			continue
		}
		if cause := blockedBy(step, stopped); cause != "" {
			report.Blocked = append(report.Blocked, Blocked{
				Artifact: step.Artifact, Name: step.Name, Cause: cause})
			stopped[step.Artifact] = cause
			continue
		}
		if queue.submit[key] {
			entry, ok := verified.Entries[step.Artifact]
			if !ok {
				report.Failures = append(report.Failures, Failure{Artifact: step.Artifact,
					Name: step.Name, Reason: "the artifact is no longer vendored"})
				stopped[step.Artifact] = step.Artifact
				continue
			}
			submitted, result, failure := queue.submitStep(ctx, root, entry, step, queue.waitFor(step), opts)
			switch {
			case failure != nil:
				report.Failures = append(report.Failures, *failure)
				stopped[step.Artifact] = step.Artifact
			case submitted != nil:
				report.Submitted = append(report.Submitted, *submitted)
			default:
				// Installed between planning and submission: adopted, not queued.
				report.Results = append(report.Results, *result)
				available[step.Artifact] = result.Path
			}
			continue
		}
		result, failure := execute(ctx, root, verified, step, plan.Match, available, scratch, opts)
		if failure != nil {
			report.Failures = append(report.Failures, *failure)
			stopped[step.Artifact] = step.Artifact
			continue
		}
		report.Results = append(report.Results, *result)
		// A destination step is a project output; nothing depends on it by
		// identity, and recording it would let it satisfy a store edge.
		if step.Destination == "" {
			available[step.Artifact] = result.Path
		}
	}
	switch {
	case len(report.Failures) > 0 || len(report.Blocked) > 0:
		return report, fmt.Errorf("%w: %d artifact(s) unavailable, %d blocked",
			ErrIncomplete, len(report.Failures), len(report.Blocked))
	case len(report.Submitted) > 0:
		return report, fmt.Errorf("%w: %d job(s)", ErrJobsSubmitted, len(report.Submitted))
	}
	return report, nil
}

// blockedBy names the artifact that failed beneath this step, or "".
//
// Only a build propagates blockage. An adopted or fetched artifact opens none of
// its dependencies, so one that could not be produced does not stop it — and
// pruning has usually dropped that dependency already.
func blockedBy(step Step, stopped map[string]string) string {
	if step.Action != ActionBuild {
		return ""
	}
	for _, dependency := range step.DependsOn {
		if cause, ok := stopped[dependency]; ok {
			return cause
		}
	}
	return ""
}

// transientRoot is the restore-scoped directory build dependencies are produced in,
// created on first use.
//
// It sits under the stable writable tmp root rather than $TMPDIR: several conda
// environments is more than node-local scratch holds, and a job sweep must not
// remove it while the restore is still mounting from it.
type transientRoot struct{ path string }

func (t *transientRoot) dir() (string, error) {
	if t.path != "" {
		return t.path, nil
	}
	base := config.GetWritableTmpDir()
	if err := utils.MkdirAllShared(base); err != nil {
		return "", err
	}
	path, err := os.MkdirTemp(base, "cnt-restore-deps-")
	if err != nil {
		return "", err
	}
	t.path = path
	return path, nil
}

func (t *transientRoot) remove() {
	if t.path != "" {
		os.RemoveAll(t.path) //nolint:errcheck
	}
}

func execute(ctx context.Context, root string, verified *lock.Verified, step Step, match Match,
	available map[string]string, scratch *transientRoot, opts Options) (*Result, *Failure) {

	entry, ok := verified.Entries[step.Artifact]
	if !ok {
		return nil, &Failure{Artifact: step.Artifact, Name: step.Name,
			Reason: "the artifact is no longer vendored"}
	}
	result := &Result{
		Artifact: step.Artifact, Name: step.Name, Identity: step.Identity,
		Found: step.Found, Destination: step.Destination, Verdict: compare.Exact,
	}
	if step.Found != "" {
		result.Verdict = compare.Equivalent
	}

	switch step.Action {
	case ActionAdopt:
		result.Outcome, result.Path, result.Layout = OutcomeAdopted, step.Path, step.Layout
		return result, nil
	case ActionFetch:
		return fetch(ctx, root, entry, step, match, result)
	}
	return rebuild(ctx, root, entry, step, match, available, scratch, opts, result)
}

// rebuild produces one artifact from its vendored sources and puts it where its
// role says: a project destination, an images root, or the restore's transient
// directory when nothing pinned it.
func rebuild(ctx context.Context, root string, entry *lock.Entry, step Step, match Match,
	available map[string]string, scratch *transientRoot, opts Options, result *Result) (*Result, *Failure) {

	fail := func(format string, args ...any) *Failure {
		return &Failure{Artifact: step.Artifact, Name: step.Name, Reason: fmt.Sprintf(format, args...)}
	}

	sources, err := lock.Sources(root, entry)
	if err != nil {
		return nil, fail("cannot read the vendored sources: %v", err)
	}
	deps, err := dependencyPaths(entry, available)
	if err != nil {
		return nil, fail("%v", err)
	}

	// A project path is a target two restores can race for, and the one a
	// submitted job was sent to produce. The guard serializes them and is what
	// the job adopts, by the job ID whoever submitted it recorded. A store
	// artifact is guarded by its own transaction instead.
	if step.Destination != "" {
		destination := filepath.Join(root, filepath.FromSlash(step.Destination))
		if err := utils.MkdirAllShared(filepath.Dir(destination)); err != nil {
			return nil, fail("%v", err)
		}
		guard, err := producer.AcquireLocal(destination)
		if err != nil {
			return nil, fail("%v", err)
		}
		defer guard.Release() //nolint:errcheck
	}

	transient := isTransient(step, opts)
	var staging string
	if transient {
		// Produced inside the restore's own directory, so it survives until every
		// dependent has mounted it and goes when the restore does. Named for the
		// artifact rather than randomly: the tree is what someone reads when a
		// build fails, and one entry per identity cannot collide within a restore.
		dir, dirErr := scratch.dir()
		if dirErr != nil {
			return nil, fail("cannot create a directory for build dependencies: %v", dirErr)
		}
		staging = filepath.Join(dir, capsule.EntryName(entry.Manifest.Name, entry.Identity.Digest()))
		if err = utils.MkdirAllShared(staging); err != nil {
			return nil, fail("cannot create a staging directory: %v", err)
		}
	} else {
		if staging, err = os.MkdirTemp(stagingDir(root, step), ".cnt-restore-"); err != nil {
			return nil, fail("cannot create a staging directory: %v", err)
		}
		defer os.RemoveAll(staging) //nolint:errcheck
	}
	output := filepath.Join(staging, "artifact.sqf")

	log := logging.FromContext(ctx)
	var buildErr error
	for attempt, condaSource := range condaSources(entry, sources, match) {
		if attempt > 0 {
			// A fresh output: the locked build refuses an occupied one, and the
			// failed attempt may have left something behind.
			output = filepath.Join(staging, fmt.Sprintf("artifact-%d.sqf", attempt))
			log.Info("retrying the Conda replay from the recorded environment", "kind", "note",
				"name", step.Name, "source", condaSource, "err", buildErr)
		}
		object, err := build.NewLockedObject(ctx, build.LockedSpec{
			Manifest:    entry.Manifest,
			Sources:     sources,
			Deps:        deps,
			Output:      output,
			Answers:     opts.Answers[step.Artifact],
			CondaSource: condaSource,
		})
		if err != nil {
			return nil, fail("%v", err)
		}
		log.Info("rebuilding from the lock", "kind", "note",
			"name", step.Name, "identity", step.Identity)
		if buildErr = object.Build(ctx, false); buildErr == nil {
			break
		}
		if ctx.Err() != nil {
			break
		}
	}
	if buildErr != nil {
		return nil, fail("rebuild failed: %v", buildErr)
	}

	verdict, diffs, err := verify(entry, output, match)
	if err != nil {
		return nil, fail("cannot verify the rebuilt artifact: %v", err)
	}
	if verdict != compare.Exact && verdict != compare.Equivalent {
		return nil, reject(ctx, step, transient, output, verdict, diffs)
	}
	result.Verdict = verdict
	result.Outcome = OutcomeBuilt

	if transient {
		result.Path, result.Transient = output, true
		return result, nil
	}
	if step.Destination != "" {
		path, err := placeAt(ctx, filepath.Join(root, filepath.FromSlash(step.Destination)), output)
		if err != nil {
			return nil, fail("%v", err)
		}
		result.Path = path
		return result, nil
	}

	candidate, err := install(entry, step, output)
	if err != nil {
		return nil, fail("cannot install: %v", err)
	}
	result.Path, result.Layout = candidate.Path, candidate.Layout
	if candidate.Identity != entry.Identity {
		result.Found = candidate.Identity.Digest()
	}
	return result, nil
}

// isTransient reports that a rebuild is scaffolding: a closure node no pin
// names, produced only so its dependent can be built.
//
// Nothing asked for it by name, and nothing needs it afterwards — an artifact
// bakes its inputs in, so a dependency image is needed to rebuild it and never
// to use it. --keep-build-deps installs it through the ordinary rule instead.
func isTransient(step Step, opts Options) bool {
	return !step.Direct && step.Destination == "" && !opts.KeepBuildDeps
}

// stagingDir picks where a rebuild is staged so publication is a rename on one
// filesystem. A project destination stages beside itself; a store-bound
// artifact stages in the configured temporary root, since where it lands is not
// known until the transaction picks it.
func stagingDir(root string, step Step) string {
	if step.Destination == "" {
		return ""
	}
	dir := filepath.Dir(filepath.Join(root, filepath.FromSlash(step.Destination)))
	if err := utils.MkdirAllShared(dir); err != nil {
		return ""
	}
	return dir
}

// dependencyPaths pairs each manifest edge with the image that satisfied it, in
// the manifest's own order.
//
// Order is the contract: the rebuild re-derives its edges from what it mounts,
// so a different order produces a different identity rather than a wrong
// artifact. Nothing is resolved by name — every path came from a completed step.
func dependencyPaths(entry *lock.Entry, available map[string]string) ([]build.LockedDep, error) {
	var deps []build.LockedDep
	for _, dep := range entry.Manifest.Dependencies {
		artifact := lock.EntryPath(capsule.EntryName(dep.Name, dep.Identity.Digest()))
		path, ok := available[artifact]
		if !ok {
			return nil, fmt.Errorf("dependency %s was not restored before its dependent", dep.Name)
		}
		deps = append(deps, build.LockedDep{Name: dep.Name, Identity: dep.Identity, Path: path})
	}
	return deps, nil
}

// condaSources is the vendored Conda exports to replay, in the order they are
// worth trying. Every other build type gets a single empty entry, meaning "the
// source this build type already has".
//
// explicit.txt is always first: it names exact package URLs and is the only
// input that reproduces the recorded identity. Those URLs rot — bioconda prunes
// old builds — and under MatchEquivalent environment.yml is the second chance.
// It is `--no-builds` output, so it pins every transitive package at an exact
// version, the interpreter included: a solve from it cannot land cutadapt on a
// different Python, only on different build strings. That is precisely what the
// equivalence key tolerates, and the file *is* that key's preimage, so a result
// matching it satisfies the mode by construction.
//
// There is no third tier. Re-solving from the recipe's original request would
// change both keys, so it could not answer the lock at all; that is a re-pin,
// not a restore.
func condaSources(entry *lock.Entry, sources map[string][]byte, match Match) []string {
	if entry.Manifest.BuildType != build.BuildTypeConda.String() {
		return []string{""}
	}
	ordered := []string{conda.ExplicitFileName}
	if match.Normalize() == MatchEquivalent {
		if _, ok := sources[conda.EnvironmentFileName]; ok {
			ordered = append(ordered, conda.EnvironmentFileName)
		}
	}
	return ordered
}

// verify reads what was produced and reports how it relates to the lock.
func verify(entry *lock.Entry, path string, match Match) (compare.Verdict, []string, error) {
	got, err := compare.Read(path)
	if err != nil {
		return compare.Unverifiable, nil, err
	}
	want := compare.Artifact{
		Name: entry.Manifest.Name, Type: entry.Manifest.Type,
		Arch: entry.Manifest.Platform.Arch, Format: entry.Manifest.BuildType,
		// Digest(), not SHA256: compare.Read reports what it regenerated as
		// "sha256:<hex>", and a bare hex here compares unequal to every artifact.
		Identity: entry.Identity.Digest(), Equiv: entry.Equiv.Digest(),
		IdentityScheme: entry.Identity.Scheme, EquivScheme: entry.Equiv.Scheme,
		Dependencies: entry.Manifest.Dependencies,
	}
	outcome := compare.Compare(want, got)
	diffs := make([]string, 0, len(outcome.Diffs))
	for _, diff := range outcome.Diffs {
		diffs = append(diffs, diff.String())
	}
	if outcome.Reason != "" {
		diffs = append(diffs, outcome.Reason)
	}
	// Under MatchIdentity an equivalent result is not acceptable, but it is
	// still worth naming as equivalent rather than as merely different.
	if outcome.Verdict == compare.Equivalent && match.Normalize() == MatchIdentity {
		return compare.Different, append(diffs,
			"the rebuild is equivalent but not identical, and --match identity was requested"), nil
	}
	return outcome.Verdict, diffs, nil
}

// reject records a result that does not answer the lock.
//
// It is never installed under the locked name and never renamed into a project
// destination — that substitution is the whole thing a lock prevents. A
// store-bound one is kept under its own true identity, because the store is
// addressed by identity and an honest entry there is never wrong: it makes the
// remedy `project pin` rather than a manual rebuild. A project
// destination has no such address, and a build dependency was never something to
// keep, so both are dropped with the staging directory.
func reject(ctx context.Context, step Step, transient bool, output string,
	verdict compare.Verdict, diffs []string) *Failure {

	failure := &Failure{Artifact: step.Artifact, Name: step.Name, Diffs: diffs,
		Reason: fmt.Sprintf("the rebuild is %s, not the locked artifact", verdict)}
	if step.Destination != "" || transient {
		return failure
	}
	produced, err := compare.Read(output)
	if err != nil {
		return failure
	}
	identity := produced.IdentityRef()
	candidate, err := store.InstallFile(produced.Name, identity, output, store.BeginOptions{})
	if err != nil {
		logging.FromContext(ctx).Debug("could not keep the rejected rebuild", "name", step.Name, "err", err)
		return failure
	}
	failure.Rejected = candidate.Path
	failure.Reason += fmt.Sprintf("; it was kept as %s, which `condatainer project pin` can pin",
		identity.Digest())
	return failure
}

// fetch acquires one artifact from the exact locations the lock recorded, in
// recorded order, and puts the result where the step's role says.
//
// It is a fetch and never a search: every repository and platform digest came
// from the lock, so nothing enumerates tags and nothing asks a registry to find
// an identity. A location that no longer resolves is skipped and the next one
// tried; every location failing is an acquisition fault reported with what each
// one said, not a silent fall back to building — `--no-prebuilt` is how a caller
// asks for the build instead.
func fetch(ctx context.Context, root string, entry *lock.Entry, step Step, match Match, result *Result) (*Result, *Failure) {
	log := logging.FromContext(ctx)
	fail := func(format string, args ...any) *Failure {
		return &Failure{Artifact: step.Artifact, Name: step.Name, Reason: fmt.Sprintf(format, args...)}
	}

	var attempts []string
	for _, remote := range step.Remotes {
		where := remote.Repository + "@" + remote.ManifestDigest
		path, err := fetchFrom(ctx, root, entry, step, match, remote)
		if err == nil {
			log.Info("fetched", "kind", "note", "name", step.Name, "from", where)
			result.Outcome = OutcomeFetched
			result.Path = path.Path
			result.Layout = path.Layout
			if !path.Identity.Empty() && path.Identity != entry.Identity {
				result.Found = path.Identity.Digest()
			}
			return result, nil
		}
		if ctx.Err() != nil {
			return nil, fail("%v", ctx.Err())
		}
		log.Debug("remote did not serve the artifact", "name", step.Name, "from", where, "err", err)
		attempts = append(attempts, fmt.Sprintf("%s: %v", where, err))
	}
	// Distinct from "nothing was recorded", which the planner turns into a build
	// rather than routing here at all — so this reports the acquisition fault and
	// names the flag that asks for the build, instead of quietly taking it.
	return nil, fail("%v: no recorded remote served %s (%s); use --no-prebuilt to build from source instead",
		ErrNotAcquirable, step.Identity, strings.Join(attempts, "; "))
}

// fetched is where one successful fetch landed.
type fetched struct {
	Path     string
	Layout   store.Layout
	Identity meta.KeyRef
}

// fetchFrom resolves one recorded location, checks it before spending the
// bytes, downloads it, and publishes it.
//
// The identity is compared against the lock at the manifest, before any blob
// moves: a digest that has been re-pointed at different content is refused here
// having transferred nothing. It is checked again after assembly by the store
// transaction, from the payload's own regenerated keys — annotations are a fast
// index over the artifact, never a substitute for it.
func fetchFrom(ctx context.Context, root string, entry *lock.Entry, step Step, match Match, remote lock.Remote) (fetched, error) {
	base, repo, err := registry.SplitCoordinate(remote.Repository)
	if err != nil {
		return fetched{}, err
	}
	desc, annotations, err := registry.ResolveArtifact(ctx, base, repo, remote.ManifestDigest)
	if err != nil {
		return fetched{}, err
	}
	if err := registry.Check(annotations, registry.Want{Name: entry.Manifest.Name}); err != nil {
		return fetched{}, err
	}
	if err := acceptable(entry, match, registry.Identity(annotations), registry.Equiv(annotations)); err != nil {
		return fetched{}, err
	}

	if step.Destination != "" {
		return fetchToDestination(ctx, root, entry, step, match, base, repo, desc, annotations)
	}
	return fetchToStore(ctx, entry, step, base, repo, desc, annotations)
}

// acceptable reports whether the keys a published manifest advertises satisfy
// the lock under the active mode, before anything is downloaded.
//
// Under identity only the recorded build will do. Under equivalent an artifact
// carrying the locked equivalence substitutes, which is the whole point of the
// second key — but the exact identity still wins where a location offers it,
// because the caller tries the recorded remotes in order and each one either
// matches exactly or is a substitution reported as such.
func acceptable(entry *lock.Entry, match Match, identity, equiv meta.KeyRef) error {
	if identity == entry.Identity {
		return nil
	}
	if match == MatchIdentity {
		return fmt.Errorf("%w: it publishes %s, the lock pins %s",
			registry.ErrMismatch, describeKey(identity), entry.Identity.Digest())
	}
	if equiv.Empty() || equiv != entry.Equiv {
		return fmt.Errorf("%w: it publishes equivalence %s, the lock pins %s",
			registry.ErrMismatch, describeKey(equiv), entry.Equiv.Digest())
	}
	return nil
}

func describeKey(k meta.KeyRef) string {
	if k.Empty() {
		return "no key"
	}
	return k.Digest()
}

// fetchToStore downloads into a store transaction's staging path and publishes
// it through the ordinary destination rule: the flat name when it is free, the
// store when it is not.
//
// Commit is what verifies the payload against the lock, from the bytes rather
// than from what the registry claimed about them.
func fetchToStore(ctx context.Context, entry *lock.Entry, step Step, base, repo string,
	desc ocispec.Descriptor, annotations map[string]string) (fetched, error) {

	tx, err := store.Begin(entry.Manifest.Name, entry.Identity,
		store.BeginOptions{Equiv: entry.Equiv, StoreOnly: step.StoreOnly})
	if err != nil {
		return fetched{}, err
	}
	if !tx.Reserved() {
		// Somebody published the exact identity between planning and here.
		candidate, err := tx.Commit()
		if err != nil {
			return fetched{}, err
		}
		return fetched{Path: candidate.Path, Layout: candidate.Layout, Identity: candidate.Identity}, nil
	}
	if err := registry.Fetch(ctx, base, repo, desc, annotations, tx.Prepared, registry.KindOverlay); err != nil {
		tx.Abort()
		return fetched{}, err
	}
	candidate, err := tx.Commit()
	if err != nil {
		return fetched{}, err
	}
	return fetched{Path: candidate.Path, Layout: candidate.Layout, Identity: candidate.Identity}, nil
}

// fetchToDestination downloads a project-path artifact beside where it belongs
// and renames it into place.
//
// It never reaches the images roots in either direction. A copy of the same
// identity sitting in flat or store/ is not a substitute, because the
// destination is where the project *mounts* it from — adopting one would report
// a successful restore while the declared path stayed empty.
func fetchToDestination(ctx context.Context, root string, entry *lock.Entry, step Step, match Match,
	base, repo string, desc ocispec.Descriptor, annotations map[string]string) (fetched, error) {

	destination := filepath.Join(root, filepath.FromSlash(step.Destination))
	if err := utils.MkdirAllShared(filepath.Dir(destination)); err != nil {
		return fetched{}, err
	}
	guard, err := producer.AcquireLocal(destination)
	if err != nil {
		return fetched{}, err
	}
	defer guard.Release() //nolint:errcheck

	staged := producer.PreparedPath(destination, guard.Info())
	_ = os.Remove(staged)
	if err := registry.Fetch(ctx, base, repo, desc, annotations, staged, registry.KindOverlay); err != nil {
		os.Remove(staged) //nolint:errcheck
		return fetched{}, err
	}
	// From the payload, not from the annotations that got us here: a publisher
	// whose manifest disagrees with its own bytes must fail before the project
	// mounts them.
	verdict, diffs, err := verify(entry, staged, match)
	if err != nil {
		os.Remove(staged) //nolint:errcheck
		return fetched{}, fmt.Errorf("cannot verify the fetched artifact: %w", err)
	}
	if verdict != compare.Exact && verdict != compare.Equivalent {
		os.Remove(staged) //nolint:errcheck
		return fetched{}, fmt.Errorf("the fetched artifact is %s: %s", verdict, strings.Join(diffs, "; "))
	}
	path, err := placeAt(ctx, destination, staged)
	if err != nil {
		os.Remove(staged) //nolint:errcheck
		return fetched{}, err
	}
	return fetched{Path: path}, nil
}

// install publishes a store-addressed artifact through the destination rule:
// the flat name when it is free, the store when it is not — or when the plan
// has already given that name to a pin.
func install(entry *lock.Entry, step Step, output string) (store.Candidate, error) {
	return store.InstallFile(entry.Manifest.Name, entry.Identity, output,
		store.BeginOptions{Equiv: entry.Equiv, StoreOnly: step.StoreOnly})
}

// placeAt publishes a project destination by rename, refusing to overwrite
// anything that is not a plain file this restore may replace.
//
// A symlink is refused rather than followed: it points somewhere the project
// does not describe, and replacing what it points at would write outside the
// checkout.
func placeAt(ctx context.Context, destination, output string) (string, error) {
	if info, err := os.Lstat(destination); err == nil {
		switch {
		case info.Mode()&os.ModeSymlink != 0:
			return "", fmt.Errorf("%s is a symlink", destination)
		case !info.Mode().IsRegular():
			return "", fmt.Errorf("%s is not a regular file", destination)
		case !strings.HasSuffix(destination, ".sqf"):
			return "", fmt.Errorf("%s is not a .sqf", destination)
		}
		// The file is replaced, not merged: say so rather than losing someone's
		// hand-placed overlay silently. Not a refusal — a drifted project could
		// then be repaired only by deleting the file first.
		logging.FromContext(ctx).Warn("replacing the file already at this project path",
			"path", destination)
	} else if !os.IsNotExist(err) {
		return "", err
	}
	if err := utils.MkdirAllShared(filepath.Dir(destination)); err != nil {
		return "", err
	}
	if err := os.Rename(output, destination); err != nil {
		return "", fmt.Errorf("cannot publish %s: %w", destination, err)
	}
	utils.ShareWithParentGroup(destination)
	return destination, nil
}
