package restore

import (
	"context"
	"errors"
	"fmt"
	"os"
	"path/filepath"
	"strings"
	"time"

	"github.com/Justype/condatainer/catalog"
	"github.com/Justype/condatainer/internal/artifact/meta"
	"github.com/Justype/condatainer/internal/build"
	"github.com/Justype/condatainer/internal/config"
	"github.com/Justype/condatainer/internal/image/producer"
	"github.com/Justype/condatainer/internal/logging"
	"github.com/Justype/condatainer/internal/project"
	"github.com/Justype/condatainer/internal/project/lock"
	"github.com/Justype/condatainer/internal/scheduler"
	"github.com/Justype/condatainer/internal/store"
	"github.com/Justype/condatainer/internal/utils"
)

// Submitted is one rebuild handed to the scheduler. The artifact is not
// available yet; a later restore adopts what the job publishes.
type Submitted struct {
	Artifact string `json:"artifact"`
	Name     string `json:"name"`
	Identity string `json:"identity"`
	JobID    string `json:"job_id"`
	// Destination is the project path the job will write, empty for a
	// store-addressed artifact.
	Destination string `json:"destination,omitempty"`
	// DependsOn are the job IDs this one waits for, sorted by step order.
	DependsOn []string `json:"depends_on,omitempty"`
	// Queued reports that the job was already in flight and this restore
	// submitted nothing. Re-running a restore is how it is resumed.
	Queued bool `json:"queued,omitempty"`
}

// jobs is one restore's scheduler partition: which steps the scheduler runs,
// which are produced inside somebody else's job, and the IDs to wait on.
type jobs struct {
	sched scheduler.Scheduler
	// submit and deferred are keyed by stepKey, and record decisions made
	// before any step ran.
	submit   map[string]bool
	deferred map[string]bool
	// specs are the parsed directives per artifact, for the steps that have any.
	specs map[string]*scheduler.ScriptSpecs
	// ids map a submitted artifact to its job ID, filled as submissions land.
	ids map[string]string
}

// stepKey identifies one step. A destination step is one output rather than one
// artifact, so two paths selecting one artifact are two jobs.
func stepKey(step Step) string { return step.Artifact + "\x00" + step.Destination }

// artifactKey is the stepKey of the step that produces an artifact for other
// steps to depend on. Only a store-addressed step is ever depended on, so it is
// the one with no destination.
func artifactKey(artifact string) string { return artifact + "\x00" }

// waitFor is the job IDs a step must run after: those of the dependencies this
// restore submitted, plus any it found already queued from an earlier one. A
// dependency already installed, or produced inside this step's own job, is
// nothing to wait on.
func (j *jobs) waitFor(step Step) []string {
	var ids []string
	for _, dependency := range step.DependsOn {
		if id, ok := j.ids[dependency]; ok && id != "" {
			ids = append(ids, id)
		}
	}
	return ids
}

// planJobs decides which of a plan's rebuilds the scheduler will run, before
// any of them starts.
//
// It reads each rebuild's directives out of its vendored recipe and hands the
// decision to partition. An empty result means everything runs here: submission
// was not asked for, no scheduler is installed, or this is already the job.
func planJobs(ctx context.Context, root string, verified *lock.Verified, plan *Plan, opts Options) *jobs {
	j := &jobs{
		submit:   map[string]bool{},
		deferred: map[string]bool{},
		specs:    map[string]*scheduler.ScriptSpecs{},
		ids:      map[string]string{},
	}
	if !opts.SubmitJobs {
		return j
	}
	log := logging.FromContext(ctx)
	if scheduler.IsInsideJob() {
		log.Info("already inside a scheduler job, restoring locally", "kind", "note")
		return j
	}
	sched := scheduler.ActiveScheduler()
	if sched == nil {
		return j
	}
	j.sched = sched

	// What each step declares for itself, and whether that alone calls for a job.
	needs := map[string]bool{}
	transient := map[string]bool{}
	for _, step := range plan.Steps {
		if step.Action != ActionBuild {
			continue
		}
		transient[stepKey(step)] = isTransient(step, opts)
		entry, ok := verified.Entries[step.Artifact]
		if !ok {
			continue
		}
		specs, err := recipeSpecs(root, entry)
		if err != nil {
			// A recipe that will not parse is the rebuild's problem to report,
			// not the partition's. Running it here surfaces the real error.
			log.Debug("cannot read scheduler directives", "name", step.Name, "err", err)
			continue
		}
		j.specs[step.Artifact] = specs
		needs[stepKey(step)] = callsForJob(entry.Manifest.Type, specs)
	}

	j.submit, j.deferred = partition(plan.Steps, needs, transient)
	return j
}

// callsForJob reports whether a rebuild of this type is a scheduler job by itself:
// its recipe carries directives, or it is data and build.always_submit_data is set.
func callsForJob(typ catalog.Type, specs *scheduler.ScriptSpecs) bool {
	return scheduler.HasSchedulerSpecs(specs) ||
		(config.Global.Build.AlwaysSubmitData && typ == catalog.TypeData)
}

// partition splits a plan into the steps the scheduler runs (submit) and the
// build dependencies a submitted job produces instead of this process (deferred).
// needs says which steps call for a job of their own (callsForJob).
//
// The split cannot be made step by step while executing. A build dependency comes
// before its dependent in the order, and whether it runs here or inside that
// dependent's job is only knowable once the dependent's fate is known.
//
// Three rules produce it:
//
//   - A build dependency is never its own job. It is produced inside the job of
//     whatever needs it, because the restore-scoped directory it lives in cannot
//     cross a job boundary. --keep-build-deps installs it instead, which makes it
//     an ordinary step and gives it its own job and its own directives.
//   - A step needs the scheduler when its own recipe carries directives, or when
//     a build dependency it will produce inline carries them. Otherwise a heavy
//     dependency would drag its light dependent onto the login node.
//   - A step whose dependency was submitted is submitted too, waiting on it. Its
//     input does not exist yet, so it cannot run here whatever it declares.
//
// It takes no scheduler, no checkout and no filesystem, which is the whole
// reason it is separate.
func partition(steps []Step, needs, transient map[string]bool) (submit, deferred map[string]bool) {
	submit, deferred = map[string]bool{}, map[string]bool{}

	// A build dependency's directives are its dependent's, since its dependent's
	// job is where it runs. Dependencies precede dependents in the order, so one
	// forward pass carries a whole inline chain up — which is why this does not
	// skip transient steps: a dependency two levels down reaches the top only by
	// being carried through the one between.
	for _, step := range steps {
		key := stepKey(step)
		for _, dependency := range step.DependsOn {
			if depKey := artifactKey(dependency); transient[depKey] && needs[depKey] {
				needs[key] = true
			}
		}
	}

	for _, step := range steps {
		key := stepKey(step)
		if step.Action != ActionBuild || transient[key] {
			continue
		}
		submit[key] = needs[key]
		for _, dependency := range step.DependsOn {
			if submit[artifactKey(dependency)] {
				submit[key] = true
			}
		}
	}

	// A build dependency is skipped here only when every dependent that would open
	// it runs in a job. One local dependent still needs it now, and the job
	// rebuilds its own copy.
	dependents := map[string][]string{}
	for _, step := range steps {
		if step.Action != ActionBuild {
			continue
		}
		for _, dependency := range step.DependsOn {
			dependents[dependency] = append(dependents[dependency], stepKey(step))
		}
	}
	for i := len(steps) - 1; i >= 0; i-- {
		step := steps[i]
		key := stepKey(step)
		if !transient[key] {
			continue
		}
		above := dependents[step.Artifact]
		if len(above) == 0 {
			continue
		}
		elsewhere := true
		for _, dependent := range above {
			if !submit[dependent] && !deferred[dependent] {
				elsewhere = false
				break
			}
		}
		deferred[key] = elsewhere
	}
	return submit, deferred
}

// recipeSpecs reads the scheduler directives a rebuild would run under.
//
// The vendored recipe is written to a scratch file because the parsers read
// paths rather than bytes. Nothing else about the build is set up, so this
// stays cheap enough to run over every step before the first one starts. A
// Conda build has no recipe and never carries directives.
func recipeSpecs(root string, entry *lock.Entry) (*scheduler.ScriptSpecs, error) {
	sources, err := lock.Sources(root, entry)
	if err != nil {
		return nil, err
	}
	text, ok := sources[meta.RecipeFileName]
	if !ok {
		return nil, nil
	}
	dir, err := os.MkdirTemp(utils.GetTmpDir(), "cnt-directives-")
	if err != nil {
		return nil, err
	}
	defer os.RemoveAll(dir) //nolint:errcheck
	path := filepath.Join(dir, meta.RecipeFileName)
	if err := os.WriteFile(path, text, utils.PermFile); err != nil {
		return nil, err
	}
	return scheduler.ReadScriptSpecsFromPath(path)
}

// submitStep hands one rebuild to the scheduler and reserves where it will land.
//
// The reservation is an ordinary producer lock carrying the job ID, so the job
// adopts its own lock when it installs and no second mechanism is needed to
// tell "queued" from "someone else is building this". A lock already held by a
// live job means this artifact was submitted by an earlier restore: that job is
// reported and nothing is submitted again, which is the whole of resume.
func (j *jobs) submitStep(ctx context.Context, root string, entry *lock.Entry, step Step, depIDs []string,
	opts Options) (*Submitted, *Result, *Failure) {

	fail := func(format string, args ...any) *Failure {
		return &Failure{Artifact: step.Artifact, Name: step.Name, Reason: fmt.Sprintf(format, args...)}
	}
	queued := func(info producer.Info) *Submitted {
		// Recorded like a fresh submission: a dependent of this step still has to
		// wait for the job that is already running it.
		j.ids[step.Artifact] = info.JobID
		return &Submitted{Artifact: step.Artifact, Name: step.Name, Identity: step.Identity,
			Destination: step.Destination, JobID: info.JobID, Queued: true}
	}

	reservation, adopted, err := reserve(root, entry, step)
	if err != nil {
		var producing *producer.ProducingError
		if errors.As(err, &producing) {
			return queued(producing.Info), nil, nil
		}
		return nil, nil, fail("cannot reserve where the job will publish: %v", err)
	}
	if adopted != nil {
		// It appeared between planning and now. Nothing to submit, and the
		// result is honest about where it came from.
		return nil, adopted, nil
	}

	jobID, err := j.launch(ctx, root, entry, step, depIDs, opts)
	if err != nil {
		reservation.release()
		return nil, nil, fail("%v", err)
	}
	info := producer.Info{
		Runner:    string(j.sched.GetType()),
		JobID:     jobID,
		Node:      producer.ShortHostname(),
		PID:       os.Getpid(),
		CreatedAt: time.Now().Format(time.RFC3339),
	}
	if err := reservation.detach(info); err != nil {
		// The job is already queued, so the reservation cannot be given back.
		// It will run and publish; the lock it could not claim only means a
		// concurrent producer is not warned off.
		logging.FromContext(ctx).Warn("could not hand the reservation to the job",
			"name", step.Name, "jobID", jobID, "err", err)
	}
	j.ids[step.Artifact] = jobID
	return &Submitted{Artifact: step.Artifact, Name: step.Name, Identity: step.Identity,
		Destination: step.Destination, JobID: jobID, DependsOn: depIDs}, nil, nil
}

// launch generates and submits the batch job that will re-enter restore for one
// artifact.
//
// The job re-runs `project restore --only`, so the worker recomputes the plan on
// the node and produces that artifact together with whatever build dependencies it
// needs. Nothing about the build is serialized into the job: the lock is the
// specification, and the node reads the same checkout.
func (j *jobs) launch(ctx context.Context, root string, entry *lock.Entry, step Step,
	depIDs []string, opts Options) (string, error) {

	// Copied because two project destinations selecting one artifact are two
	// jobs off one parse, and each fills in its own control settings.
	specs := &scheduler.ScriptSpecs{}
	if parsed := j.specs[step.Artifact]; parsed != nil {
		copied := *parsed
		specs = &copied
	}
	if scheduler.IsPassthrough(specs) {
		return "", fmt.Errorf("the recipe has scheduler directives that cannot be normalized")
	}
	// The same priority chain a catalog build resolves: build defaults, then the
	// recipe's own directives, then the resources of the job doing the submitting.
	specs.Spec = build.EffectiveResourceSpec(specs)
	// Parsed from a scratch copy of the vendored recipe, which is gone. The job
	// header echoes it, and a path that no longer exists is worse than none.
	specs.ScriptPath = ""
	workdir, err := project.WorkDir(root, specs.Control.WorkDir)
	if err != nil {
		return "", err
	}
	specs.Control.WorkDir = workdir
	if specs.Control.JobName == "" {
		name := entry.Manifest.Name
		if idx := strings.LastIndex(name, "/"); idx != -1 {
			name = name[:idx]
		}
		specs.Control.JobName = "cnt-" + name
	}
	if config.Global.ProxyPerJob {
		if host, err := os.Hostname(); err == nil && host != "" {
			specs.ProxyVia = host
		}
	}

	// A destination is part of the job's name so two of them do not write one
	// log, which OverrideOutput derives from it.
	jobName := entry.Manifest.Name
	if step.Destination != "" {
		jobName += "-" + strings.TrimSuffix(filepath.Base(step.Destination), ".sqf")
	}
	jobSpec := &scheduler.JobSpec{
		Name:           jobName,
		Command:        restoreCommand(root, step, opts),
		Specs:          specs,
		DepJobIDs:      depIDs,
		OverrideOutput: true,
		Metadata: map[string]string{
			"Target":   entry.Manifest.Name,
			"Identity": step.Identity,
			"Project":  root,
		},
	}
	scriptPath, err := j.sched.CreateScriptWithSpec(jobSpec, config.Global.LogsDir)
	if err != nil {
		return "", fmt.Errorf("cannot create the batch script: %w", err)
	}
	var deps []scheduler.Dependency
	if len(depIDs) > 0 {
		deps = []scheduler.Dependency{{Type: scheduler.DependencyAfterOK, JobIDs: depIDs}}
	}
	jobID, err := j.sched.Submit(ctx, scriptPath, deps)
	if err != nil {
		return "", fmt.Errorf("cannot submit: %w", err)
	}
	logging.FromContext(ctx).Info("submitted scheduler job", "kind", "note",
		"type", j.sched.GetType(), "jobID", jobID, "name", step.Name)
	return jobID, nil
}

// restoreCommand is what the job runs: the same restore, narrowed to one
// artifact and carrying the policy flags that decided this plan.
//
// --project is passed rather than relied on: a project is anchored on the
// working directory, and the job's is set from the root, but naming it keeps
// the job correct under a scheduler that resolves the directory differently.
func restoreCommand(root string, step Step, opts Options) string {
	var cmd strings.Builder
	cmd.WriteString("condatainer project restore --project ")
	cmd.WriteString(shellQuote(root))
	cmd.WriteString(" --only ")
	cmd.WriteString(shellQuote(step.Artifact))
	if opts.Match.Normalize() == MatchIdentity {
		cmd.WriteString(" --match identity")
	}
	if opts.SkipPrebuilt {
		cmd.WriteString(" --no-prebuilt")
	}
	if opts.KeepBuildDeps {
		cmd.WriteString(" --keep-build-deps")
	}
	return cmd.String()
}

// shellQuote wraps a value for the generated job script. A project root and a
// vendored artifact path are both user-chosen and may hold anything a
// filesystem allows.
func shellQuote(value string) string {
	return "'" + strings.ReplaceAll(value, "'", `'\''`) + "'"
}

// reservation is the claim a submitting process leaves on where a job will
// publish, whether that is the store or a project path.
type reservation struct {
	tx    *store.Transaction
	guard *producer.Guard
	path  string
}

func (r *reservation) release() {
	if r.tx != nil {
		r.tx.Abort()
	}
	if r.guard != nil {
		r.guard.Release() //nolint:errcheck
	}
}

func (r *reservation) detach(info producer.Info) error {
	if r.tx != nil {
		return r.tx.Detach(info)
	}
	r.guard = nil // the job owns the lock now; releasing it here would free the target
	return producer.Overwrite(producer.Path(r.path), info)
}

// reserve claims the pathname a submitted job will write, so a second restore
// reports the queued job instead of submitting a duplicate.
//
// It returns an adopted result when an exact copy turned up between planning and
// submission, which is a race worth winning rather than an error.
func reserve(root string, entry *lock.Entry, step Step) (*reservation, *Result, error) {
	if step.Destination != "" {
		target := filepath.Join(root, filepath.FromSlash(step.Destination))
		if err := utils.MkdirAllShared(filepath.Dir(target)); err != nil {
			return nil, nil, err
		}
		guard, err := producer.AcquireLocal(target)
		if err != nil {
			return nil, nil, err
		}
		return &reservation{guard: guard, path: target}, nil, nil
	}
	tx, err := store.Begin(entry.Manifest.Name, entry.Identity,
		store.BeginOptions{Equiv: entry.Equiv, StoreOnly: step.StoreOnly})
	if err != nil {
		return nil, nil, err
	}
	if !tx.Reserved() {
		candidate, err := tx.Commit()
		if err != nil {
			return nil, nil, err
		}
		return nil, &Result{
			Artifact: step.Artifact, Name: step.Name, Identity: step.Identity,
			Found: step.Found, Verdict: step.verdict(), Outcome: OutcomeAdopted,
			Path: candidate.Path, Layout: candidate.Layout,
		}, nil
	}
	return &reservation{tx: tx, path: tx.TargetPath}, nil, nil
}
