package build

import (
	"context"
	"fmt"
	"os"
	"regexp"
	"strings"
	"time"

	"github.com/Justype/condatainer/catalog"
	"github.com/Justype/condatainer/internal/config"
	"github.com/Justype/condatainer/internal/logging"
	"github.com/Justype/condatainer/internal/scheduler"
)

// BuildGraph turns a solved dependency plan into build work: a BuildObject per
// missing node, in dependency-first order, split into localBuilds and
// schedulerBuilds.
type BuildGraph struct {
	graph           map[string]*BuildObject // All build objects by name/version
	localBuilds     []*BuildObject          // Builds to run locally (no scheduler)
	schedulerBuilds []*BuildObject          // Builds to submit via scheduler
	jobIDs          map[string]string       // Job IDs for scheduler builds (name/version -> job ID)
	scheduler       scheduler.Scheduler     // Active scheduler (SLURM, PBS, etc.)

	// Config
	ctx        context.Context
	imagesDir  string
	submitJobs bool // Whether to actually submit scheduler jobs
	update     bool // If true, rebuild all overlays even if already installed

	// jobFlags are the create flags that change what a build produces or where it
	// lands, repeated on every submitted job's command line.
	jobFlags []string
}

// SetJobFlags sets the flags every submitted job's create command carries, as
// separate arguments.
func (bg *BuildGraph) SetJobFlags(flags []string) { bg.jobFlags = flags }

// NewBuildGraph creates a BuildGraph from a list of BuildObjects
// All overlays are stored in imagesDir regardless of type
func NewBuildGraph(ctx context.Context, buildObjects []*BuildObject, imagesDir string, submitJobs bool, update bool) (*BuildGraph, error) {
	bg := &BuildGraph{
		graph:           make(map[string]*BuildObject),
		localBuilds:     []*BuildObject{},
		schedulerBuilds: []*BuildObject{},
		jobIDs:          make(map[string]string),
		ctx:             ctx,
		imagesDir:       imagesDir,
		submitJobs:      submitJobs,
		update:          update,
	}

	log := logging.FromContext(ctx)

	// Assign scheduler if job submission is enabled and we're not inside a job
	if submitJobs {
		if scheduler.IsInsideJob() {
			log.Info("already inside a scheduler job, all builds will run locally", "kind", "note")
		} else if sched := scheduler.ActiveScheduler(); sched != nil {
			bg.scheduler = sched
			log.Debug("using scheduler", "type", sched.GetType(), "binary", sched.GetBinary())
		} else {
			log.Warn("no scheduler detected, all builds will run locally")
		}
	}

	// Seed graph with provided build objects
	roots := make([]string, 0, len(buildObjects))
	hidden := map[string]bool{}
	for _, obj := range buildObjects {
		bg.graph[obj.NameVersion()] = obj
		roots = append(roots, obj.NameVersion())
		if update || obj.StoreOverflow() {
			// A root being rebuilt, or filed beside the installed name, must
			// resolve as missing, or the walk stops at it and never reaches what
			// it needs.
			hidden[obj.NameVersion()] = true
		} else if obj.IsInstalled() {
			log.Info("overlay already installed, skipping", "name", obj.NameVersion())
		}
	}

	if err := bg.resolvePlan(ctx, roots, hidden); err != nil {
		return nil, err
	}

	return bg, nil
}

// InstalledVersions reports the versions of name already built, as the Have the
// catalog resolver — and catalog.SolveName — asks with. The version is the
// last segment and nothing deeper counts. hidden drops entries the graph
// intends to rebuild, as --update does; pass nil outside a build graph.
func InstalledVersions(hidden map[string]bool) catalog.Have {
	return func(name string) []string {
		prefix := name + "/"
		var out []string
		for key := range getInstalledOverlays() {
			version, ok := strings.CutPrefix(key, prefix)
			if !ok || version == "" || strings.Contains(version, "/") || hidden[key] {
				continue
			}
			out = append(out, version)
		}
		return out
	}
}

// resolvePlan expands the seeded roots into a dependency-first build order and
// splits it into local and scheduler work, noting which each build gets. Resolution runs over the catalog
// index alone, before any recipe is fetched or temp file written.
func (bg *BuildGraph) resolvePlan(ctx context.Context, roots []string, hidden map[string]bool) error {
	cat, err := config.OpenCatalog(ctx)
	if err != nil {
		return err
	}
	plan, err := cat.Resolve(ctx, roots, InstalledVersions(hidden))
	if err != nil {
		return err
	}

	order := make([]*BuildObject, 0, len(plan.Order))
	for _, node := range plan.Order {
		name := node.Name()
		obj, seeded := bg.graph[name]
		if !seeded {
			// Installed and not a root: nothing to build. It is mounted into
			// the dependent's build, which is the only thing a #DEP: does.
			if node.Installed != "" {
				continue
			}
			obj, err = NewBuildObject(ctx, name, false, bg.imagesDir, false)
			if err != nil {
				return fmt.Errorf("failed to create BuildObject for dependency '%s': %w", name, err)
			}
			bg.graph[name] = obj
		}
		order = append(order, obj)
	}

	// Settled before the split below: a node with a prebuilt to pull is local
	// work, whatever its recipe would ask of a scheduler.
	for _, obj := range order {
		if bg.installedAlready(obj) {
			continue
		}
		planned, known := bg.dependencyPlan(obj)
		obj.plannedDeps = planned
		if err := obj.planPrebuilt(ctx, known); err != nil {
			return err
		}
	}
	order = bg.pruneDependencies(order, roots)

	log := logging.FromContext(ctx)
	for _, obj := range order {
		pull := obj.prebuilt.choice == prebuiltPull
		scheduled := !pull && bg.submitJobs && bg.scheduler != nil && obj.RequiresScheduler()
		if scheduled {
			bg.schedulerBuilds = append(bg.schedulerBuilds, obj)
		} else {
			bg.localBuilds = append(bg.localBuilds, obj)
		}
		if obj.IsInstalled() {
			continue
		}
		verb := "installed"
		if obj.prebuilt.choice == prebuiltNone {
			verb = "built"
		}
		switch {
		case pull:
			log.Info(obj.NameVersion()+" will be pulled from "+obj.prebuilt.candidate.endpoint, "kind", "note")
			if len(obj.prunedDeps) > 0 {
				log.Info(fmt.Sprintf("its build dependencies are not installed: %s; `condatainer install` adds them",
					strings.Join(obj.prunedDeps, ", ")), "kind", "note")
			}
		case scheduled:
			log.Info(obj.NameVersion()+" will be "+verb+" by a scheduler job", "kind", "note")
		default:
			log.Info(obj.NameVersion()+" will be "+verb+" locally", "kind", "note")
		}
	}
	return nil
}

// dependencyPlan returns what obj's dependencies will contribute to its
// equivalence, for those not installed yet, and whether every dependency's
// contribution is known: installed with keys, or planned. A dependency neither
// installed nor planned leaves obj undecided.
func (bg *BuildGraph) dependencyPlan(obj *BuildObject) (map[string]plannedDep, bool) {
	planned := make(map[string]plannedDep)
	for _, raw := range obj.Dependencies() {
		requested, constrained := raw, raw
		if !catalog.IsPathDep(raw) {
			parsed, err := catalog.ParseDep(raw)
			if err != nil {
				continue
			}
			requested, constrained = parsed.NameVersion(), parsed.String()
		}
		if depObj, inGraph := bg.graph[requested]; inGraph && !bg.installedAlready(depObj) {
			dep, ok := depObj.plannedAs()
			if !ok {
				return nil, false
			}
			planned[requested] = dep
			continue
		}
		if _, _, err := readDependencyManifest(constrained); err != nil {
			return nil, false
		}
	}
	return planned, true
}

// pruneDependencies drops the nodes only a pulled node depended on. A pull opens
// none of its build dependencies, so installing them is left to the user; a
// dependency a root or a node that builds needs stays.
func (bg *BuildGraph) pruneDependencies(order []*BuildObject, roots []string) []*BuildObject {
	needed := make(map[string]bool, len(order))
	for _, name := range roots {
		needed[name] = true
	}
	for i := len(order) - 1; i >= 0; i-- {
		obj := order[i]
		if !needed[obj.NameVersion()] || obj.prebuilt.choice == prebuiltPull {
			continue
		}
		for _, dep := range bg.graphDependencies(obj) {
			needed[dep] = true
		}
	}
	kept := make([]*BuildObject, 0, len(order))
	for _, obj := range order {
		if !needed[obj.NameVersion()] {
			obj.pruned = true
			continue
		}
		if obj.prebuilt.choice == prebuiltPull {
			for _, dep := range bg.graphDependencies(obj) {
				if !needed[dep] {
					obj.prunedDeps = append(obj.prunedDeps, dep)
				}
			}
		}
		kept = append(kept, obj)
	}
	return kept
}

// graphDependencies names the dependencies of obj that are nodes of this graph.
func (bg *BuildGraph) graphDependencies(obj *BuildObject) []string {
	var out []string
	for _, raw := range obj.Dependencies() {
		name := raw
		if parsed, err := catalog.ParseDep(raw); err == nil {
			name = parsed.NameVersion()
		}
		if depObj, ok := bg.graph[name]; ok && !bg.installedAlready(depObj) {
			out = append(out, name)
		}
	}
	return out
}

// resolveBase satisfies the implicit edge from every script and Conda build to
// the base image it runs inside, building it first if missing. Runs before any
// node; a definition-only plan resolves nothing here.
func (bg *BuildGraph) resolveBase(ctx context.Context) error {
	var dependents []*BuildObject
	for _, obj := range bg.graph {
		// A pull, and a node dropped for one, never runs in the base; a pull
		// that falls back to a build resolves it then.
		if obj.pruned || obj.prebuilt.choice == prebuiltPull {
			continue
		}
		if obj.BuildType() != BuildTypeDef && (bg.update || !obj.IsInstalled()) {
			dependents = append(dependents, obj)
		}
	}
	if len(dependents) == 0 {
		return nil
	}

	base, err := ResolveBase(ctx)
	if err != nil {
		return err
	}
	for _, obj := range dependents {
		obj.spec.Base = base
	}
	return nil
}

// Run executes the build graph
// First runs local builds, then submits scheduler jobs
func (bg *BuildGraph) Run(ctx context.Context) error {
	if err := bg.resolveBase(ctx); err != nil {
		return err
	}
	if err := bg.runLocalStep(ctx); err != nil {
		return err
	}
	if err := bg.runSchedulerStep(); err != nil {
		return err
	}

	// Check if any apptainer jobs were run
	hasDefBuilds := false
	for _, obj := range bg.schedulerBuilds {
		if obj.BuildType() == BuildTypeDef {
			hasDefBuilds = true
			break
		}
	}
	if !hasDefBuilds {
		for _, obj := range bg.localBuilds {
			if obj.BuildType() == BuildTypeDef {
				hasDefBuilds = true
				break
			}
		}
	}

	if hasDefBuilds {
		logging.FromContext(ctx).Debug("apptainer was used; run 'apptainer cache clean' to free up space")
	}

	return nil
}

// installedAlready reports whether a build has nothing to do because its name is
// installed. A store build is filed beside the installed name, so it never has.
func (bg *BuildGraph) installedAlready(obj *BuildObject) bool {
	return !bg.update && !obj.StoreOverflow() && obj.IsInstalled()
}

// runLocalStep executes builds that don't require scheduler
func (bg *BuildGraph) runLocalStep(ctx context.Context) error {
	for _, obj := range bg.localBuilds {
		if bg.installedAlready(obj) {
			continue
		}
		logging.FromContext(ctx).Debug("processing overlay (local build)", "name", obj.NameVersion())
		if err := obj.Build(ctx, len(obj.prunedDeps) > 0); err != nil {
			return fmt.Errorf("failed to build %s: %w", obj.NameVersion(), err)
		}
	}
	return nil
}

// runSchedulerStep submits builds that require scheduler
func (bg *BuildGraph) runSchedulerStep() error {
	if bg.scheduler == nil {
		return nil
	}

	for _, obj := range bg.schedulerBuilds {
		if bg.installedAlready(obj) {
			continue
		}
		logging.FromContext(bg.ctx).Debug("processing overlay (scheduler job)", "name", obj.NameVersion())

		// Collect dependency job IDs
		depIDs := []string{}
		for _, rawDep := range obj.Dependencies() {
			dep := rawDep
			if parsed, err := catalog.ParseDep(rawDep); err == nil {
				dep = parsed.NameVersion()
			}
			// Absent from the graph means the resolver satisfied it with an
			// installed version, so there is no job to wait on.
			if _, inGraph := bg.graph[dep]; !inGraph {
				continue
			}
			if jobID, exists := bg.jobIDs[dep]; exists {
				// Dependency was submitted as a scheduler job
				depIDs = append(depIDs, jobID)
			} else {
				// Check if dependency is already installed
				depObj, exists := bg.graph[dep]
				if !exists || !depObj.IsInstalled() {
					// Dependency should either be installed or have a job ID
					return fmt.Errorf("dependency %s for %s is not installed and was not submitted via scheduler",
						dep, obj.NameVersion())
				}
				// Dependency is installed, no need to add to depIDs
			}
		}

		// Submit scheduler job
		jobID, err := bg.submitJob(obj, depIDs)
		if err != nil {
			return fmt.Errorf("failed to submit job for %s: %w", obj.NameVersion(), err)
		}

		bg.jobIDs[obj.NameVersion()] = jobID
	}
	return nil
}

// submitJob creates and submits a scheduler job for the build, provisioning
// micromamba first when the build needs it.
func (bg *BuildGraph) submitJob(obj *BuildObject, depIDs []string) (string, error) {
	log := logging.FromContext(bg.ctx)
	log.Debug("submitting scheduler job", "type", bg.scheduler.GetType(), "name", obj.NameVersion(), "deps", depIDs)

	// The job runs the build on a compute node, which may have no outbound
	// access, so the toolchain is provisioned here on the submitting host.
	if obj.BuildType() != BuildTypeDef {
		if err := ensureMicromamba(bg.ctx); err != nil {
			return "", err
		}
	}

	// Acquire lock before submitting to prevent duplicate scheduler submissions.
	// The lock is created with an empty job_id and updated after Submit() returns.
	lockPath := obj.LockPath()
	pendingLock := BuildLockInfo{
		Runner:    string(bg.scheduler.GetType()), // e.g. "slurm", "pbs", "lsf", "htcondor"
		Node:      shortHostname(),
		PID:       os.Getpid(),
		CreatedAt: time.Now().Format(time.RFC3339),
	}
	if err := acquireBuildLockFile(lockPath, pendingLock); err != nil {
		if os.IsExist(err) {
			return "", fmt.Errorf("build already queued or running for %s (lock exists at %s)",
				obj.NameVersion(), lockPath)
		}
		return "", fmt.Errorf("failed to create build lock for %s: %w", obj.NameVersion(), err)
	}
	if err := obj.adoptWorkspace(pendingLock); err != nil {
		os.Remove(lockPath) // release the lock whose workspace could not be adopted
		return "", err
	}

	// Get script specs; when always_submit_data forces submission without directives, synthesize empty specs
	specs := obj.ScriptSpecs()
	if specs == nil {
		effRS := EffectiveResourceSpec(nil)
		specs = &scheduler.ScriptSpecs{Spec: effRS}
	}
	if config.Global.ProxyPerJob {
		if h, err := os.Hostname(); err == nil && h != "" {
			specs.ProxyVia = h
		}
	}

	// Derive job name from name/version if not set in script
	if specs.Control.JobName == "" {
		name := obj.NameVersion()
		if idx := strings.LastIndex(name, "/"); idx != -1 {
			name = name[:idx]
		}
		specs.Control.JobName = "cnt-" + name
	}

	// Create job specification
	jobSpec := &scheduler.JobSpec{
		Name:           obj.NameVersion(),
		Command:        buildSchedulerCreateCommand(obj.jobTarget(), bg.jobFlags, bg.update, obj.StoreOverflow(), config.Global.Build.SkipPrebuilt || obj.prebuilt.choice == prebuiltNone, obj.InputAnswers()),
		Specs:          specs,
		DepJobIDs:      depIDs,
		OverrideOutput: true,
		Metadata: map[string]string{
			"Target":     obj.NameVersion(),
			"Build Type": obj.BuildType().String(),
		},
	}

	// Create batch script
	scriptPath, err := bg.scheduler.CreateScriptWithSpec(jobSpec, config.Global.LogsDir)
	if err != nil {
		return "", fmt.Errorf("failed to create batch script: %w", err)
	}

	// Submit job (build chain always uses afterok)
	var deps []scheduler.Dependency
	if len(depIDs) > 0 {
		deps = []scheduler.Dependency{{Type: scheduler.DependencyAfterOK, JobIDs: depIDs}}
	}
	jobID, err := bg.scheduler.Submit(bg.ctx, scriptPath, deps)
	if err != nil {
		os.Remove(lockPath) // release lock on submission failure
		obj.Cleanup(true)   //nolint:errcheck
		return "", fmt.Errorf("failed to submit job: %w", err)
	}

	// Update lock with the actual job ID now that we have it.
	pendingLock.JobID = jobID
	_ = overwriteBuildLockFile(lockPath, pendingLock) // best-effort; we already hold the lock
	obj.Cleanup(false)                                //nolint:errcheck

	log.Info("submitted scheduler job", "type", bg.scheduler.GetType(), "jobID", jobID, "name", obj.NameVersion())
	return jobID, nil
}

// plainToken matches an argument that needs no quoting in the job script.
var plainToken = regexp.MustCompile(`^[A-Za-z0-9_./:=+@%-]+$`)

// jobToken renders one argument of the job's create command, quoted unless it is
// plainly safe: a name, a path or a flag.
func jobToken(arg string) string {
	if plainToken.MatchString(arg) {
		return arg
	}
	return shellQuote(arg)
}

// buildSchedulerCreateCommand returns the condatainer create command for a
// scheduler job, naming the build by target (a catalog name, or the arguments of
// an external script), and propagating flags, --update, --store and --no-prebuilt and
// embedding any input answers as a heredoc so the node needs no TTY.
func buildSchedulerCreateCommand(target, flags []string, update, store, noPrebuilt bool, inputAnswers []string) string {
	var cmd strings.Builder
	cmd.WriteString("condatainer create")
	for _, flag := range flags {
		cmd.WriteString(" ")
		cmd.WriteString(jobToken(flag))
	}
	if update {
		cmd.WriteString(" --update")
	}
	if store {
		cmd.WriteString(" --store")
	}
	if noPrebuilt {
		cmd.WriteString(" --no-prebuilt")
	}
	for _, token := range target {
		cmd.WriteString(" ")
		cmd.WriteString(jobToken(token))
	}
	if len(inputAnswers) > 0 {
		cmd.WriteString(" << 'CNT_INPUTS_EOF'")
		for _, input := range inputAnswers {
			cmd.WriteString("\n")
			cmd.WriteString(input)
		}
		cmd.WriteString("\nCNT_INPUTS_EOF")
	}
	return cmd.String()
}

// GetJobIDs returns the map of job IDs for scheduler builds
func (bg *BuildGraph) GetJobIDs() map[string]string {
	return bg.jobIDs
}
