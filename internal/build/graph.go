package build

import (
	"context"
	"fmt"
	"os"
	"strings"
	"time"

	"github.com/Justype/condatainer/catalog"
	"github.com/Justype/condatainer/internal/config"
	"github.com/Justype/condatainer/internal/logging"
	"github.com/Justype/condatainer/internal/scheduler"
)

// BuildGraph turns a solved dependency plan into build work.
// The catalog resolves the graph and detects cycles; this type creates a
// BuildObject per missing node, in dependency-first order,
// and separates items into two ordered lists:
//   - localBuilds: BuildObjects without scheduler requirements (build locally)
//   - schedulerBuilds: BuildObjects that require scheduler submission
type BuildGraph struct {
	graph           map[string]*BuildObject // All build objects by name/version
	localBuilds     []*BuildObject          // Builds to run locally (no scheduler)
	schedulerBuilds []*BuildObject          // Builds to submit via scheduler
	jobIDs          map[string]string       // Job IDs for scheduler builds (name/version -> job ID)
	scheduler       scheduler.Scheduler     // Active scheduler (SLURM, PBS, etc.)

	// Config
	ctx        context.Context
	imagesDir  string
	tmpDir     string
	submitJobs bool // Whether to actually submit scheduler jobs
	update     bool // If true, rebuild all overlays even if already installed
}

// NewBuildGraph creates a BuildGraph from a list of BuildObjects
// All overlays are stored in imagesDir regardless of type
func NewBuildGraph(ctx context.Context, buildObjects []*BuildObject, imagesDir, tmpDir string, submitJobs bool, update bool) (*BuildGraph, error) {
	bg := &BuildGraph{
		graph:           make(map[string]*BuildObject),
		localBuilds:     []*BuildObject{},
		schedulerBuilds: []*BuildObject{},
		jobIDs:          make(map[string]string),
		ctx:             ctx,
		imagesDir:       imagesDir,
		tmpDir:          tmpDir,
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
		if update {
			// A root being rebuilt must resolve as missing, or the walk stops
			// at it and never reaches what it needs.
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

// installedVersions reports the versions of name already built, as the Have the
// catalog resolver asks with. Names carry slashes, so the version is the last
// segment and nothing deeper counts.
//
// hidden drops entries the graph intends to rebuild, which is how --update
// reaches past a root that is technically installed.
func installedVersions(hidden map[string]bool) catalog.Have {
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
// splits it into local and scheduler work.
//
// Resolution runs over the catalog index alone, so a cycle or an unresolvable
// name fails before any recipe is fetched or any temp file written, and an
// already-installed dependency costs a map lookup instead of a build object.
func (bg *BuildGraph) resolvePlan(ctx context.Context, roots []string, hidden map[string]bool) error {
	cat, err := config.OpenCatalog(ctx)
	if err != nil {
		return err
	}
	plan, err := cat.Resolve(ctx, roots, installedVersions(hidden))
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
			obj, err = NewBuildObject(ctx, name, false, bg.imagesDir, bg.tmpDir, false)
			if err != nil {
				return fmt.Errorf("failed to create BuildObject for dependency '%s': %w", name, err)
			}
			bg.graph[name] = obj
		}
		order = append(order, obj)
	}

	for _, obj := range order {
		if bg.submitJobs && bg.scheduler != nil && obj.RequiresScheduler() {
			bg.schedulerBuilds = append(bg.schedulerBuilds, obj)
		} else {
			bg.localBuilds = append(bg.localBuilds, obj)
		}
	}
	return nil
}

// Run executes the build graph
// First runs local builds, then submits scheduler jobs
func (bg *BuildGraph) Run(ctx context.Context) error {
	if err := bg.runLocalStep(ctx); err != nil {
		return err
	}
	if err := bg.runSchedulerStep(); err != nil {
		return err
	}

	// Check if any apptainer jobs were run
	hasDefBuilds := false
	for _, obj := range bg.schedulerBuilds {
		if obj.Type() == BuildTypeDef {
			hasDefBuilds = true
			break
		}
	}
	if !hasDefBuilds {
		for _, obj := range bg.localBuilds {
			if obj.Type() == BuildTypeDef {
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

// runLocalStep executes builds that don't require scheduler
func (bg *BuildGraph) runLocalStep(ctx context.Context) error {
	for _, meta := range bg.localBuilds {
		if !bg.update && meta.IsInstalled() {
			continue
		}
		logging.FromContext(ctx).Debug("processing overlay (local build)", "name", meta.NameVersion())
		if err := meta.Build(ctx, false); err != nil {
			return fmt.Errorf("failed to build %s: %w", meta.NameVersion(), err)
		}
	}
	return nil
}

// runSchedulerStep submits builds that require scheduler
func (bg *BuildGraph) runSchedulerStep() error {
	if bg.scheduler == nil {
		return nil
	}

	for _, meta := range bg.schedulerBuilds {
		if !bg.update && meta.IsInstalled() {
			continue
		}
		logging.FromContext(bg.ctx).Debug("processing overlay (scheduler job)", "name", meta.NameVersion())

		// Collect dependency job IDs
		depIDs := []string{}
		for _, rawDep := range meta.Dependencies() {
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
						dep, meta.NameVersion())
				}
				// Dependency is installed, no need to add to depIDs
			}
		}

		// Submit scheduler job
		jobID, err := bg.submitJob(meta, depIDs)
		if err != nil {
			return fmt.Errorf("failed to submit job for %s: %w", meta.NameVersion(), err)
		}

		bg.jobIDs[meta.NameVersion()] = jobID
	}
	return nil
}

// submitJob creates and submits a scheduler job for the build
func (bg *BuildGraph) submitJob(meta *BuildObject, depIDs []string) (string, error) {
	log := logging.FromContext(bg.ctx)
	log.Debug("submitting scheduler job", "type", bg.scheduler.GetType(), "name", meta.NameVersion(), "deps", depIDs)

	// Acquire lock before submitting to prevent duplicate scheduler submissions.
	// The lock is created with an empty job_id and updated after Submit() returns.
	lockPath := meta.LockPath()
	pendingLock := BuildLockInfo{
		Type:      string(bg.scheduler.GetType()), // e.g. "slurm", "pbs", "lsf", "htcondor"
		CreatedAt: time.Now().Format(time.RFC3339),
	}
	if err := acquireBuildLockFile(lockPath, pendingLock); err != nil {
		if os.IsExist(err) {
			return "", fmt.Errorf("build already queued or running for %s (lock exists at %s)",
				meta.NameVersion(), lockPath)
		}
		return "", fmt.Errorf("failed to create build lock for %s: %w", meta.NameVersion(), err)
	}

	// Get script specs; when always_submit forces submission without directives, synthesize empty specs
	specs := meta.ScriptSpecs()
	if specs == nil {
		effRS := buildEffectiveResourceSpec(nil)
		specs = &scheduler.ScriptSpecs{Spec: effRS}
	}
	if config.Global.ProxyPerJob {
		if h, err := os.Hostname(); err == nil && h != "" {
			specs.ProxyVia = h
		}
	}

	// Derive job name from name/version if not set in script
	if specs.Control.JobName == "" {
		name := meta.NameVersion()
		if idx := strings.LastIndex(name, "/"); idx != -1 {
			name = name[:idx]
		}
		specs.Control.JobName = "cnt-" + name
	}

	// Create job specification
	jobSpec := &scheduler.JobSpec{
		Name:           meta.NameVersion(),
		Command:        buildSchedulerCreateCommand(meta.NameVersion(), bg.update, meta.InteractiveInputs()),
		Specs:          specs,
		DepJobIDs:      depIDs,
		OverrideOutput: true,
		Metadata: map[string]string{
			"Target": meta.NameVersion(),
			"Type":   meta.Type().String(),
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
		return "", fmt.Errorf("failed to submit job: %w", err)
	}

	// Update lock with the actual job ID now that we have it.
	pendingLock.JobID = jobID
	_ = overwriteBuildLockFile(lockPath, pendingLock) // best-effort; we already hold the lock

	log.Info("submitted scheduler job", "type", bg.scheduler.GetType(), "jobID", jobID, "name", meta.NameVersion())
	return jobID, nil
}

// buildSchedulerCreateCommand returns the condatainer create command for scheduler jobs,
// propagating the --update flag if active.
// If interactiveInputs is non-empty, the inputs are embedded as a heredoc so the
// scheduler node does not need a TTY.
func buildSchedulerCreateCommand(nameVersion string, update bool, interactiveInputs []string) string {
	var cmd strings.Builder
	cmd.WriteString("condatainer create")
	if update {
		cmd.WriteString(" --update")
	}
	cmd.WriteString(" ")
	cmd.WriteString(nameVersion)
	if len(interactiveInputs) > 0 {
		cmd.WriteString(" << 'CNT_INPUTS_EOF'")
		for _, input := range interactiveInputs {
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
