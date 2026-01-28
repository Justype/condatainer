package build

import (
	"fmt"
	"strings"

	"github.com/Justype/condatainer/internal/scheduler"
	"github.com/Justype/condatainer/internal/utils"
)

// BuildGraph manages dependency resolution and build ordering
// It expands BuildObject items to include all transitive dependencies,
// detects cycles, produces a topologically-sorted order (dependencies before dependents),
// and separates items into two ordered lists:
//   - localBuilds: BuildObjects without scheduler requirements (build locally)
//   - schedulerBuilds: BuildObjects that require scheduler submission
type BuildGraph struct {
	graph           map[string]BuildObject // All build objects by name/version
	localBuilds     []BuildObject          // Builds to run locally (no scheduler)
	schedulerBuilds []BuildObject          // Builds to submit via scheduler
	jobIDs          map[string]string      // Job IDs for scheduler builds (name/version -> job ID)
	scheduler       scheduler.Scheduler    // Active scheduler (SLURM, PBS, etc.)

	// Config
	imagesDir  string
	tmpDir     string
	submitJobs bool // Whether to actually submit scheduler jobs
}

// NewBuildGraph creates a BuildGraph from a list of BuildObjects
// All overlays are stored in imagesDir regardless of type
func NewBuildGraph(buildObjects []BuildObject, imagesDir, tmpDir string, submitJobs bool) (*BuildGraph, error) {
	bg := &BuildGraph{
		graph:           make(map[string]BuildObject),
		localBuilds:     []BuildObject{},
		schedulerBuilds: []BuildObject{},
		jobIDs:          make(map[string]string),
		imagesDir:       imagesDir,
		tmpDir:          tmpDir,
		submitJobs:      submitJobs,
	}

	// Detect and initialize scheduler if jobs should be submitted
	if submitJobs {
		sched, err := scheduler.DetectScheduler()
		if err != nil {
			utils.PrintWarning("No scheduler detected, all builds will run locally")
		} else if !sched.IsAvailable() {
			utils.PrintWarning("Scheduler not available, all builds will run locally")
		} else {
			bg.scheduler = sched
			info := sched.GetInfo()
			utils.PrintDebug("Using %s scheduler: %s", info.Type, info.Binary)
		}
	}

	// Seed graph with provided build objects
	for _, obj := range buildObjects {
		bg.graph[obj.NameVersion()] = obj
		if obj.IsInstalled() {
			utils.PrintMessage("Overlay %s is already installed. Skipping.", utils.StyleName(obj.NameVersion()))
		}
	}

	// Topologically sort the graph
	if err := bg.topologicalSort(); err != nil {
		return nil, err
	}

	return bg, nil
}

// topologicalSort performs topological sort with cycle detection
// and separates builds into local and sbatch lists
func (bg *BuildGraph) topologicalSort() error {
	visiting := make(map[string]bool)
	visited := make(map[string]bool)
	order := []string{}

	var visit func(string) error
	visit = func(node string) error {
		if visited[node] {
			return nil
		}
		if visiting[node] {
			return fmt.Errorf("circular dependency detected involving '%s'", node)
		}

		visiting[node] = true
		nodeMeta, exists := bg.graph[node]
		if !exists {
			// Expand on-the-fly if needed
			newObj, err := NewBuildObject(node, false, bg.imagesDir, bg.tmpDir)
			if err != nil {
				delete(visiting, node)
				return fmt.Errorf("failed to create BuildObject for dependency '%s': %w", node, err)
			}
			bg.graph[node] = newObj
			nodeMeta = newObj
		}

		// Visit dependencies first
		for _, dep := range nodeMeta.Dependencies() {
			if err := visit(dep); err != nil {
				return err
			}
		}

		delete(visiting, node)
		visited[node] = true
		order = append(order, node)
		return nil
	}

	// Visit all nodes in graph
	for node := range bg.graph {
		if !visited[node] {
			if err := visit(node); err != nil {
				return err
			}
		}
	}

	// Separate into local and scheduler builds based on sorted order
	for _, nameVersion := range order {
		meta := bg.graph[nameVersion]
		if bg.submitJobs && bg.scheduler != nil && meta.RequiresScheduler() {
			bg.schedulerBuilds = append(bg.schedulerBuilds, meta)
		} else {
			bg.localBuilds = append(bg.localBuilds, meta)
		}
	}

	return nil
}

// Run executes the build graph
// First runs local builds, then submits scheduler jobs
func (bg *BuildGraph) Run() error {
	if err := bg.runLocalStep(); err != nil {
		return err
	}
	if err := bg.runSchedulerStep(); err != nil {
		return err
	}

	// Check if any apptainer jobs were run
	hasDefBuilds := false
	for _, obj := range bg.schedulerBuilds {
		if obj.IsDef() {
			hasDefBuilds = true
			break
		}
	}
	if !hasDefBuilds {
		for _, obj := range bg.localBuilds {
			if obj.IsDef() {
				hasDefBuilds = true
				break
			}
		}
	}

	if hasDefBuilds {
		utils.PrintDebug("Apptainer is used. You can run 'apptainer cache clean' to free up space.")
	}

	return nil
}

// runLocalStep executes builds that don't require scheduler
func (bg *BuildGraph) runLocalStep() error {
	for _, meta := range bg.localBuilds {
		if meta.IsInstalled() {
			continue
		}
		utils.PrintDebug("Processing overlay %s (local build)...", meta.NameVersion())
		if err := meta.Build(false); err != nil {
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
		utils.PrintDebug("Processing overlay %s (scheduler job)...", meta.NameVersion())

		// Collect dependency job IDs
		depIDs := []string{}
		for _, dep := range meta.Dependencies() {
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
func (bg *BuildGraph) submitJob(meta BuildObject, depIDs []string) (string, error) {
	info := bg.scheduler.GetInfo()
	utils.PrintDebug("Submitting %s job for %s with dependencies: %s",
		info.Type, meta.NameVersion(), strings.Join(depIDs, ", "))

	// Create job specification
	jobSpec := &scheduler.JobSpec{
		Name:      meta.NameVersion(),
		Command:   fmt.Sprintf("condatainer create %s", meta.NameVersion()),
		Specs:     meta.ScriptSpecs(),
		DepJobIDs: depIDs,
	}

	// Create batch script
	scriptPath, err := bg.scheduler.CreateScriptWithSpec(jobSpec, bg.tmpDir)
	if err != nil {
		return "", fmt.Errorf("failed to create batch script: %w", err)
	}

	// Submit job
	jobID, err := bg.scheduler.Submit(scriptPath, depIDs)
	if err != nil {
		return "", fmt.Errorf("failed to submit job: %w", err)
	}

	utils.PrintMessage("Submitted %s job %s for %s", info.Type, jobID, meta.NameVersion())
	return jobID, nil
}

// GetLocalBuilds returns the list of builds to run locally
func (bg *BuildGraph) GetLocalBuilds() []BuildObject {
	return bg.localBuilds
}

// GetSchedulerBuilds returns the list of builds to submit via scheduler
func (bg *BuildGraph) GetSchedulerBuilds() []BuildObject {
	return bg.schedulerBuilds
}

// GetJobIDs returns the map of job IDs for scheduler builds
func (bg *BuildGraph) GetJobIDs() map[string]string {
	return bg.jobIDs
}

// GetAllBuilds returns all builds in topologically sorted order
func (bg *BuildGraph) GetAllBuilds() []BuildObject {
	all := make([]BuildObject, 0, len(bg.localBuilds)+len(bg.schedulerBuilds))
	all = append(all, bg.localBuilds...)
	all = append(all, bg.schedulerBuilds...)
	return all
}
