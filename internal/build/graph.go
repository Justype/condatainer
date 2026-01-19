package build

import (
	"fmt"
	"strings"

	"github.com/Justype/condatainer/internal/utils"
)

// BuildGraph manages dependency resolution and build ordering
// It expands BuildObject items to include all transitive dependencies,
// detects cycles, produces a topologically-sorted order (dependencies before dependents),
// and separates items into two ordered lists:
//   - localBuilds: BuildObjects with sbatch == false (build locally)
//   - sbatchBuilds: BuildObjects with sbatch == true (submit to SLURM)
type BuildGraph struct {
	graph        map[string]BuildObject   // All build objects by name/version
	localBuilds  []BuildObject            // Builds to run locally (no sbatch)
	sbatchBuilds []BuildObject            // Builds to submit via sbatch
	jobIDs       map[string]string        // Job IDs for sbatch builds (name/version -> job ID)

	// Config
	imagesDir  string
	tmpDir     string
	submitJobs bool // Whether to actually submit sbatch jobs
}

// NewBuildGraph creates a BuildGraph from a list of BuildObjects
// All overlays are stored in imagesDir regardless of type
func NewBuildGraph(buildObjects []BuildObject, imagesDir, tmpDir string, submitJobs bool) (*BuildGraph, error) {
	bg := &BuildGraph{
		graph:        make(map[string]BuildObject),
		localBuilds:  []BuildObject{},
		sbatchBuilds: []BuildObject{},
		jobIDs:       make(map[string]string),
		imagesDir:    imagesDir,
		tmpDir:       tmpDir,
		submitJobs:   submitJobs,
	}

	// Seed graph with provided build objects
	for _, obj := range buildObjects {
		bg.graph[obj.NameVersion()] = obj
		if obj.IsInstalled() {
			utils.PrintDebug("Overlay %s is already installed. Skipping.", obj.NameVersion())
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

	// Separate into local and sbatch builds based on sorted order
	for _, nameVersion := range order {
		meta := bg.graph[nameVersion]
		// TODO: Check if sbatch is available
		if bg.submitJobs && meta.NeedsSbatch() {
			bg.sbatchBuilds = append(bg.sbatchBuilds, meta)
		} else {
			bg.localBuilds = append(bg.localBuilds, meta)
		}
	}

	return nil
}

// Run executes the build graph
// First runs local builds, then submits sbatch jobs
func (bg *BuildGraph) Run() error {
	if err := bg.runLocalStep(); err != nil {
		return err
	}
	if err := bg.runSbatchStep(); err != nil {
		return err
	}

	// Check if any apptainer jobs were run
	hasDefBuilds := false
	for _, obj := range bg.sbatchBuilds {
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

// runLocalStep executes builds that don't require sbatch
func (bg *BuildGraph) runLocalStep() error {
	for _, meta := range bg.localBuilds {
		if meta.IsInstalled() {
			continue
		}
		utils.PrintDebug("Processing overlay %s (no sbatch)...", meta.NameVersion())
		if err := meta.Build(false); err != nil {
			return fmt.Errorf("failed to build %s: %w", meta.NameVersion(), err)
		}
	}
	return nil
}

// runSbatchStep submits builds that require sbatch
func (bg *BuildGraph) runSbatchStep() error {
	for _, meta := range bg.sbatchBuilds {
		utils.PrintDebug("Processing overlay %s (with sbatch)...", meta.NameVersion())

		// Collect dependency job IDs
		depIDs := []string{}
		for _, dep := range meta.Dependencies() {
			if jobID, exists := bg.jobIDs[dep]; exists {
				depIDs = append(depIDs, jobID)
			} else if !meta.IsInstalled() {
				// Dependency should either be installed or have a job ID
				return fmt.Errorf("dependency %s for %s is not installed and was not submitted via sbatch",
					dep, meta.NameVersion())
			}
		}

		// Submit sbatch job
		jobID, err := bg.submitSbatch(meta, depIDs)
		if err != nil {
			return fmt.Errorf("failed to submit sbatch for %s: %w", meta.NameVersion(), err)
		}

		bg.jobIDs[meta.NameVersion()] = jobID
	}
	return nil
}

// submitSbatch creates and submits a SLURM batch script for the build
func (bg *BuildGraph) submitSbatch(meta BuildObject, depIDs []string) (string, error) {
	// TODO: Implement actual sbatch script creation and submission
	// This will involve:
	// 1. Creating a temporary sbatch script with scheduler directives
	// 2. Adding the condatainer create command
	// 3. Submitting via sbatch with dependency flags
	// 4. Parsing the job ID from output
	// 5. Returning the job ID

	utils.PrintDebug("Submitting sbatch for %s with dependencies: %s",
		meta.NameVersion(), strings.Join(depIDs, ", "))

	// Placeholder - return error for now
	return "", fmt.Errorf("sbatch submission not yet implemented")
}

// GetLocalBuilds returns the list of builds to run locally
func (bg *BuildGraph) GetLocalBuilds() []BuildObject {
	return bg.localBuilds
}

// GetSbatchBuilds returns the list of builds to submit via sbatch
func (bg *BuildGraph) GetSbatchBuilds() []BuildObject {
	return bg.sbatchBuilds
}

// GetJobIDs returns the map of job IDs for sbatch builds
func (bg *BuildGraph) GetJobIDs() map[string]string {
	return bg.jobIDs
}

// GetAllBuilds returns all builds in topologically sorted order
func (bg *BuildGraph) GetAllBuilds() []BuildObject {
	all := make([]BuildObject, 0, len(bg.localBuilds)+len(bg.sbatchBuilds))
	all = append(all, bg.localBuilds...)
	all = append(all, bg.sbatchBuilds...)
	return all
}
