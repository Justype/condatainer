package build

import (
	"errors"
	"fmt"
	"os"
	"path/filepath"
	"runtime"
	"strings"

	"github.com/Justype/condatainer/internal/config"
	"github.com/Justype/condatainer/internal/overlay"
	"github.com/Justype/condatainer/internal/scheduler"
	"github.com/Justype/condatainer/internal/utils"
)

var ErrTmpOverlayExists = errors.New("temporary overlay already exists")
var ErrBuildCancelled = errors.New("build cancelled by user")

// ScriptSpecs mirrors the scheduler module's job spec metadata.
type ScriptSpecs = scheduler.ScriptSpecs

func defaultBuildNcpus() int {
	if config.Global.Build.DefaultCPUs > 0 {
		return config.Global.Build.DefaultCPUs
	}
	if count := runtime.NumCPU(); count > 0 {
		return count
	}
	return 1
}

// BuildObject represents a build target with its metadata and dependencies
type BuildObject interface {
	// Core properties
	NameVersion() string
	BuildSource() string
	Dependencies() []string
	IsInstalled() bool

	// Type checks (implemented by concrete types)
	IsConda() bool
	IsDef() bool
	IsShell() bool
	IsRef() bool

	// Path management
	TmpOverlayPath() string
	TargetOverlayPath() string
	CntDirPath() string

	// Scheduler metadata
	ScriptSpecs() *ScriptSpecs
	RequiresScheduler() bool

	// Build operations
	Build(buildDeps bool) error
	GetMissingDependencies() ([]string, error)
	CreateTmpOverlay(force bool) error
	Cleanup(failed bool) error

	// String representation
	String() string
}

// BaseBuildObject provides common functionality for all build types
type BaseBuildObject struct {
	nameVersion       string
	buildSource       string
	dependencies      []string
	tmpOverlayPath    string
	targetOverlayPath string
	cntDirPath        string
	ncpus             int
	isRemote          bool // Whether build source was downloaded
	scriptSpecs       *scheduler.ScriptSpecs

	// Interactive inputs for shell scripts
	interactiveInputs []string
}

// Common interface implementations for BaseBuildObject

func (b *BaseBuildObject) NameVersion() string       { return b.nameVersion }
func (b *BaseBuildObject) BuildSource() string       { return b.buildSource }
func (b *BaseBuildObject) Dependencies() []string    { return b.dependencies }
func (b *BaseBuildObject) TmpOverlayPath() string    { return b.tmpOverlayPath }
func (b *BaseBuildObject) TargetOverlayPath() string { return b.targetOverlayPath }
func (b *BaseBuildObject) CntDirPath() string        { return b.cntDirPath }
func (b *BaseBuildObject) ScriptSpecs() *ScriptSpecs { return b.scriptSpecs }
func (b *BaseBuildObject) RequiresScheduler() bool {
	if specs := b.scriptSpecs; specs != nil {
		return len(specs.RawFlags) > 0
	}
	return false
}

func (b *BaseBuildObject) IsInstalled() bool {
	// Check if target overlay exists
	_, err := os.Stat(b.targetOverlayPath)
	return err == nil
}

func (b *BaseBuildObject) GetMissingDependencies() ([]string, error) {
	// TODO: Implement check for which dependencies are not installed
	missing := []string{}
	for _, dep := range b.dependencies {
		// For now, just return empty - needs integration with overlay checking
		_ = dep
	}
	return missing, nil
}

func (b *BaseBuildObject) CreateTmpOverlay(force bool) error {
	// Check if tmp overlay already exists
	if _, err := os.Stat(b.tmpOverlayPath); err == nil {
		if !force {
			return fmt.Errorf("%w: %s", ErrTmpOverlayExists, b.tmpOverlayPath)
		}
		// Remove existing if force=true
		if err := os.Remove(b.tmpOverlayPath); err != nil {
			return fmt.Errorf("failed to remove existing tmp overlay: %w", err)
		}
	}

	// Create actual ext3 overlay using overlay package
	utils.PrintDebug("Creating temporary overlay at %s", utils.StylePath(b.tmpOverlayPath))

	// Use overlay package to create ext3 overlay
	// For build overlays, we use a temporary size from config (default 20GB)
	// Use "conda" profile for small files, sparse=true for faster creation
	if err := overlay.CreateForCurrentUser(b.tmpOverlayPath, config.Global.Build.TmpSizeMB, "default", true, config.Global.Build.OverlayType); err != nil {
		return fmt.Errorf("failed to create temporary overlay: %w", err)
	}

	return nil
}

func (b *BaseBuildObject) Cleanup(failed bool) error {
	// Remove remote build source if downloaded
	if b.isRemote && b.buildSource != "" {
		if err := os.Remove(b.buildSource); err != nil && !os.IsNotExist(err) {
			utils.PrintDebug("Failed to remove remote build source %s: %v", b.buildSource, err)
		}
	}

	// Remove tmp overlay
	if b.tmpOverlayPath != "" {
		if err := os.Remove(b.tmpOverlayPath); err != nil && !os.IsNotExist(err) {
			utils.PrintDebug("Failed to remove tmp overlay %s: %v", b.tmpOverlayPath, err)
		}
	}

	// Remove cnt directory
	if b.cntDirPath != "" {
		cntBaseDir := filepath.Dir(b.cntDirPath)
		if err := os.RemoveAll(cntBaseDir); err != nil && !os.IsNotExist(err) {
			utils.PrintDebug("Failed to remove cnt dir %s: %v", cntBaseDir, err)
		}
	}

	// If failed, also remove target overlay
	if failed && b.targetOverlayPath != "" {
		if err := os.Remove(b.targetOverlayPath); err != nil && !os.IsNotExist(err) {
			utils.PrintDebug("Failed to remove target overlay %s: %v", b.targetOverlayPath, err)
		}
	}

	return nil
}

// parseScriptMetadata extracts dependencies, sbatch flags, and interactive prompts from shell scripts
func (b *BaseBuildObject) parseScriptMetadata() error {
	if b.buildSource == "" {
		return nil
	}

	if !config.Global.SubmitJob {
		return nil
	}

	sched := scheduler.ActiveScheduler()
	if sched == nil {
		return nil
	}

	b.dependencies = []string{}
	b.interactiveInputs = []string{}

	specs, err := sched.ReadScriptSpecs(b.buildSource)
	if err != nil {
		return err
	}

	if specs.Ncpus <= 0 {
		specs.Ncpus = defaultBuildNcpus()
	}

	b.scriptSpecs = specs
	b.ncpus = specs.Ncpus

	return nil
}

// Helper functions

// getCntDirPath returns the container directory path for a name/version
// Format: <tmpDir>/build_<nameVersion>/cnt
func getCntDirPath(nameVersion, tmpDir string) string {
	buildDirName := "build_" + strings.ReplaceAll(nameVersion, "/", "_")
	return filepath.Join(tmpDir, buildDirName, "cnt")
}

// getTmpOverlayPath returns the temporary overlay path
// Format: <tmpDir>/<nameVersion>.sqf (with / replaced by --)
func getTmpOverlayPath(nameVersion, tmpDir string) string {
	filename := strings.ReplaceAll(nameVersion, "/", "--") + ".img"
	return filepath.Join(tmpDir, filename)
}

// NewBuildObject creates a BuildObject from a name/version string
// Format: "name/version" for conda/shell, "name" for def, "prefix/name/version" for ref
// All overlays are stored in imagesDir regardless of type
func NewBuildObject(nameVersion string, external bool, imagesDir, tmpDir string) (BuildObject, error) {
	normalized := utils.NormalizeNameVersion(nameVersion)
	slashCount := strings.Count(normalized, "/")

	// Determine build type based on slash count
	// 0 slashes: def file (e.g., "ubuntu")
	// 1 slash: conda or shell (e.g., "numpy/1.24")
	// 2+ slashes: ref shell script (e.g., "genomes/hg38/full")
	isDef := slashCount == 0
	isRef := slashCount > 1

	// Create base object
	base := &BaseBuildObject{
		nameVersion:       normalized,
		ncpus:             defaultBuildNcpus(),
		cntDirPath:        getCntDirPath(normalized, tmpDir),
		tmpOverlayPath:    getTmpOverlayPath(normalized, tmpDir),
		targetOverlayPath: filepath.Join(imagesDir, strings.ReplaceAll(normalized, "/", "--")+".sqf"),
	}

	if external {
		// External builds don't need to resolve build source
		return createConcreteType(base, isDef, isRef, tmpDir)
	}

	// Check if already installed
	if base.IsInstalled() {
		return createConcreteType(base, isDef, isRef, tmpDir)
	}

	// Resolve build source and determine concrete type
	return createConcreteType(base, isDef, isRef, tmpDir)
}

// NewCondaObjectWithSource creates a CondaBuildObject with custom buildSource
// This is used for the -n flag to create a single sqf with multiple packages or YAML
// buildSource can be:
//   - Path to YAML file (e.g., "/path/to/environment.yml")
//   - Comma-separated package list (e.g., "nvim,nodejs,samtools/1.16")
func NewCondaObjectWithSource(nameVersion, buildSource string, imagesDir, tmpDir string) (BuildObject, error) {
	normalized := utils.NormalizeNameVersion(nameVersion)

	base := &BaseBuildObject{
		nameVersion:       normalized,
		buildSource:       buildSource,
		ncpus:             defaultBuildNcpus(),
		cntDirPath:        getCntDirPath(normalized, tmpDir),
		tmpOverlayPath:    getTmpOverlayPath(normalized, tmpDir),
		targetOverlayPath: filepath.Join(imagesDir, strings.ReplaceAll(normalized, "/", "--")+".sqf"),
	}

	return newCondaBuildObject(base)
}

// FromExternalSource creates a BuildObject from an external build script or def file
// All overlays are stored in imagesDir regardless of type
func FromExternalSource(targetPrefix, source string, isApptainer bool, imagesDir string) (BuildObject, error) {
	nameVersion := filepath.Base(targetPrefix)
	nameVersion = utils.NormalizeNameVersion(nameVersion)

	// Determine build type from source file extension
	isDef := isApptainer || strings.HasSuffix(source, ".def")
	isShell := strings.HasSuffix(source, ".sh") || strings.HasSuffix(source, ".bash")

	base := &BaseBuildObject{
		nameVersion:       nameVersion,
		buildSource:       source,
		ncpus:             defaultBuildNcpus(),
		targetOverlayPath: targetPrefix + ".sqf",
	}

	// Parse script metadata if it's a shell script
	if isShell {
		if err := base.parseScriptMetadata(); err != nil {
			return nil, err
		}
	}

	var obj BuildObject
	var err error

	if isDef {
		obj, err = newDefBuildObject(base)
	} else if isShell {
		obj, err = newScriptBuildObject(base, false) // External sources are not ref overlays
	} else {
		return nil, fmt.Errorf("unknown source type for %s", source)
	}

	if err != nil {
		return nil, err
	}

	utils.PrintDebug("[EXTERNAL] Created BuildObject from external path: %s", obj.String())
	return obj, nil
}

// createConcreteType creates the appropriate concrete BuildObject type
// It determines whether to create a conda, def, or script build based on:
// - isRef flag (from slash count > 1)
// - build source resolution: .def -> DefBuildObject, shell script -> ScriptBuildObject, not found -> conda
func createConcreteType(base *BaseBuildObject, isDef, isRef bool, tmpDir string) (BuildObject, error) {
	// Resolve build source - this determines the actual type based on file extension
	isConda, isContainer, err := resolveBuildSource(base, tmpDir)
	if err != nil {
		return nil, err
	}

	if isConda {
		return newCondaBuildObject(base)
	}

	if isContainer {
		// It's a .def file
		return newDefBuildObject(base)
	}

	// It's a shell script (no extension or .sh/.bash)
	// Parse script metadata
	if err := base.parseScriptMetadata(); err != nil {
		return nil, err
	}

	return newScriptBuildObject(base, isRef)
}

// resolveBuildSource finds the build script and determines the build type
// Returns (isConda, isContainer, error):
// - isConda=true: no build script found, use conda
// - isContainer=true: found .def file
// - both false: found shell script (no extension)
// First checks local and remote build scripts, falls back to conda if not found
func resolveBuildSource(base *BaseBuildObject, tmpDir string) (isConda bool, isContainer bool, err error) {
	// Check if already has a build source
	if base.buildSource != "" {
		// Determine type from extension
		isContainer = strings.HasSuffix(base.buildSource, ".def")
		return false, isContainer, nil
	}

	// Look for build script (local first, then remote)
	info, found := FindBuildScript(base.nameVersion)
	if !found {
		// No build script found, use conda
		utils.PrintDebug("No build script found for %s, using conda", utils.StyleName(base.nameVersion))
		return true, false, nil
	}

	// If remote, download to tmp directory
	if info.IsRemote {
		localPath, err := DownloadRemoteScript(info, tmpDir)
		if err != nil {
			return false, false, fmt.Errorf("failed to download remote build script: %w", err)
		}
		base.buildSource = localPath
		base.isRemote = true
		utils.PrintDebug("Downloaded remote build script to %s", utils.StylePath(localPath))
	} else {
		base.buildSource = info.Path
		utils.PrintDebug("Using local build script at %s", utils.StylePath(info.Path))
	}

	return false, info.IsContainer, nil
}
