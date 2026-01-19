package build

import (
	"fmt"
	"os"
	"path/filepath"
	"strings"

	"github.com/Justype/condatainer/internal/utils"
)

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

	// SLURM integration
	NeedsSbatch() bool
	SbatchFlags() []string
	Ncpus() int

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
	sbatchFlags       []string
	ncpus             int
	isRemote          bool // Whether build source was downloaded

	// SLURM
	needsSbatch bool

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
func (b *BaseBuildObject) NeedsSbatch() bool         { return b.needsSbatch }
func (b *BaseBuildObject) SbatchFlags() []string     { return b.sbatchFlags }
func (b *BaseBuildObject) Ncpus() int                { return b.ncpus }

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
			return fmt.Errorf("temporary overlay already exists at %s", b.tmpOverlayPath)
		}
		// Remove existing if force=true
		if err := os.Remove(b.tmpOverlayPath); err != nil {
			return fmt.Errorf("failed to remove existing tmp overlay: %w", err)
		}
	}

	// TODO: Create actual ext3 overlay using overlay package
	// For now, just a placeholder
	utils.PrintDebug("Creating temporary overlay at %s", utils.StylePath(b.tmpOverlayPath))
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

	// TODO: Implement parsing of:
	// - #DEP directives for dependencies
	// - #SBATCH directives for scheduler flags
	// - #PROMPT directives for interactive input
	// For now, just set empty values
	b.dependencies = []string{}
	b.sbatchFlags = []string{}
	b.interactiveInputs = []string{}
	b.needsSbatch = len(b.sbatchFlags) > 0

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
	filename := strings.ReplaceAll(nameVersion, "/", "--") + ".sqf"
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
		ncpus:             1, // Default, can be overridden from sbatch flags
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
// - isDef flag (from slash count == 0)
// - isRef flag (from slash count > 1)
// - build source resolution (conda if no build script found)
func createConcreteType(base *BaseBuildObject, isDef, isRef bool, tmpDir string) (BuildObject, error) {
	if isDef {
		// No slashes = def file
		if err := resolveBuildSourceForDef(base, tmpDir); err != nil {
			return nil, err
		}
		return newDefBuildObject(base)
	}

	// Has slashes = conda or shell script
	// Resolve build source to determine if it's conda or shell
	isConda, err := resolveBuildSourceForPackage(base, tmpDir)
	if err != nil {
		return nil, err
	}

	if isConda {
		return newCondaBuildObject(base)
	}

	// It's a shell script
	// Parse script metadata
	if err := base.parseScriptMetadata(); err != nil {
		return nil, err
	}

	return newScriptBuildObject(base, isRef)
}

// resolveBuildSourceForDef finds the build script for def files
func resolveBuildSourceForDef(base *BaseBuildObject, tmpDir string) error {
	// TODO: Implement local and remote script discovery
	if base.buildSource == "" {
		return fmt.Errorf("build script for %s not found locally or remotely", base.nameVersion)
	}
	return nil
}

// resolveBuildSourceForPackage determines if a package is conda or has a build script
// Returns true if it's a conda package, false if it's a shell script
func resolveBuildSourceForPackage(base *BaseBuildObject, tmpDir string) (bool, error) {
	// TODO: Implement local and remote script discovery
	// For now, if no build source found, assume conda package
	if base.buildSource == "" {
		return true, nil // It's a conda package
	}
	return false, nil // It has a build script, so it's a shell script
}
