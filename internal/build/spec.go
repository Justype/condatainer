package build

import (
	"path/filepath"
	"strings"

	"github.com/Justype/condatainer/catalog"
	"github.com/Justype/condatainer/internal/image/meta"
	"github.com/Justype/condatainer/internal/scheduler"
)

// A build has four parts, split by lifetime — Spec, Options, Workspace, Target —
// and BuildObject holds them. See the README's Key Types.

// Spec is the complete description of the image to produce: what is being built,
// never how this invocation was launched.
type Spec struct {
	Image        ImageSpec
	Source       SourceSpec
	Base         string // normalized type=base image this build runs inside; empty when none
	Dependencies []string
	Runtime      meta.Runtime
}

// ImageSpec is the identity the manifest records.
type ImageSpec struct {
	Name        string
	Type        catalog.Type
	Description string
	URL         string
}

// SourceSpec is what the build runs. Exactly one field is set, and which one it
// is *is* the build type — see BuildType.
type SourceSpec struct {
	Script     *ScriptSource
	Definition *DefinitionSource
	Conda      *CondaSource
}

// SourceFile is a build input captured during resolution, already expanded.
// Data rather than a path, so the origin stops mattering once resolution is done.
type SourceFile struct {
	Name string
	Data []byte
}

// ScriptSource is a recipe run as a shell script.
type ScriptSource struct {
	File    SourceFile
	Prompts []string // #INPUT: declarations, in order; the answers are execution input
}

// DefinitionSource is an Apptainer definition.
type DefinitionSource struct {
	File SourceFile
}

// CondaSource is a micromamba environment. Exactly one of the three is set.
type CondaSource struct {
	Package  *CondaPackage // a single name/version
	Packages []string      // an explicit list
	File     *SourceFile   // a YAML or explicit-spec input
}

// CondaPackage is one primary package.
type CondaPackage struct {
	Name           string
	Version        string
	ChannelPackage string // channel-annotated spec, e.g. "bioconda::star"
}

// BuildType reports how this image is produced, derived from which source is set
// rather than stored beside it — storing it would need a rule that the two agree.
func (s SourceSpec) BuildType() BuildType {
	switch {
	case s.Conda != nil:
		return BuildTypeConda
	case s.Definition != nil:
		return BuildTypeDef
	case s.Script != nil:
		return BuildTypeScript
	}
	return 0
}

// File returns the source file this build materializes into the workspace, and
// whether there is one. A single Conda package or package list has none.
func (s SourceSpec) File() (SourceFile, bool) {
	switch {
	case s.Script != nil:
		return s.Script.File, true
	case s.Definition != nil:
		return s.Definition.File, true
	case s.Conda != nil && s.Conda.File != nil:
		return *s.Conda.File, true
	}
	return SourceFile{}, false
}

// Options is what this invocation chose, as opposed to what the image is.
type Options struct {
	Update      bool
	ScriptSpecs *scheduler.ScriptSpecs // resolved resource spec; never nil once built
}

// Workspace is where a build does its work. All of it is removed afterwards.
type Workspace struct {
	Root     string // scratch root this build works under
	BuildDir string // Root/build_<name>, holding the three below
	CntDir   string // BuildDir/cnt — payload root, bound as /cnt
	TmpDir   string // BuildDir/tmp — scratch, bound as /cnt_tmp
	MetaDir  string // BuildDir/.cnt — manifest, staged before packing
	Source   string // materialized script, definition, or Conda input file
	Overlay  string // Root/<name-->.img|.sif; "" in directory mode
}

// UsesImage reports whether the build runs inside a scratch image rather than
// host directories. It is the mode, decided once at construction: reading the
// config again later is how the two came apart.
func (w Workspace) UsesImage() bool { return w.Overlay != "" }

// HostPayload reports whether the payload is written to the host rather than
// inside the scratch image — which is exactly "there is no image", since only an
// app build gets one and an app build fills it. See appExt3ScratchExt.
func (w Workspace) HostPayload() bool { return !w.UsesImage() }

// Target is where the finished image lands.
type Target struct {
	Path     string // the installed image
	Prepared string // temporary output, renamed over Path on success
	Lock     string // Path + ".lock"
}

// ScratchPath is where a build's scratch space appears inside the container,
// whichever workspace mode is in use — it is what $CNT_TMP points at.
const ScratchPath = "/cnt_tmp"

// buildEnv is the environment a script recipe runs with: the image from Spec,
// the machine from Options. NCPUS, MEM and MEM_GB are appended unprefixed — the
// exception to CNT_, so a recipe copied from a cluster's docs still works.
func buildEnv(spec Spec, opts Options) []string {
	effRS := buildEffectiveResourceSpec(opts.ScriptSpecs)
	cpus := effRS.CpusPerTask
	if effRS.TasksPerNode > 1 {
		cpus *= effRS.TasksPerNode
	}
	buildRS := &scheduler.ResourceSpec{
		Nodes:        1,
		TasksPerNode: 1,
		CpusPerTask:  cpus,
		MemPerCpuMB:  effRS.MemPerCpuMB,
		MemPerNodeMB: effRS.MemPerNodeMB,
	}

	typ := spec.Image.Type
	if typ == "" {
		typ = catalog.TypeApp
	}
	prefix := spec.Runtime.Prefix
	if prefix == "" {
		prefix = meta.Prefix(spec.Image.Name, typ)
	}

	// CNT_NAME is the complete name and there is no CNT_VERSION: not every image
	// has one version axis. A recipe that varies by version uses a #PH: instead.
	return append(scheduler.ResourceEnvVars(buildRS),
		"CNT_NAME="+spec.Image.Name,
		"CNT_TYPE="+string(typ),
		"CNT_PREFIX="+prefix,
		"CNT_TMP="+ScratchPath,
		"TMPDIR="+ScratchPath,
		"IN_CONDATAINER=1",
	)
}

// Manifest renders the metadata this build embeds. A projection of Spec and
// nothing else, so no host, job ID, local path or #INPUT: answer can reach it —
// none of them is in Spec to begin with.
func (s Spec) Manifest() meta.Manifest {
	return meta.Manifest{
		SchemaVersion: meta.SchemaVersion,
		Name:          s.Image.Name,
		Type:          s.Image.Type,
		BuildType:     s.Source.BuildType().String(),
		Description:   s.Image.Description,
		URL:           s.Image.URL,
		Runtime:       s.Runtime,
	}
}

// runtimeFromRecipe builds the runtime block from a recipe's #ENV: declarations.
// Values keep {prefix} intact, for the manifest to substitute at load time.
func runtimeFromRecipe(name string, typ catalog.Type, env []catalog.EnvVar) meta.Runtime {
	rt := meta.Runtime{Prefix: meta.Prefix(name, typ)}
	for _, e := range env {
		rt.Env = append(rt.Env, meta.EnvVar{
			Key:   e.Key,
			Value: e.Value(nil),
			Note:  e.Note,
		})
	}
	return rt
}

// sourceFileName is the workspace filename for a materialized recipe.
func sourceFileName(name string, isDef bool) string {
	base := "cnt--" + strings.ReplaceAll(name, "/", "--")
	if isDef {
		return base + ".def"
	}
	return base + ".sh"
}

// workspaceFor derives the whole path set for a build under root; every
// constructor goes through it. An empty overlayExt means directory mode: no
// scratch image, so the payload lands on the host.
func workspaceFor(name, root, overlayExt string) Workspace {
	cntDir := getCntDirPath(name, root)
	buildDir := filepath.Dir(cntDir)
	ws := Workspace{
		Root:     root,
		BuildDir: buildDir,
		CntDir:   cntDir,
		TmpDir:   filepath.Join(buildDir, "tmp"),
		MetaDir:  filepath.Join(buildDir, meta.DirName),
	}
	if overlayExt != "" {
		ws.Overlay = filepath.Join(root, strings.ReplaceAll(name, "/", "--")+overlayExt)
	}
	return ws
}

// targetFor derives the target paths for an installed image. Prepared is left
// empty: it derives from the lock owner, so it is unknowable until the lock is
// held. Not yet adopted — see workspaceFor.
func targetFor(path string) Target {
	if abs, err := filepath.Abs(path); err == nil {
		path = abs
	}
	return Target{Path: path, Lock: path + ".lock"}
}
