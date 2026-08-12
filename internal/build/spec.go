package build

import (
	"path/filepath"
	"strings"

	"github.com/Justype/condatainer/catalog"
	"github.com/Justype/condatainer/internal/artifact/meta"
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
}

// ImageSpec is what the image is: the identity both documents record, plus what
// loading it does to the environment.
type ImageSpec struct {
	Name        string
	Type        catalog.Type
	Description string
	URL         string
	// Prefix is where the payload goes, /cnt/<name> for app and data and empty
	// for base and os. Env holds the recipe's #ENV: contributions with {prefix}
	// intact, substituted when the image is loaded.
	Prefix string
	Env    []meta.EnvVar
	// Arch is the recipe's #ARCH: assertion. Empty means the default: the
	// artifact runs only where it was built.
	Arch catalog.Arch
}

// SourceSpec is what the build runs. Exactly one of the three kinds is set, and
// which one it is *is* the build type — see BuildType. The remaining fields
// describe the source whichever kind it turned out to be.
type SourceSpec struct {
	Script     *ScriptSource
	Definition *DefinitionSource
	Conda      *CondaSource

	// Placeholders are the selected #PH: values. The embedded recipe keeps its
	// {placeholder} tokens, so this is the only record of which variant of a
	// template was built.
	Placeholders map[string]string
	// TargetTemplate is the #TARGET: the name was rendered from; empty when the
	// recipe is not a template.
	TargetTemplate string
	// RequiresInput reports that the recipe declared #INPUT: prompts. The
	// answers are execution input and are never recorded — this says only that a
	// rebuild needs a human.
	RequiresInput bool
	// Repository is the collection's repository, recorded as manifest
	// build.source. Empty for a Conda build, a local file, or a collection that
	// declares none.
	Repository string
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
	// From is the upstream image the definition bootstraps from, resolved before
	// the build runs. Nil when there is no upstream.
	From *meta.From
}

// CondaSource is a micromamba environment. Exactly one of the three inputs is set.
type CondaSource struct {
	Package  *CondaPackage // a single name/version
	Packages []string      // an explicit list
	File     *SourceFile   // a YAML or explicit-spec input

	// Channels are the configured channels the solve was offered, in priority
	// order. Captured here rather than read from config when the manifest is
	// rendered, so the manifest stays a projection of Spec and nothing else.
	Channels []string
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

// RecipeFile returns the recipe this build embeds at /.cnt/recipe, and whether
// there is one. A Conda build has none: its inputs are the exports captured from
// the environment it installed, not the request that produced it.
func (s SourceSpec) RecipeFile() (SourceFile, bool) {
	switch {
	case s.Script != nil:
		return s.Script.File, true
	case s.Definition != nil:
		return s.Definition.File, true
	}
	return SourceFile{}, false
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
	MetaDir  string // BuildDir/.cnt — runtime and manifest, staged before packing
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
	prefix := spec.Image.Prefix
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

// Manifest renders what the image records about itself and where it came from. A
// projection of Spec and nothing else, so no host, job ID, local path or #INPUT:
// answer can reach it — none of them is in Spec to begin with.
func (s Spec) Manifest() meta.Manifest {
	return meta.Manifest{
		SchemaVersion: meta.SchemaVersion,
		Name:          s.Image.Name,
		Type:          s.Image.Type,
		BuildType:     s.Source.BuildType().String(),
		Description:   s.Image.Description,
		URL:           s.Image.URL,
		Platform:      s.platform(),
		Source:        s.sourceBlock(),
		Build:         s.buildBlock(),
	}
}

// platform is where this artifact was built and where it may run. Only a recipe
// that declared #ARCH:noarch is portable — nothing infers portability, because
// getting it wrong does not crash, it silently returns wrong answers.
func (s Spec) platform() meta.Platform {
	if s.Image.Arch == catalog.ArchNoarch {
		return meta.Platform{OS: "linux", Arch: meta.ArchNone}
	}
	return meta.NativePlatform()
}

// sourceBlock describes what the image was built from. Files is left to the
// caller: what is embedded is decided when it is staged, and a manifest that
// named a file the image does not carry would be worse than one that named none.
func (s Spec) sourceBlock() meta.Source {
	return meta.Source{
		Placeholders:   s.Source.Placeholders,
		TargetTemplate: s.Source.TargetTemplate,
		RequiresInput:  s.Source.RequiresInput,
	}
}

// buildBlock is what the build knew and the recipe does not say. Created is not
// set here: a timestamp is not in Spec, so stageMetadata stamps it.
func (s Spec) buildBlock() meta.Build {
	block := meta.Build{Source: s.Source.Repository}
	switch {
	case s.Source.Conda != nil:
		block.Channels = s.Source.Conda.Channels
	case s.Source.Definition != nil:
		block.From = s.Source.Definition.From
	}
	return block
}

// UpstreamDigest returns what the definition's bootstrap reference resolved to.
// Empty means no upstream at all, which the record distinguishes from an
// unresolved one by omitting the line rather than writing meta.Unrecorded.
func (s SourceSpec) UpstreamDigest() string {
	if s.Definition == nil || s.Definition.From == nil {
		return ""
	}
	return s.Definition.From.Digest
}

// Runtime renders the mount-time contract: what the image is called, where its
// payload sits, and what it contributes to the environment. Container setup
// reads this and nothing else, so it repeats the few identity fields it needs
// rather than sending a reader to the manifest for them.
func (s Spec) Runtime() meta.Runtime {
	return meta.Runtime{
		SchemaVersion: meta.SchemaVersion,
		Name:          s.Image.Name,
		Type:          s.Image.Type,
		Description:   s.Image.Description,
		Platform:      s.platform(),
		Prefix:        s.Image.Prefix,
		Env:           s.Image.Env,
	}
}

// envFromRecipe converts a recipe's #ENV: declarations. Values keep {prefix}
// intact, for the loader to substitute at mount time.
func envFromRecipe(env []catalog.EnvVar) []meta.EnvVar {
	var out []meta.EnvVar
	for _, e := range env {
		out = append(out, meta.EnvVar{
			Key:   e.Key,
			Value: e.Value(nil),
			Note:  e.Note,
		})
	}
	return out
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
