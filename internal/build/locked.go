package build

import (
	"bytes"
	"context"
	"errors"
	"fmt"
	"path/filepath"
	"slices"

	"github.com/Justype/condatainer/catalog"
	"github.com/Justype/condatainer/internal/artifact/meta"
	"github.com/Justype/condatainer/internal/conda"
	"github.com/Justype/condatainer/internal/config"
	"github.com/Justype/condatainer/internal/utils"
)

// ErrLockedInvalid reports that a rebuild specification cannot be started: its
// vendored records disagree with each other, or the caller supplied a
// destination or dependency set the records do not describe.
var ErrLockedInvalid = errors.New("locked rebuild specification is invalid")

// LockedDep is one build dependency that has already been materialized, given
// as a path rather than a name.
type LockedDep struct {
	// Name and Identity are what the dependent's manifest edge records. They are
	// reported in errors and never used to look anything up.
	Name     string
	Identity meta.KeyRef
	// Path is the image to mount, absolute.
	Path string
}

// LockedSpec is a rebuild described entirely by a project lock: vendored
// records, exact dependency paths, and a destination the caller owns.
type LockedSpec struct {
	// Manifest is the vendored manifest, whose keys the result must reproduce.
	Manifest meta.Manifest
	// Sources are the vendored build inputs by manifest.source.files name.
	Sources map[string][]byte
	// Deps are the dependency images to mount, in the manifest's dependency
	// order. The rebuild re-derives its own edges from what it mounts, so an
	// order that does not match the manifest's produces a key mismatch rather
	// than a wrong artifact.
	Deps []LockedDep
	// Output is where the finished image lands. The caller owns it: nothing
	// here installs, and a flat name is never claimed.
	Output string
	// Answers are the recipe's #INPUT: answers, in declaration order. Supplied
	// per invocation because they are never recorded.
	Answers []string
	// TmpRoot overrides where the build works. Empty uses the configured root
	// for the artifact's type.
	TmpRoot string
	// CondaSource names which vendored Conda export to replay: explicit.txt by
	// default, or environment.yml when the caller has accepted an equivalent
	// result. Ignored by every other build type.
	CondaSource string
}

// NewLockedObject builds one artifact from a lock rather than from the catalog.
//
// It reuses the ordinary build implementations, workspace, base resolution and
// scheduler parsing, and differs in what it refuses to do: no catalog lookup, no
// fuzzy version choice, no prebuilt pull by name, no flat install. Dependencies
// are mounted exactly as supplied and the output goes exactly where the caller
// said.
//
// Nothing here checks the result. The rebuild derives its own keys from what it
// actually mounted and built, and the caller compares them against the lock —
// which is why an inconsistency shows up as a key mismatch rather than as a
// plausible artifact under the right name.
func NewLockedObject(ctx context.Context, spec LockedSpec) (*BuildObject, error) {
	manifest := spec.Manifest
	manifest.Normalize()
	if err := meta.ValidateManifest(manifest); err != nil {
		return nil, fmt.Errorf("%w: %v", ErrLockedInvalid, err)
	}
	if manifest.Keys.Identity.Empty() || manifest.Keys.Equiv.Empty() {
		return nil, fmt.Errorf("%w: %s records no complete keys, so nothing could verify a rebuild",
			ErrLockedInvalid, manifest.Name)
	}
	output, err := lockedOutput(spec.Output, manifest.Name)
	if err != nil {
		return nil, err
	}
	for _, dep := range spec.Deps {
		if !filepath.IsAbs(dep.Path) {
			return nil, fmt.Errorf("%w: dependency %s must be an absolute path, got %q",
				ErrLockedInvalid, dep.Name, dep.Path)
		}
		if !utils.FileExists(dep.Path) {
			return nil, fmt.Errorf("%w: dependency %s is not at %s", ErrLockedInvalid, dep.Name, dep.Path)
		}
	}

	tmpRoot := spec.TmpRoot
	if tmpRoot == "" {
		tmpRoot = tmpRootForType(manifest.Type)
	}
	if abs, err := filepath.Abs(tmpRoot); err == nil {
		tmpRoot = abs
	}

	b := &BuildObject{
		spec: Spec{Image: ImageSpec{
			Name:   manifest.Name,
			Type:   manifest.Type,
			Prefix: meta.Prefix(manifest.Name, manifest.Type),
		}},
		ws:           workspaceFor(manifest.Name, tmpRoot, appExt3ScratchExt(manifest.Type)),
		tgt:          targetFor(output),
		submitJob:    config.Global.SubmitJob,
		locked:       true,
		lockedKeys:   manifest.Keys,
		inputAnswers: slices.Clone(spec.Answers),
	}
	// A path is what makes the mount exact. Everything downstream already reads
	// a dependency's own manifest for its name and keys, so recording the edge
	// is unaffected by having located it this way.
	for _, dep := range spec.Deps {
		b.spec.Dependencies = append(b.spec.Dependencies, dep.Path)
	}
	b.spec.Source.Repository = manifest.Build.Source

	switch manifest.BuildType {
	case BuildTypeScript.String(), BuildTypeDef.String():
		err = b.lockRecipeSource(manifest, spec.Sources)
	case BuildTypeConda.String():
		err = b.lockCondaSource(manifest, spec.Sources, spec.CondaSource)
	default:
		err = fmt.Errorf("%w: %s records build type %q", ErrLockedInvalid, manifest.Name, manifest.BuildType)
	}
	if err != nil {
		return nil, err
	}
	if err := b.resolveResourceSpec(); err != nil {
		return nil, err
	}
	return b, nil
}

// lockedOutput validates a caller-owned destination: an absolute path in an
// existing directory, and nothing that is already there. Restore builds into a
// temporary sibling and renames, so an occupied output is a caller bug rather
// than a race worth tolerating.
func lockedOutput(output, name string) (string, error) {
	if output == "" {
		return "", fmt.Errorf("%w: rebuilding %s needs an output path", ErrLockedInvalid, name)
	}
	absolute, err := filepath.Abs(output)
	if err != nil {
		return "", err
	}
	if !utils.DirExists(filepath.Dir(absolute)) {
		return "", fmt.Errorf("%w: output directory %s does not exist", ErrLockedInvalid, filepath.Dir(absolute))
	}
	if utils.FileExists(absolute) {
		return "", fmt.Errorf("%w: output %s already exists", ErrLockedInvalid, absolute)
	}
	return absolute, nil
}

// lockRecipeSource fills the Spec from a vendored recipe, exactly as a catalog
// build fills it from a materialized one.
//
// A template is expanded with the placeholders the manifest recorded and the
// expansion is what runs, while what the image embeds is the template itself —
// the same split the catalog path makes, and the reason the recipe's bytes
// still regenerate the recorded identity.
func (b *BuildObject) lockRecipeSource(manifest meta.Manifest, sources map[string][]byte) error {
	text, ok := sources[meta.RecipeFileName]
	if !ok {
		return fmt.Errorf("%w: %s vendors no %s", ErrLockedInvalid, manifest.Name, meta.RecipeFileName)
	}
	recipe, err := catalog.ParseRecipe(meta.RecipeFileName, bytes.NewReader(text))
	if err != nil {
		return fmt.Errorf("%w: cannot parse the vendored recipe for %s: %v", ErrLockedInvalid, manifest.Name, err)
	}
	runnable := recipe
	if recipe.IsTemplate {
		if runnable, err = catalog.Expand(recipe, manifest.Source.Placeholders); err != nil {
			return fmt.Errorf("%w: cannot expand the vendored recipe for %s: %v", ErrLockedInvalid, manifest.Name, err)
		}
	}
	if len(runnable.Deps) != len(b.spec.Dependencies) {
		return fmt.Errorf("%w: %s declares %d dependencies but %d paths were supplied",
			ErrLockedInvalid, manifest.Name, len(runnable.Deps), len(b.spec.Dependencies))
	}

	isDef := manifest.BuildType == BuildTypeDef.String()
	file := SourceFile{Name: meta.RecipeFileName, Data: recipe.Text}
	source := SourceSpec{Script: &ScriptSource{File: file, Prompts: runnable.Inputs}}
	if isDef {
		source = SourceSpec{Definition: &DefinitionSource{File: file, From: manifest.Build.From}}
	}
	// Placeholders and the target template come from the manifest, not from the
	// expansion: they are what was recorded, and a rebuild has to reproduce the
	// record rather than re-derive it.
	source.Placeholders = manifest.Source.Placeholders
	source.TargetTemplate = manifest.Source.TargetTemplate
	source.RequiresInput = manifest.Source.RequiresInput
	source.Repository = manifest.Build.Source
	b.spec.Source = source
	// The recipe is authoritative for everything it declares, as it is on the
	// catalog path. Anything taken from the manifest instead would be a second
	// copy that could disagree with the bytes the identity was derived from.
	b.spec.Image.Description = runnable.Description
	b.spec.Image.URL = runnable.URL
	b.spec.Image.Env = envFromRecipe(runnable.Env)
	b.spec.Image.Arch = runnable.Arch
	b.embedSource(file)

	b.inputPrompts = runnable.Inputs
	if b.spec.Source.RequiresInput && len(b.inputAnswers) != len(b.inputPrompts) {
		return fmt.Errorf("%w: %s declares %d #INPUT: prompts but %d answers were supplied",
			ErrLockedInvalid, manifest.Name, len(b.inputPrompts), len(b.inputAnswers))
	}

	path, err := writeRecipeFile(runnable, b.ws.Source)
	if err != nil {
		return err
	}
	b.buildSource = path
	b.tempSource = true
	if isDef {
		// Re-sites the workspace as a .sif build. The Spec is already set, so
		// this does not go looking for a local source to capture one from.
		b.asDefinitionBuild()
		// Re-siting moved ws.Source, so the recipe has to follow it.
		if path, err = writeRecipeFile(runnable, b.ws.Source); err != nil {
			return err
		}
		b.buildSource = path
		return nil
	}
	b.buildType = BuildTypeScript
	return nil
}

// lockCondaSource fills the Spec from the vendored Conda export named by which,
// defaulting to explicit.txt.
//
// explicit.txt is the only form that reproduces the recorded identity, because
// that identity *is* its bytes. environment.yml re-solves and reproduces the
// equivalence instead — it is the equivalence preimage — which is what a caller
// asks for by accepting an equivalent result.
func (b *BuildObject) lockCondaSource(manifest meta.Manifest, sources map[string][]byte, which string) error {
	if which == "" {
		which = conda.ExplicitFileName
	}
	if which != conda.ExplicitFileName && which != conda.EnvironmentFileName {
		return fmt.Errorf("%w: %q is not a Conda export", ErrLockedInvalid, which)
	}
	export, ok := sources[which]
	if !ok {
		return fmt.Errorf("%w: %s vendors no %s, so its package set cannot be reproduced",
			ErrLockedInvalid, manifest.Name, which)
	}
	if len(b.spec.Dependencies) > 0 {
		return fmt.Errorf("%w: %s is a Conda build and cannot have dependencies", ErrLockedInvalid, manifest.Name)
	}
	path := filepath.Join(b.ws.Root, which)
	if err := utils.MkdirAllShared(filepath.Dir(path)); err != nil {
		return err
	}
	file, err := utils.CreateFileWritable(path)
	if err != nil {
		return err
	}
	if _, err := file.Write(export); err != nil {
		file.Close() //nolint:errcheck
		return err
	}
	if err := file.Close(); err != nil {
		return err
	}
	utils.ShareWithParentGroup(path)

	b.buildSource = path
	b.tempSource = true
	if err := b.setupCondaFields(); err != nil {
		return fmt.Errorf("%w: %v", ErrLockedInvalid, err)
	}
	b.buildType = BuildTypeConda
	b.setCondaSpec()
	// Channels are diagnostic, but the manifest is a projection of the Spec, so
	// taking them from live config would rewrite what the lock recorded.
	b.spec.Source.Conda.Channels = slices.Clone(manifest.Build.Channels)
	b.spec.Source.Repository = manifest.Build.Source
	// A Conda build has no recipe to take a description from — the original
	// fetched one from anaconda.org, and the manifest is where it survives.
	b.spec.Image.Description = manifest.Description
	b.spec.Image.URL = manifest.URL
	return nil
}
