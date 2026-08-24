package compare

import (
	"encoding/json"
	"fmt"
	"os"
	"path/filepath"
	"reflect"
	"strings"

	"github.com/Justype/condatainer/catalog"
	"github.com/Justype/condatainer/internal/artifact/key"
	"github.com/Justype/condatainer/internal/artifact/meta"
	"github.com/Justype/condatainer/internal/image"
)

// Artifact is one side of a comparison: what an image says it is, and the keys
// it can be held to. Read builds one from an image; a build can construct one
// from what it is about to produce.
type Artifact struct {
	Name string
	Type catalog.Type
	// Arch is the recorded architecture, or "noarch".
	Arch string
	// Format is the build type. Two artifacts of different formats are never
	// the same build, whatever their digests: the digests are over different
	// kinds of file.
	Format string
	// Identity and Equiv are recomputed from the files manifest.keys names,
	// never taken on trust from the manifest. Both carry the "sha256:" prefix,
	// which meta.KeyRef does not — cross that boundary with IdentityRef and
	// EquivRef rather than by hand.
	Identity       string
	Equiv          string
	IdentityScheme string
	EquivScheme    string
	// Dependencies is the direct adjacency list, for naming which one moved.
	Dependencies []meta.Dependency

	// runtimeEnv is what the image contributes at mount time.
	runtimeEnv []meta.EnvVar
	// identity is the canonical identity model, nil for a Conda app, whose keys
	// derive directly from exports.
	identity *key.Model
	// unusable is set when the image cannot be held to a key at all.
	unusable string
}

// IdentityRef is the identity as a meta.KeyRef.
//
// An Artifact digest is "sha256:<hex>" and a meta.KeyRef holds the bare hex, so
// assigning one to the other directly yields a key that compares equal to
// nothing and fails validation. That mistake has been made at four separate
// call sites; this is the only correct crossing.
func (a Artifact) IdentityRef() meta.KeyRef {
	return meta.KeyRef{Scheme: a.IdentityScheme, SHA256: strings.TrimPrefix(a.Identity, "sha256:")}
}

// EquivRef is the equivalence key as a meta.KeyRef. See IdentityRef.
func (a Artifact) EquivRef() meta.KeyRef {
	return meta.KeyRef{Scheme: a.EquivScheme, SHA256: strings.TrimPrefix(a.Equiv, "sha256:")}
}

// usable reports why this artifact cannot be compared, or "".
func (a Artifact) usable() string { return a.unusable }

// runtimeAgrees reports why runtime.json disagrees with the regenerated identity model,
// or "".
//
// This is free integrity: runtime.json was never part of the hashed preimage, so
// an edited one cannot agree with the regenerated model. A Conda app has no
// env= lines and must therefore contribute no variables.
func (a Artifact) runtimeAgrees() string {
	var recorded []meta.EnvVar
	if a.identity != nil {
		for _, e := range a.identity.Env {
			recorded = append(recorded, meta.EnvVar{Key: e.Key, Value: e.Value})
		}
	}
	live := make([]meta.EnvVar, 0, len(a.runtimeEnv))
	for _, e := range a.runtimeEnv {
		live = append(live, meta.EnvVar{Key: e.Key, Value: e.Value}) // notes are not recorded
	}
	if len(recorded) == 0 && len(live) == 0 {
		return ""
	}
	if !reflect.DeepEqual(recorded, live) {
		return fmt.Sprintf("%s: runtime.json contributes %v, its identity model says %v",
			a.Name, keysOf(live), keysOf(recorded))
	}
	return ""
}

func keysOf(env []meta.EnvVar) []string {
	out := make([]string, 0, len(env))
	for _, e := range env {
		out = append(out, e.Key)
	}
	return out
}

// Read builds an Artifact from an image, in one archive extraction.
//
// Every file-backed hash is recomputed from the file manifest.keys names rather
// than trusted: a manifest is what the publisher wrote, and a key that does not
// match its own file is exactly what verification exists to catch. What this
// does *not* prove is that the payload matches the derived keys — two images with
// identical keys and different payloads take no cleverness to produce, and
// that needs a signature policy rather than a hash.
func Read(imagePath string) (Artifact, error) {
	dir, err := os.MkdirTemp("", "cnt-compare-")
	if err != nil {
		return Artifact{}, fmt.Errorf("failed to create scratch dir: %w", err)
	}
	defer os.RemoveAll(dir) //nolint:errcheck

	// One extraction, not one read per file: every archive read spawns a process.
	if err := image.ExtractDir(imagePath, "/"+meta.DirName, dir); err != nil {
		return Artifact{}, fmt.Errorf("%s carries no CondaTainer metadata: %w", imagePath, err)
	}
	return fromDir(filepath.Join(dir, meta.DirName), imagePath)
}

// fromDir assembles an Artifact from an extracted /.cnt.
func fromDir(dir, imagePath string) (Artifact, error) {
	var rt meta.Runtime
	if err := readJSON(filepath.Join(dir, meta.RuntimeFileName), &rt); err != nil {
		return Artifact{}, fmt.Errorf("%s: %w", imagePath, err)
	}
	rt.Normalize()
	if err := meta.ValidateRuntime(rt); err != nil {
		return Artifact{}, fmt.Errorf("%s: %w", imagePath, err)
	}

	artifact := Artifact{
		Name:       rt.Name,
		Type:       rt.Type,
		Arch:       rt.Platform.Arch,
		runtimeEnv: rt.Env,
	}

	var m meta.Manifest
	if err := readJSON(filepath.Join(dir, meta.FileName), &m); err != nil {
		artifact.unusable = fmt.Sprintf("%s has no readable manifest", rt.Name)
		return artifact, nil
	}
	m.Normalize()
	if err := meta.ValidateManifest(m); err != nil {
		artifact.unusable = fmt.Sprintf("%s has an invalid manifest: %v", rt.Name, err)
		return artifact, nil
	}
	artifact.Format = m.BuildType
	artifact.Dependencies = m.Dependencies

	// Imported, or built without a complete source capture. Unverifiable rather
	// than different: the question was not answered in the negative.
	if m.Keys.Identity.Empty() {
		artifact.unusable = fmt.Sprintf("%s records no keys", rt.Name)
		return artifact, nil
	}

	derived, err := key.VerifyDir(dir, m)
	if err != nil {
		artifact.unusable = fmt.Sprintf("%s: %v", rt.Name, err)
		return artifact, nil
	}
	artifact.Identity = derived.Identity.Ref.Digest()
	artifact.Equiv = derived.Equiv.Ref.Digest()
	artifact.IdentityScheme = derived.Identity.Ref.Scheme
	artifact.EquivScheme = derived.Equiv.Ref.Scheme
	artifact.identity = derived.IdentityModel
	return artifact, nil
}

func readJSON(path string, into any) error {
	data, err := os.ReadFile(path)
	if err != nil {
		return fmt.Errorf("cannot read %s: %w", filepath.Base(path), err)
	}
	if err := json.Unmarshal(data, into); err != nil {
		return fmt.Errorf("cannot decode %s: %w", filepath.Base(path), err)
	}
	return nil
}
