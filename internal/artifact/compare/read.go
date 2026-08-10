package compare

import (
	"encoding/json"
	"fmt"
	"os"
	"path/filepath"
	"reflect"

	"github.com/Justype/condatainer/catalog"
	"github.com/Justype/condatainer/internal/artifact/meta"
	"github.com/Justype/condatainer/internal/artifact/record"
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
	// never taken on trust from the manifest.
	Identity string
	Equiv    string
	// Dependencies is the direct adjacency list, for naming which one moved.
	Dependencies []meta.Dependency

	// runtimeEnv is what the image contributes at mount time.
	runtimeEnv []meta.EnvVar
	// identity is the parsed identity record, nil for a Conda app, whose keys
	// name exports rather than records.
	identity *record.Record
	// unusable is set when the image cannot be held to a key at all.
	unusable string
}

// usable reports why this artifact cannot be compared, or "".
func (a Artifact) usable() string { return a.unusable }

// runtimeAgrees reports why runtime.json disagrees with the record beside it,
// or "".
//
// This is free integrity: runtime.json was never part of the hashed preimage, so
// an edited one cannot agree with a record it did not produce. A Conda app has no
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
		return fmt.Sprintf("%s: runtime.json contributes %v, its identity record says %v",
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
// does *not* prove is that the payload matches the records — two images with
// identical records and different payloads take no cleverness to produce, and
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

	// Imported, or built before this format. Unverifiable rather than different:
	// the question was not answered, not answered in the negative.
	if m.Keys.Identity.File == "" || m.Keys.Equiv.File == "" {
		artifact.unusable = fmt.Sprintf("%s records no keys", rt.Name)
		return artifact, nil
	}

	identity, err := verify(dir, m.Keys.Identity)
	if err != nil {
		artifact.unusable = fmt.Sprintf("%s: %v", rt.Name, err)
		return artifact, nil
	}
	equivBytes, err := verify(dir, m.Keys.Equiv)
	if err != nil {
		artifact.unusable = fmt.Sprintf("%s: %v", rt.Name, err)
		return artifact, nil
	}
	artifact.Identity = record.Digest(identity)
	artifact.Equiv = record.Digest(equivBytes)

	// A record parses; an export does not, and does not need to.
	if filepath.Ext(m.Keys.Identity.File) == ".record" {
		parsed, err := record.Parse(identity)
		if err != nil {
			artifact.unusable = fmt.Sprintf("%s has an unreadable identity record: %v", rt.Name, err)
			return artifact, nil
		}
		artifact.identity = &parsed
	}
	return artifact, nil
}

// verify reads the file a key names and checks that hashing it reproduces the
// recorded value.
func verify(dir string, ref meta.KeyRef) ([]byte, error) {
	data, err := os.ReadFile(filepath.Join(dir, ref.File))
	if err != nil {
		return nil, fmt.Errorf("manifest names %s, which the image does not carry", ref.File)
	}
	if got := record.Sum(data); got != ref.SHA256 {
		return nil, fmt.Errorf("%s hashes to %s, but the manifest says %s", ref.File, short(got), short(ref.SHA256))
	}
	return data, nil
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
