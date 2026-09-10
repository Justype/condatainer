package container

import (
	"path/filepath"
	"strings"
	"testing"

	"github.com/Justype/condatainer/catalog"
	"github.com/Justype/condatainer/internal/artifact/meta"
)

// A bare env-typed .sqf (no paired .img) must still put /cnt_env/bin on
// PATH — it is an honest TypeEnv contribution, not merely a mislabeled app.
func TestBuildPathEnvIncludesBareEnvSqf(t *testing.T) {
	requireSquashfsTools(t)
	dir := t.TempDir()
	envSqf := filepath.Join(dir, "env.sqf")
	packRuntimeSqf(t, envSqf, envRuntime())

	got := BuildPathEnv([]string{envSqf})
	if !containsPathEntry(got, EnvPrefix+"/bin") {
		t.Errorf("PATH missing %s/bin for a bare env.sqf mount: %s", EnvPrefix, got)
	}
}

// A writable .img and its paired env.sqf snapshot both claim EnvPrefix but
// merge into the one physical /cnt_env directory at mount time, so PATH
// needs only one entry for it, not one per overlay that claims the prefix.
func TestBuildPathEnvDedupesPairedImgAndSnapshot(t *testing.T) {
	requireSquashfsTools(t)
	dir := t.TempDir()
	envSqf := filepath.Join(dir, "env.sqf")
	packRuntimeSqf(t, envSqf, envRuntime())
	img := filepath.Join(dir, "env.img")

	got := BuildPathEnv([]string{envSqf, img})
	if n := strings.Count(got, EnvPrefix+"/bin"); n != 1 {
		t.Errorf("want exactly one %s/bin entry, got %d: %s", EnvPrefix, n, got)
	}
}

// An app overlay still contributes its own <prefix>/bin, unaffected by
// ContributesBin() also accepting TypeEnv.
func TestBuildPathEnvStillIncludesApp(t *testing.T) {
	requireSquashfsTools(t)
	dir := t.TempDir()
	appSqf := filepath.Join(dir, "myapp.sqf")
	packRuntimeSqf(t, appSqf, meta.Runtime{
		SchemaVersion: meta.SchemaVersion,
		Name:          "myapp/1.0",
		Type:          catalog.TypeApp,
		Platform:      meta.NativePlatform(),
		Prefix:        "/cnt/myapp",
	})

	got := BuildPathEnv([]string{appSqf})
	if !containsPathEntry(got, "/cnt/myapp/bin") {
		t.Errorf("PATH missing /cnt/myapp/bin: %s", got)
	}
}

func containsPathEntry(pathEnv, entry string) bool {
	for _, p := range strings.Split(pathEnv, ":") {
		if p == entry {
			return true
		}
	}
	return false
}
