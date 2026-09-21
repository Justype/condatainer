package meta

import (
	"os"
	"os/exec"
	"path/filepath"
	"strings"
	"testing"

	"github.com/Justype/condatainer/catalog"
	artifactcache "github.com/Justype/condatainer/internal/artifact/cache"
)

// validRuntime returns a runtime document that passes ValidateRuntime, for a
// test to then break in exactly one way.
func validRuntime() Runtime {
	return Runtime{
		SchemaVersion: SchemaVersion,
		Name:          "samtools/1.23.1",
		Type:          catalog.TypeApp,
		Platform:      NativePlatform(),
		Prefix:        "/cnt/samtools/1.23.1",
		Env:           []EnvVar{{Key: "SAMTOOLS_HOME", Value: "{prefix}"}},
	}
}

// validManifest returns a manifest that passes ValidateManifest.
func validManifest() Manifest {
	return Manifest{
		SchemaVersion: SchemaVersion,
		Name:          "samtools/1.23.1",
		Type:          catalog.TypeApp,
		BuildType:     "conda",
		Platform:      NativePlatform(),
	}
}

// withTempCache points the runtime cache at a scratch file so a test never reads
// or writes the user's real cache, and never sees another test's entries.
func withTempCache(t *testing.T) {
	t.Helper()
	dir := t.TempDir()
	prev := runtimeCache
	runtimeCache = artifactcache.New(func() string { return filepath.Join(dir, artifactcache.FileName) })
	t.Cleanup(func() { runtimeCache = prev })
}

// requireSquashfsTools skips when the host cannot build or read a .sqf.
func requireSquashfsTools(t *testing.T) {
	t.Helper()
	for _, bin := range []string{"mksquashfs", "unsquashfs"} {
		if _, err := exec.LookPath(bin); err != nil {
			t.Skipf("%s not available", bin)
		}
	}
}

// marshalRuntime renders a runtime document the way a build stages it, for the
// tests that only care what DecodeRuntime makes of the bytes.
func marshalRuntime(t *testing.T, rt Runtime) []byte {
	t.Helper()
	data, err := MarshalRuntime(rt)
	if err != nil {
		t.Fatalf("MarshalRuntime: %v", err)
	}
	return data
}

// marshalManifest renders a manifest document, as marshalRuntime does.
func marshalManifest(t *testing.T, m Manifest) []byte {
	t.Helper()
	data, err := MarshalManifest(m)
	if err != nil {
		t.Fatalf("MarshalManifest: %v", err)
	}
	return data
}

// packSqf builds a .sqf whose root is dir.
func packSqf(t *testing.T, dir string) string {
	t.Helper()
	out := filepath.Join(t.TempDir(), "image.sqf")
	cmd := exec.Command("mksquashfs", dir, out, "-no-progress", "-noappend", "-quiet", "-no-xattrs")
	if output, err := cmd.CombinedOutput(); err != nil {
		t.Fatalf("mksquashfs: %v\n%s", err, output)
	}
	return out
}

// stagedImage stages both documents into a payload root and packs it, the way a
// build backend does.
func stagedImage(t *testing.T, rt Runtime, m Manifest) string {
	t.Helper()
	root := t.TempDir()
	dir := filepath.Join(root, DirName)
	if err := StageRuntime(dir, rt); err != nil {
		t.Fatalf("StageRuntime: %v", err)
	}
	if err := StageManifest(dir, m); err != nil {
		t.Fatalf("StageManifest: %v", err)
	}
	payload := filepath.Join(root, "cnt", "samtools", "1.23.1", "bin")
	if err := os.MkdirAll(payload, 0o755); err != nil {
		t.Fatalf("mkdir payload: %v", err)
	}
	if err := os.WriteFile(filepath.Join(payload, "samtools"), []byte("#!/bin/sh\n"), 0o755); err != nil {
		t.Fatalf("write payload: %v", err)
	}
	return packSqf(t, root)
}

func TestPrefix(t *testing.T) {
	tests := []struct {
		name string
		typ  catalog.Type
		want string
	}{
		{"samtools/1.23.1", catalog.TypeApp, "/cnt/samtools/1.23.1"},
		{"grch38/star/2.7.11b/gencode47", catalog.TypeData, "/cnt/grch38/star/2.7.11b/gencode47"},
		{"ubuntu24/base", catalog.TypeBase, ""},
		{"ubuntu24/r-essential", catalog.TypeOS, ""},
	}
	for _, tt := range tests {
		if got := Prefix(tt.name, tt.typ); got != tt.want {
			t.Errorf("Prefix(%q, %s) = %q, want %q", tt.name, tt.typ, got, tt.want)
		}
	}
}

func TestEnvVarResolved(t *testing.T) {
	v := EnvVar{Key: "PATH_EXTRA", Value: "{prefix}/bin:{prefix}/libexec"}
	if got, want := v.Resolved("/cnt/x/1"), "/cnt/x/1/bin:/cnt/x/1/libexec"; got != want {
		t.Errorf("Resolved() = %q, want %q", got, want)
	}
	// A value with no token survives untouched.
	plain := EnvVar{Key: "K", Value: "/absolute"}
	if got := plain.Resolved("/cnt/x/1"); got != "/absolute" {
		t.Errorf("Resolved() = %q, want /absolute", got)
	}
	// An os artifact's empty prefix is not a special case: its files really are
	// at the root, so {prefix}/bin is /bin, which is where they are.
	root := EnvVar{Key: "K", Value: "{prefix}/bin"}
	if got := root.Resolved(""); got != "/bin" {
		t.Errorf("Resolved() = %q, want /bin", got)
	}
}

// The recorded architecture is Go's spelling.
func TestNativeArchIsGoForm(t *testing.T) {
	switch got := NativeArch(); got {
	case "amd64", "arm64":
	default:
		t.Logf("unrecognized host architecture %q", got)
	}
	if NativePlatform().OS != "linux" {
		t.Errorf("platform OS = %q, want linux", NativePlatform().OS)
	}
}

func TestStagePayloadKeyRecordsTheKeyInTheStagedManifest(t *testing.T) {
	dir := t.TempDir()
	if err := StageManifest(dir, validManifest()); err != nil {
		t.Fatal(err)
	}
	ref := KeyRef{Scheme: "payload-tree-v1", SHA256: strings.Repeat("ab", 32)}
	if err := StagePayloadKey(dir, ref); err != nil {
		t.Fatal(err)
	}
	data, err := os.ReadFile(filepath.Join(dir, FileName))
	if err != nil {
		t.Fatal(err)
	}
	m, err := DecodeManifest(data, "staged")
	if err != nil {
		t.Fatal(err)
	}
	if m.Keys.Payload != ref {
		t.Errorf("payload key = %+v, want %+v", m.Keys.Payload, ref)
	}
}

func TestFrozenEnvironmentRefusesAPayloadKey(t *testing.T) {
	m := validManifest()
	m.Name, m.Type, m.BuildType = EnvName, catalog.TypeEnv, BuildTypeSnapshot
	m.Keys = Keys{Payload: KeyRef{Scheme: "payload-tree-v1", SHA256: strings.Repeat("ab", 32)}}
	if err := ValidateManifest(m); err == nil {
		t.Error("a frozen environment with a payload key validated")
	}
}
