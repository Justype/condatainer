package freeze

import (
	"context"
	"errors"
	"os"
	"path/filepath"
	"testing"

	"github.com/Justype/condatainer/catalog"
	"github.com/Justype/condatainer/internal/artifact/key"
	"github.com/Justype/condatainer/internal/artifact/meta"
)

// readManifest pulls the embedded manifest back out of a frozen artifact,
// through the same reader the rest of the tool uses.
func readManifest(t *testing.T, sqf string) meta.Manifest {
	t.Helper()
	m, err := meta.ReadManifest(sqf)
	if err != nil {
		t.Fatalf("reading the embedded manifest: %v", err)
	}
	return m
}

func freezeImage(t *testing.T, img string) Result {
	t.Helper()
	res, err := Freeze(context.Background(), Options{
		Image: img, Target: filepath.Join(artifactDir(t), "frozen.sqf"),
		Base: repoBase(t), CompressArgs: "-comp zstd",
	})
	if err != nil {
		t.Fatalf("freeze: %v", err)
	}
	return res
}

// The artifact carries a manifest that says what it is and what it contributes.
// The manifest cannot describe the archive it sits inside, so it is packed in a
// second pass once the payload has been identified.
func TestFreezeEmbedsMetadata(t *testing.T) {
	requireSquashfs(t)
	t.Parallel()
	img := newImage(t, `cd upper
mkdir cnt_env
cd /upper
mkdir opt`)
	writeInto(t, img, "/upper/cnt_env", "f", "payload\n")

	res := freezeImage(t, img)
	m := readManifest(t, res.Path)

	if m.Name != meta.EnvName {
		t.Errorf("name = %q, want %q", m.Name, meta.EnvName)
	}
	if m.Type != catalog.TypeEnv {
		t.Errorf("type = %q, want %q", m.Type, catalog.TypeEnv)
	}
	if m.BuildType != meta.BuildTypeSnapshot {
		t.Errorf("build_type = %q, want %q", m.BuildType, meta.BuildTypeSnapshot)
	}
	if m.Keys.Identity.Empty() {
		t.Error("a frozen environment embedded no identity key")
	}
	if m.Snapshot == nil {
		t.Fatal("no snapshot block")
	}
	if err := meta.ValidateManifest(m); err != nil {
		t.Errorf("the embedded manifest does not validate: %v", err)
	}
}

// The identity is hashed from the payload and embedded, so the artifact answers
// for itself without anyone re-deriving it.
func TestFreezeEmbedsTheIdentity(t *testing.T) {
	requireSquashfs(t)
	t.Parallel()
	img := newImage(t, "cd upper\nmkdir cnt_env")
	writeInto(t, img, "/upper/cnt_env", "f", "payload\n")

	res := freezeImage(t, img)
	if res.Identity.Scheme != string(key.PayloadTreeV1) {
		t.Errorf("scheme = %q, want %q", res.Identity.Scheme, key.PayloadTreeV1)
	}
	if len(res.Identity.SHA256) != 64 {
		t.Errorf("sha256 = %q", res.Identity.SHA256)
	}
	// Equivalence is the same value: a snapshot has no inputs to abstract away.
	if res.Manifest.Keys.Identity != res.Identity || res.Manifest.Keys.Equiv != res.Identity {
		t.Errorf("manifest keys = %+v, want both %v", res.Manifest.Keys, res.Identity)
	}

	// What the file says, read back the way every other artifact is read.
	got, err := meta.ReadManifest(res.Path)
	if err != nil {
		t.Fatalf("read back the embedded manifest: %v", err)
	}
	if got.Keys.Identity != res.Identity {
		t.Errorf("embedded identity = %v, want %v", got.Keys.Identity, res.Identity)
	}
}

// The identity is the payload's, not the file's: repacking the same overlay
// yields a different file and the same identity. That is the whole reason it is
// hashed from the tree rather than from the bytes.
func TestFreezeIdentitySurvivesARepack(t *testing.T) {
	requireSquashfs(t)
	t.Parallel()
	img := newImage(t, "cd upper\nmkdir cnt_env")
	writeInto(t, img, "/upper/cnt_env", "f", "payload\n")

	first := freezeImage(t, img)
	second := freezeImage(t, img)
	if first.Identity != second.Identity {
		t.Errorf("two freezes of one overlay disagreed: %v vs %v", first.Identity, second.Identity)
	}
}

// The recorded convention is what the overlay used before translation, and the
// count is what the artifact carries after it.
func TestFreezeRecordsTheConvention(t *testing.T) {
	requireSquashfs(t)
	t.Parallel()
	img := newImage(t, `cd upper
mkdir cnt_env
cd /upper
mkdir usr
cd usr
mkdir bin`)
	writeInto(t, img, "/upper/cnt_env", "f", "payload\n")
	writeInto(t, img, "/upper/usr/bin", ".wh.micromamba", "")

	m := readManifest(t, freezeImage(t, img).Path)
	if m.Snapshot.Convention != string(ConventionWhiteFile) {
		t.Errorf("convention = %q, want %q", m.Snapshot.Convention, ConventionWhiteFile)
	}
	if m.Snapshot.Whiteouts != 1 {
		t.Errorf("whiteouts = %d, want 1", m.Snapshot.Whiteouts)
	}
}

// An environment records EnvPrefix when its payload has one, so PATH resolves as
// it did for the .img — and nothing when it does not, rather than claiming a path
// it cannot provide.
func TestFreezeRecordsPrefixOnlyWhenPresent(t *testing.T) {
	requireSquashfs(t)
	t.Parallel()
	with := newImage(t, "cd upper\nmkdir cnt_env")
	writeInto(t, with, "/upper/cnt_env", "f", "x\n")
	without := newImage(t, "cd upper\nmkdir opt")
	writeInto(t, without, "/upper/opt", "f", "x\n")

	for _, tc := range []struct {
		name string
		img  string
		want string
	}{
		{"with a conda prefix", with, meta.EnvPrefix},
		{"without one", without, ""},
	} {
		t.Run(tc.name, func(t *testing.T) {
			rt, err := meta.ReadRuntime(freezeImage(t, tc.img).Path)
			if err != nil {
				t.Fatal(err)
			}
			if rt.Prefix != tc.want {
				t.Errorf("prefix = %q, want %q", rt.Prefix, tc.want)
			}
		})
	}
}

// An overlay nobody has written to has nothing to freeze.
func TestFreezeRefusesAnEmptyOverlay(t *testing.T) {
	t.Parallel()
	_, err := Freeze(context.Background(), Options{
		Image: newImage(t, ""), Target: filepath.Join(artifactDir(t), "x.sqf"),
	})
	if !errors.Is(err, ErrEmptyOverlay) {
		t.Fatalf("error = %v, want ErrEmptyOverlay", err)
	}
}

// The architecture is recorded in uname form, which is what MountAllowed
// compares against. GOARCH here would make every frozen artifact look as though
// it were built for another machine, and the contribution — roots included —
// would be dropped at mount with only a warning.
func TestFreezeRecordsTheNativeArch(t *testing.T) {
	requireSquashfs(t)
	t.Parallel()
	img := newImage(t, "cd upper\nmkdir cnt_env")
	writeInto(t, img, "/upper/cnt_env", "f", "x\n")

	m := readManifest(t, freezeImage(t, img).Path)
	if m.Platform.Arch != meta.NativeArch() {
		t.Errorf("arch = %q, want %q", m.Platform.Arch, meta.NativeArch())
	}
	rt, err := meta.ReadRuntime(freezeImage(t, img).Path)
	if err != nil {
		t.Fatal(err)
	}
	if rt.Platform.Arch != meta.NativeArch() {
		t.Errorf("runtime arch = %q, want %q", rt.Platform.Arch, meta.NativeArch())
	}
}

// A writable overlay keeps its environment in a sidecar; a frozen one cannot,
// because the mount path ignores a sidecar beside a .sqf. The variables have to
// move inside the artifact or they are lost at the freeze.
func TestFreezeCarriesTheSidecar(t *testing.T) {
	requireSquashfs(t)
	t.Parallel()
	img := newImage(t, "cd upper\nmkdir cnt_env")
	writeInto(t, img, "/upper/cnt_env", "f", "x\n")
	sidecar := "GOROOT={prefix}/go  ## Go installation\nGOPATH=/home/go\n# a comment\n"
	if err := os.WriteFile(img+".env", []byte(sidecar), 0o644); err != nil {
		t.Fatal(err)
	}

	rt, err := meta.ReadRuntime(freezeImage(t, img).Path)
	if err != nil {
		t.Fatal(err)
	}
	got := map[string]meta.EnvVar{}
	for _, v := range rt.Env {
		got[v.Key] = v
	}
	if len(got) != 2 {
		t.Fatalf("carried %d variables, want 2: %+v", len(got), rt.Env)
	}
	// {prefix} stays a token: it is substituted when the image is loaded, which
	// is the same contract every other artifact's environment has.
	if got["GOROOT"].Value != "{prefix}/go" {
		t.Errorf("GOROOT = %q, want the token kept", got["GOROOT"].Value)
	}
	if got["GOROOT"].Note != "Go installation" {
		t.Errorf("note = %q", got["GOROOT"].Note)
	}
	if got["GOPATH"].Value != "/home/go" {
		t.Errorf("GOPATH = %q", got["GOPATH"].Value)
	}
}
