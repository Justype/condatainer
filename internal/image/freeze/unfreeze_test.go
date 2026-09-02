package freeze

import (
	"context"
	"errors"
	"os"
	"os/exec"
	"path/filepath"
	"strconv"
	"strings"
	"testing"

	"github.com/Justype/condatainer/internal/artifact/meta"
)

// unfreezeInto rebuilds a writable overlay from an artifact.
func unfreezeInto(t *testing.T, sqf string, opts UnfreezeOptions) UnfreezeResult {
	t.Helper()
	if _, err := exec.LookPath("unsquashfs"); err != nil {
		t.Skip("unsquashfs not available")
	}
	opts.Artifact = sqf
	if opts.Target == "" {
		opts.Target = filepath.Join(artifactDir(t), "env.img")
	}
	if opts.UID == 0 && opts.GID == 0 {
		opts.UID, opts.GID = os.Getuid(), os.Getgid()
	}
	opts.Sparse = true
	if opts.Base == "" {
		opts.Base = repoBase(t)
	}
	res, err := Unfreeze(context.Background(), opts)
	if err != nil {
		t.Fatalf("unfreeze: %v", err)
	}
	return res
}

// imgEntry returns the debugfs listing fields for one path inside an image.
func imgEntry(t *testing.T, img, dir, name string) (mode string, found bool) {
	t.Helper()
	for _, line := range strings.Split(dbgLs(t, img, dir), "\n") {
		f := strings.Split(strings.ReplaceAll(strings.TrimSpace(line), " ", ""), "/")
		if len(f) >= 7 && f[5] == name {
			return f[2], true
		}
	}
	return "", false
}

func dbgLs(t *testing.T, img, dir string) string {
	t.Helper()
	out, err := exec.Command("debugfs", "-R", "ls -p "+dir, img).Output()
	if err != nil {
		t.Fatalf("debugfs ls %s: %v", dir, err)
	}
	return string(out)
}

// The round trip puts the payload back where it was, under upper/, with work/
// beside it and the metadata directory gone — a writable overlay carries no
// embedded metadata by design.
func TestUnfreezeRoundTrip(t *testing.T) {
	t.Parallel()
	img := newImage(t, `cd upper
mkdir cnt_env
cd cnt_env
mkdir bin`)
	writeInto(t, img, "/upper/cnt_env/bin", "tool", "payload\n")

	res := unfreezeInto(t, freezeImage(t, img).Path, UnfreezeOptions{})
	if res.From != meta.EnvName {
		t.Errorf("From = %q", res.From)
	}
	for _, dir := range []string{"/upper", "/work", "/work/work", "/upper/cnt_env/bin"} {
		if strings.Contains(dbgLs(t, res.Path, dir), "File not found") {
			t.Errorf("%s missing from the rebuilt overlay", dir)
		}
	}
	if _, found := imgEntry(t, res.Path, "/upper/cnt_env/bin", "tool"); !found {
		t.Error("the payload did not come back")
	}
	if _, found := imgEntry(t, res.Path, "/upper", meta.DirName); found {
		t.Error(".cnt came back into upper/; a writable overlay carries no embedded metadata")
	}
}

// Deletions survive because the image is built straight from the mounted
// artifact: mke2fs writes filesystem metadata rather than calling mknod, so it
// records a character device the caller could not have created. Extraction
// cannot — unsquashfs drops every one and says so only in its summary — which is
// why the payload is never extracted.
func TestUnfreezeKeepsWhiteouts(t *testing.T) {
	t.Parallel()
	img := newImage(t, `cd upper
mkdir cnt_env
cd /upper
mkdir usr
cd usr
mkdir bin`)
	writeInto(t, img, "/upper/cnt_env", "f", "x\n")
	writeInto(t, img, "/upper/usr/bin", ".wh.micromamba", "")

	res := unfreezeInto(t, freezeImage(t, img).Path, UnfreezeOptions{})
	mode, found := imgEntry(t, res.Path, "/upper/usr/bin", "micromamba")
	if !found {
		t.Fatal("the whiteout did not survive the round trip; the deleted file is back")
	}
	if !strings.HasPrefix(mode, "02") {
		t.Errorf("replayed as mode %s, want a char device", mode)
	}
}

// An archive is packed -all-root, so the payload has to arrive owned by the
// user or the overlay mounts and then denies every write into it. squashfuse is
// told the uid and gid, which is why no walk of the image is needed afterwards.
func TestUnfreezeOwnsThePayload(t *testing.T) {
	t.Parallel()
	img := newImage(t, "cd upper\nmkdir cnt_env")
	writeInto(t, img, "/upper/cnt_env", "f", "x\n")

	res := unfreezeInto(t, freezeImage(t, img).Path, UnfreezeOptions{})
	out, err := exec.Command("debugfs", "-R", "ls -l /upper/cnt_env", res.Path).Output()
	if err != nil {
		t.Fatal(err)
	}
	for _, line := range strings.Split(string(out), "\n") {
		f := strings.Fields(line)
		if len(f) < 4 || f[len(f)-1] == "." || f[len(f)-1] == ".." {
			continue
		}
		if f[3] != strconv.Itoa(os.Getuid()) {
			t.Errorf("payload entry %q is owned by uid %s, not %d", f[len(f)-1], f[3], os.Getuid())
		}
	}
}

// A size the payload does not fit in is refused before mke2fs runs, because
// mke2fs -d accepts one and only discovers the problem while writing.
func TestUnfreezeRefusesAnUndersizedImage(t *testing.T) {
	t.Parallel()
	img := newImage(t, "cd upper\nmkdir cnt_env")
	writeInto(t, img, "/upper/cnt_env", "f", "x\n")
	sqf := freezeImage(t, img).Path

	target := filepath.Join(artifactDir(t), "env.img")
	_, err := Unfreeze(context.Background(), UnfreezeOptions{
		Artifact: sqf, Target: target, SizeMB: 1, UID: os.Getuid(), GID: os.Getgid(),
		Sparse: true, Base: repoBase(t),
	})
	if !errors.Is(err, ErrTooSmall) {
		t.Fatalf("error = %v, want ErrTooSmall", err)
	}
	if _, statErr := os.Stat(target); statErr == nil {
		t.Error("a half-built image was left behind")
	}
}

// Unfreeze is the inverse of freeze and takes what freeze produced.
func TestUnfreezeRefusesANonSnapshot(t *testing.T) {
	t.Parallel()
	_, err := Unfreeze(context.Background(), UnfreezeOptions{
		Artifact: repoAppArtifact(t), Target: filepath.Join(artifactDir(t), "env.img"),
		UID: os.Getuid(), GID: os.Getgid(), Sparse: true, Base: repoBase(t),
	})
	if !errors.Is(err, ErrNotFrozen) {
		t.Fatalf("error = %v, want ErrNotFrozen", err)
	}
}

// repoAppArtifact is a recipe-built artifact the repository ships.
func repoAppArtifact(t *testing.T) string {
	t.Helper()
	wd, err := os.Getwd()
	if err != nil {
		t.Fatal(err)
	}
	p := filepath.Join(wd, "..", "..", "..", "images", "samtools--1.23.1.sqf")
	if _, err := os.Stat(p); err != nil {
		t.Skipf("no recipe-built artifact to test against: %v", err)
	}
	return p
}

// The sidecar says what this overlay came from, as a comment: the parser skips
// those, so provenance cannot be mistaken for an environment variable.
func TestUnfreezeWritesTheSidecar(t *testing.T) {
	t.Parallel()
	img := newImage(t, "cd upper\nmkdir cnt_env")
	writeInto(t, img, "/upper/cnt_env", "f", "x\n")

	res := unfreezeInto(t, freezeImage(t, img).Path, UnfreezeOptions{})
	data, err := os.ReadFile(res.Path + ".env")
	if err != nil {
		t.Fatalf("no sidecar: %v", err)
	}
	body := string(data)
	if !strings.Contains(body, meta.EnvName) {
		t.Errorf("sidecar does not name the artifact: %q", body)
	}
	for _, line := range strings.Split(strings.TrimSpace(body), "\n") {
		if !strings.HasPrefix(strings.TrimSpace(line), "#") {
			t.Errorf("sidecar line is not a comment and would parse as a variable: %q", line)
		}
	}
}
