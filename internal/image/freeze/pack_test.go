package freeze

import (
	"context"
	"fmt"
	"os"
	"os/exec"
	"path/filepath"
	"strings"
	"testing"
	"time"

	"github.com/Justype/condatainer/internal/config"
)

// repoBase is the base image the repository ships for its own tests, used
// only where a test exercises opaque-directory translation against a real
// base's contents (see BaseDirLister) — Pack and Unfreeze need no base at all.
func repoBase(t *testing.T) string {
	t.Helper()
	wd, err := os.Getwd()
	if err != nil {
		t.Fatal(err)
	}
	base := filepath.Join(wd, "..", "..", "..", "images", "ubuntu24--base.sqf")
	if _, err := os.Stat(base); err != nil {
		t.Skipf("base image not built: %v", err)
	}
	return base
}

// artifactDir is a directory for pack output, cleaned up with a retry.
//
// Apptainer's FUSE helper can still hold the artifact open when Run returns, and
// on a shared filesystem unlinking a file somebody has open leaves a .nfs*
// placeholder behind, which makes a plain RemoveAll fail. Retry briefly and then
// give up quietly: a leftover temporary directory is not worth a failed test, and
// t.TempDir has no way to express "try again in a moment".
func artifactDir(t *testing.T) string {
	t.Helper()
	dir, err := os.MkdirTemp("", "cnt-freeze-test-")
	if err != nil {
		t.Fatal(err)
	}
	t.Cleanup(func() {
		for range 20 {
			if err := os.RemoveAll(dir); err == nil {
				return
			}
			time.Sleep(100 * time.Millisecond)
		}
	})
	return dir
}

// packImage freezes img and returns the artifact path.
func packImage(t *testing.T, img string) string {
	t.Helper()
	ctx := context.Background()

	entries, err := Walk(ctx, img)
	if err != nil {
		t.Fatal(err)
	}
	tr, err := Translate(ctx, img, entries, noBase)
	if err != nil {
		t.Fatal(err)
	}
	target := filepath.Join(artifactDir(t), "frozen.sqf")
	if _, err := Pack(ctx, PackOptions{
		Image: img, Target: target, CompressArgs: "-comp zstd",
	}, entries, tr); err != nil {
		t.Fatalf("pack: %v", err)
	}
	return target
}

// listArchive returns the archive's paths, with unsquashfs's display root removed.
func listArchive(t *testing.T, sqf string) []string {
	t.Helper()
	if _, err := exec.LookPath("unsquashfs"); err != nil {
		t.Skip("unsquashfs not available")
	}
	out, err := exec.Command("unsquashfs", "-l", sqf).CombinedOutput()
	if err != nil {
		t.Fatalf("unsquashfs: %v\n%s", err, out)
	}
	var paths []string
	for _, line := range strings.Split(string(out), "\n") {
		if p, ok := strings.CutPrefix(strings.TrimSpace(line), "squashfs-root/"); ok && p != "" {
			paths = append(paths, p)
		}
	}
	return paths
}

// The archive holds the payload with the upper/ wrapper stripped, and holds no
// base content — the regression being a pack from the merged view rather than
// from the image's own upper layer.
func TestPackStripsUpperNotTheMergedView(t *testing.T) {
	requireSquashfs(t)
	t.Parallel()
	img := newImage(t, `cd upper
mkdir cnt_env
cd cnt_env
mkdir bin
cd /upper
mkdir opt`)
	writeInto(t, img, "/upper/cnt_env/bin", "tool", "payload\n")

	paths := listArchive(t, packImage(t, img))
	for _, want := range []string{"cnt_env", "cnt_env/bin", "cnt_env/bin/tool", "opt"} {
		if !contains(paths, want) {
			t.Errorf("archive is missing %q; has %v", want, paths)
		}
	}
	for _, leaked := range []string{"bin", "lib", "etc", "usr", "var"} {
		if contains(paths, leaked) {
			t.Errorf("archive contains base directory %q: the pack read the merged view", leaked)
		}
	}
}

// A .wh. marker becomes a char 0:0 node in the archive and the marker itself is
// not packed — whiteout translation, end to end.
func TestPackTranslatesWhiteouts(t *testing.T) {
	requireSquashfs(t)
	t.Parallel()
	img := newImage(t, `cd upper
mkdir usr
cd usr
mkdir bin`)
	writeInto(t, img, "/upper/usr/bin", ".wh.micromamba", "")
	sqf := packImage(t, img)

	if contains(listArchive(t, sqf), "usr/bin/.wh.micromamba") {
		t.Error("the marker file was packed; it would surface as an ordinary file at mount")
	}
	out, err := exec.Command("unsquashfs", "-ll", sqf).CombinedOutput()
	if err != nil {
		t.Fatalf("unsquashfs -ll: %v\n%s", err, out)
	}
	var line string
	for _, l := range strings.Split(string(out), "\n") {
		if strings.HasSuffix(strings.TrimSpace(l), "usr/bin/micromamba") {
			line = strings.TrimSpace(l)
		}
	}
	if line == "" {
		t.Fatalf("no whiteout for the deleted file:\n%s", out)
	}
	if !strings.HasPrefix(line, "c") || !strings.Contains(line, "0,   0") && !strings.Contains(line, "0,  0") {
		t.Errorf("whiteout is not a char 0:0 node: %q", line)
	}
}

// The source is byte-identical afterwards, and its permissions are restored.
func TestPackLeavesTheSourceUntouched(t *testing.T) {
	requireSquashfs(t)
	t.Parallel()
	img := newImage(t, "cd upper\nmkdir cnt_env")
	writeInto(t, img, "/upper/cnt_env", "f", "x\n")

	before, err := os.ReadFile(img)
	if err != nil {
		t.Fatal(err)
	}
	info, err := os.Stat(img)
	if err != nil {
		t.Fatal(err)
	}
	packImage(t, img)

	after, err := os.ReadFile(img)
	if err != nil {
		t.Fatal(err)
	}
	if string(before) != string(after) {
		t.Error("the pack modified the source image")
	}
	now, err := os.Stat(img)
	if err != nil {
		t.Fatal(err)
	}
	if now.Mode().Perm() != info.Mode().Perm() {
		t.Errorf("permissions left as %v, were %v", now.Mode().Perm(), info.Mode().Perm())
	}
}

// An image whose write bit is already clear is packed as it stands and stays
// pinned: clearing that bit is how an artifact is protected, and freeze must not
// quietly unpin one.
func TestPackKeepsAProtectedImageProtected(t *testing.T) {
	requireSquashfs(t)
	t.Parallel()
	img := newImage(t, "cd upper\nmkdir cnt_env")
	writeInto(t, img, "/upper/cnt_env", "f", "x\n")
	if err := os.Chmod(img, 0o444); err != nil {
		t.Fatal(err)
	}
	packImage(t, img)

	info, err := os.Stat(img)
	if err != nil {
		t.Fatal(err)
	}
	if info.Mode().Perm()&0o222 != 0 {
		t.Errorf("the pack made a protected image writable: %v", info.Mode().Perm())
	}
}

// An overlay nobody has written to has nothing to freeze, and says so rather
// than producing an empty artifact.
func TestPackRefusesAnEmptyOverlay(t *testing.T) {
	requireSquashfs(t)
	t.Parallel()
	img := newImage(t, "")
	ctx := context.Background()
	entries, err := Walk(ctx, img)
	if err != nil {
		t.Fatal(err)
	}
	_, err = Pack(ctx, PackOptions{
		Image: img, Target: filepath.Join(artifactDir(t), "x.sqf"),
	}, entries, Translation{})
	if err == nil {
		t.Fatal("an empty overlay was packed")
	}
	if !strings.Contains(err.Error(), "nothing to freeze") {
		t.Errorf("unhelpful refusal: %v", err)
	}
}

// A pack always names its cpu budget: without -processors mksquashfs takes every
// core on the machine, so an unset one falls back to the build default.
func TestPackScriptAlwaysBudgetsCpus(t *testing.T) {
	t.Parallel()
	for _, tc := range []struct {
		name  string
		given int
		want  string
	}{
		{"unset falls back to the build default", 0, fmt.Sprintf("-processors %d", config.DefaultNcpus)},
		{"negative falls back too", -1, fmt.Sprintf("-processors %d", config.DefaultNcpus)},
		{"a budget is honoured", 8, "-processors 8"},
	} {
		t.Run(tc.name, func(t *testing.T) {
			script := packScript([]string{"/src"}, "/out.sqf", nil, PackOptions{Processors: tc.given}, "mksquashfs")
			if !strings.Contains(script, tc.want) {
				t.Errorf("packScript(Processors=%d) = %q, want it to contain %q", tc.given, script, tc.want)
			}
		})
	}
}

// mksquashfs announces a collision and still exits 0, so the announcement is
// what AppendMeta has to fail on.
func TestCollidedReadsMksquashfsAnnouncement(t *testing.T) {
	const said = "Source directory entry .cnt already used! - trying .cnt_1"
	if !collided(said) {
		t.Errorf("collided(%q) = false, want true", said)
	}
	if collided("Number of directories 4\nNumber of hard-links 0") {
		t.Error("ordinary mksquashfs output read as a collision")
	}
}
