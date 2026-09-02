package freeze

import (
	"context"
	"os"
	"os/exec"
	"path/filepath"
	"strings"
	"testing"

	"github.com/Justype/condatainer/internal/utils"
)

// listArchiveLong returns one line per entry: type, mode and path, with the sizes
// and timestamps that legitimately differ between two packs dropped. It is what
// "the same archive" means here — the same things, in the same shapes.
func listArchiveLong(t *testing.T, sqf string) []string {
	t.Helper()
	out := runTool(t, "unsquashfs", "-ll", sqf)
	var lines []string
	for _, l := range strings.Split(out, "\n") {
		fields := strings.Fields(l)
		if len(fields) < 4 || !strings.Contains(l, "squashfs-root") {
			continue
		}
		idx := strings.Index(l, "squashfs-root")
		path := strings.TrimPrefix(l[idx:], "squashfs-root")
		if strings.HasPrefix(path, "/.cnt") {
			continue // the manifest records a build time, so it differs by design
		}
		lines = append(lines, fields[0]+" "+path)
	}
	return lines
}

func runTool(t *testing.T, name string, args ...string) string {
	t.Helper()
	out, err := exec.Command(name, args...).CombinedOutput()
	if err != nil {
		t.Fatalf("%s: %v\n%s", name, err, out)
	}
	return string(out)
}

// freezeVia freezes img by one route and returns the artifact.
func freezeVia(t *testing.T, img string, useTmp bool) Result {
	t.Helper()
	target := filepath.Join(artifactDir(t), "frozen.sqf")
	res, err := Freeze(context.Background(), Options{
		Image: img, Target: target, Base: repoBase(t),
		CompressArgs: "-comp zstd", UseTmp: useTmp,
	})
	if err != nil {
		t.Fatalf("freeze (useTmp=%v): %v", useTmp, err)
	}
	return res
}

// The two routes must produce the same archive. They read the payload by
// completely different means — a FUSE mount of the image, or a staged rdump copy
// — and rdump cannot create a device node, so without ForCopy the copy route
// drops every char 0:0 whiteout silently. This is the test that catches that.
func TestBothRoutesProduceTheSameArchive(t *testing.T) {
	t.Setenv("CNT_TMPDIR", t.TempDir())
	img := newImage(t, `cd upper
mkdir cnt_env
cd cnt_env
mkdir bin
cd /upper
mkdir usr
cd usr
mkdir bin
cd bin
mknod deleted c 0 0`)
	writeInto(t, img, "/upper/cnt_env/bin", "tool", "payload\n")
	writeInto(t, img, "/upper/usr/bin", ".wh.gone", "")

	mounted := freezeVia(t, img, false)
	copied := freezeVia(t, img, true)

	want := listArchiveLong(t, mounted.Path)
	got := listArchiveLong(t, copied.Path)
	if strings.Join(want, "\n") != strings.Join(got, "\n") {
		t.Errorf("the routes disagree.\nmount route:\n%s\n\ncopy route:\n%s",
			strings.Join(want, "\n"), strings.Join(got, "\n"))
	}

	// Both whiteouts have to be there, whichever route produced them: the one
	// translated from a marker file and the one the image already held as a node.
	var devices int
	for _, l := range want {
		if strings.HasPrefix(l, "c") {
			devices++
		}
	}
	if devices != 2 {
		t.Errorf("archive carries %d char devices, want 2 (one per deletion):\n%s",
			devices, strings.Join(want, "\n"))
	}
	if mounted.Translation.Deletions() != copied.Translation.Deletions() {
		t.Errorf("deletion counts differ: mount=%d copy=%d",
			mounted.Translation.Deletions(), copied.Translation.Deletions())
	}

	// Neither route may leave its working copy behind: the copy route stages the
	// whole payload, and an abandoned stage is a second copy of an environment
	// sitting in someone's scratch.
	left, err := os.ReadDir(utils.GetTmpDir())
	if err != nil {
		t.Fatal(err)
	}
	for _, e := range left {
		t.Errorf("left behind in scratch: %s", e.Name())
	}
}

// Scratch that cannot hold the payload is refused before anything is copied.
func TestCopyRouteRefusesShortScratch(t *testing.T) {
	if err := checkStageSpace(t.TempDir(), 1<<30); err == nil {
		t.Fatal("a petabyte-sized payload was accepted")
	}
}
