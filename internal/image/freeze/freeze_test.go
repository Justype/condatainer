package freeze

import (
	"context"
	"fmt"
	"os"
	"os/exec"
	"path/filepath"
	"strings"
	"testing"

	"github.com/Justype/condatainer/internal/toolpath"
)

// requireTools skips when the ext3 toolchain is absent, which is the case on a
// machine that could not run a freeze either.
func requireTools(t *testing.T) {
	t.Helper()
	for _, bin := range []string{"debugfs", "mke2fs", "dd"} {
		if _, err := exec.LookPath(bin); err != nil {
			t.Skipf("%s not available", bin)
		}
	}
}

// requireSquashfs skips when a pack could not run: mksquashfs as the freeze
// resolves it, and unsquashfs off PATH for the tests that read the archive back.
func requireSquashfs(t *testing.T) {
	t.Helper()
	if _, err := toolpath.Resolve("mksquashfs"); err != nil {
		t.Skipf("mksquashfs not available: %v", err)
	}
	if _, err := exec.LookPath("unsquashfs"); err != nil {
		t.Skip("unsquashfs not available")
	}
}

// newImage builds an overlay image and runs script through debugfs to populate
// upper/. Building the payload with debugfs rather than through a container
// keeps the test hermetic: no apptainer, no base image, no privilege.
func newImage(t *testing.T, script string) string {
	t.Helper()
	requireTools(t)
	path := filepath.Join(t.TempDir(), "overlay.img")

	run := func(name string, args ...string) {
		t.Helper()
		if out, err := exec.Command(name, args...).CombinedOutput(); err != nil {
			t.Fatalf("%s: %v\n%s", name, err, out)
		}
	}
	run("dd", "if=/dev/zero", "of="+path, "bs=1M", "count=0", "seek=16", "status=none")
	run("mke2fs", "-q", "-t", "ext3", "-F", path)

	cmd := exec.Command("debugfs", "-w", path)
	cmd.Stdin = strings.NewReader("mkdir upper\nmkdir work\n" + script + "\nquit\n")
	if out, err := cmd.CombinedOutput(); err != nil {
		t.Fatalf("debugfs: %v\n%s", err, out)
	}
	return path
}

// writeInto stages a file and links it into the image at dir.
func writeInto(t *testing.T, img, dir, name, content string) {
	t.Helper()
	src := filepath.Join(t.TempDir(), name)
	if err := os.WriteFile(src, []byte(content), 0o644); err != nil {
		t.Fatal(err)
	}
	cmd := exec.Command("debugfs", "-w", img)
	cmd.Stdin = strings.NewReader(fmt.Sprintf("cd %s\nwrite %s %s\nquit\n", dir, src, name))
	if out, err := cmd.CombinedOutput(); err != nil {
		t.Fatalf("write %s: %v\n%s", name, err, out)
	}
}

func paths(entries []Entry) []string {
	out := make([]string, 0, len(entries))
	for _, e := range entries {
		out = append(out, e.Path)
	}
	return out
}

func contains(list []string, want string) bool {
	for _, s := range list {
		if s == want {
			return true
		}
	}
	return false
}

// The walk descends the whole tree and reports paths relative to upper/, which
// is what they will be in the archive once the wrapper is stripped.
func TestWalkFindsTheWholeTree(t *testing.T) {
	t.Parallel()
	img := newImage(t, `cd upper
mkdir cnt_env
cd cnt_env
mkdir bin
cd /upper
mkdir opt`)
	writeInto(t, img, "/upper/cnt_env/bin", "tool", "payload\n")

	entries, err := Walk(context.Background(), img)
	if err != nil {
		t.Fatal(err)
	}
	for _, want := range []string{"cnt_env", "cnt_env/bin", "cnt_env/bin/tool", "opt"} {
		if !contains(paths(entries), want) {
			t.Errorf("walk missed %q; found %v", want, paths(entries))
		}
	}
	for _, e := range entries {
		if e.Path == "cnt_env/bin/tool" {
			if e.IsDir() || e.IsCharDev() {
				t.Errorf("tool reported as dir=%v chardev=%v", e.IsDir(), e.IsCharDev())
			}
			if e.Size != int64(len("payload\n")) {
				t.Errorf("size = %d, want %d", e.Size, len("payload\n"))
			}
		}
	}
}

// The walk never touches the image: enumeration has to be safe against an
// overlay whose write bit has not been cleared.
func TestWalkDoesNotModifyTheImage(t *testing.T) {
	t.Parallel()
	img := newImage(t, "cd upper\nmkdir cnt_env")
	before, err := os.ReadFile(img)
	if err != nil {
		t.Fatal(err)
	}
	if _, err := Walk(context.Background(), img); err != nil {
		t.Fatal(err)
	}
	after, err := os.ReadFile(img)
	if err != nil {
		t.Fatal(err)
	}
	if string(before) != string(after) {
		t.Fatal("walking the image changed its bytes")
	}
}

func noBase(context.Context, []string) (map[string][]string, error) { return nil, nil }

// fuse-overlayfs records a deletion as a plain .wh.<name> file, which is inert in
// a read-only layer. It must become a char 0:0 node, and the marker must not be
// packed — it would otherwise surface in the merged view as an ordinary file.
func TestTranslateWhiteoutFile(t *testing.T) {
	t.Parallel()
	img := newImage(t, `cd upper
mkdir usr
cd usr
mkdir bin`)
	writeInto(t, img, "/upper/usr/bin", ".wh.micromamba", "")

	entries, err := Walk(context.Background(), img)
	if err != nil {
		t.Fatal(err)
	}
	tr, err := Translate(context.Background(), img, entries, noBase)
	if err != nil {
		t.Fatal(err)
	}
	if tr.Convention != ConventionWhiteFile {
		t.Errorf("convention = %q, want %q", tr.Convention, ConventionWhiteFile)
	}
	if !contains(tr.Exclude, "usr/bin/.wh.micromamba") {
		t.Errorf("marker not excluded: %v", tr.Exclude)
	}
	if !contains(tr.Pseudo, "usr/bin/micromamba c 0 0 0 0 0") {
		t.Errorf("no whiteout emitted for the deleted file: %v", tr.Pseudo)
	}
}

// A char 0:0 node is already the form the artifact needs, so it passes through:
// not excluded, and not duplicated as a pseudo definition.
func TestTranslateCharDevPassesThrough(t *testing.T) {
	t.Parallel()
	img := newImage(t, `cd upper
mkdir usr
cd usr
mkdir bin
cd bin
mknod micromamba c 0 0`)

	entries, err := Walk(context.Background(), img)
	if err != nil {
		t.Fatal(err)
	}
	tr, err := Translate(context.Background(), img, entries, noBase)
	if err != nil {
		t.Fatal(err)
	}
	if tr.Convention != ConventionCharDev {
		t.Errorf("convention = %q, want %q", tr.Convention, ConventionCharDev)
	}
	if len(tr.Exclude) != 0 || len(tr.Pseudo) != 0 {
		t.Errorf("a char device was rewritten: exclude=%v pseudo=%v", tr.Exclude, tr.Pseudo)
	}
}

// An opaque directory hides everything the base has at that path, which has no
// pseudo-file representation — so it becomes one whiteout per base entry, minus
// whatever the overlay provides itself.
func TestTranslateOpaqueMarkerFile(t *testing.T) {
	t.Parallel()
	img := newImage(t, `cd upper
mkdir etc`)
	writeInto(t, img, "/upper/etc", opqMarker, "")
	writeInto(t, img, "/upper/etc", "mine.conf", "kept\n")

	base := func(_ context.Context, dirs []string) (map[string][]string, error) {
		if len(dirs) != 1 || dirs[0] != "/etc" {
			t.Errorf("base listed %v, want [/etc]", dirs)
		}
		return map[string][]string{"/etc": {"hosts", "passwd", "mine.conf"}}, nil
	}

	entries, err := Walk(context.Background(), img)
	if err != nil {
		t.Fatal(err)
	}
	tr, err := Translate(context.Background(), img, entries, base)
	if err != nil {
		t.Fatal(err)
	}
	if !contains(tr.Exclude, "etc/"+opqMarker) {
		t.Errorf("opaque marker not excluded: %v", tr.Exclude)
	}
	for _, want := range []string{"etc/hosts c 0 0 0 0 0", "etc/passwd c 0 0 0 0 0"} {
		if !contains(tr.Pseudo, want) {
			t.Errorf("missing %q; got %v", want, tr.Pseudo)
		}
	}
	if contains(tr.Pseudo, "etc/mine.conf c 0 0 0 0 0") {
		t.Error("whiteouted a file the overlay provides itself")
	}
	if !contains(tr.Opaque, "etc") {
		t.Errorf("opaque directory not reported: %v", tr.Opaque)
	}
}

// The kernel-overlayfs opaque convention has no marker file at all. A name-only
// scan finds nothing and the replaced directory silently reappears, so the xattr
// pass is what makes freeze correct on that driver.
func TestTranslateOpaqueXattrOnly(t *testing.T) {
	t.Parallel()
	img := newImage(t, `cd upper
mkdir etc
ea_set /upper/etc trusted.overlay.opaque y`)

	base := func(_ context.Context, dirs []string) (map[string][]string, error) {
		out := map[string][]string{}
		for _, d := range dirs {
			out[d] = []string{"hosts"}
		}
		return out, nil
	}
	entries, err := Walk(context.Background(), img)
	if err != nil {
		t.Fatal(err)
	}
	tr, err := Translate(context.Background(), img, entries, base)
	if err != nil {
		t.Fatal(err)
	}
	if !contains(tr.Opaque, "etc") {
		t.Fatalf("xattr-marked opaque directory not found: %+v", tr)
	}
	if !contains(tr.Pseudo, "etc/hosts c 0 0 0 0 0") {
		t.Errorf("no whiteout for the hidden base entry: %v", tr.Pseudo)
	}
}

// An overlay with no deletions needs no translation and must not invent one.
func TestTranslateNoDeletions(t *testing.T) {
	t.Parallel()
	img := newImage(t, "cd upper\nmkdir cnt_env")
	entries, err := Walk(context.Background(), img)
	if err != nil {
		t.Fatal(err)
	}
	tr, err := Translate(context.Background(), img, entries, noBase)
	if err != nil {
		t.Fatal(err)
	}
	if tr.Convention != ConventionNone || tr.Deletions() != 0 || len(tr.Exclude) != 0 {
		t.Errorf("invented a translation for a clean overlay: %+v", tr)
	}
}

// debugfs echoes each command through readline, which redraws one longer than the
// terminal width: the path arrives split by a carriage return and cursor
// sequences, or with its front cut off. Attribution must survive that, because
// the payload this exists for — a conda prefix — is full of paths long enough to
// trigger it, and reading the echo silently mislabelled whole subtrees: the
// whiteout markers under them stopped matching their exclude entries and were
// packed as ordinary files.
func TestParseListingIgnoresTheEchoedPath(t *testing.T) {
	dirs := []string{"cnt_env", "cnt_env/deep/nested/path"}
	// The second echo is what readline actually produces for a long line.
	out := "debugfs:  ls -p /upper/cnt_env\n" +
		"/15/040755/10295/10295/.//\n" +
		"/12/040755/10295/10295/..//\n" +
		"/16/100700/10295/10295/.wh..wh..opq/0/\n" +
		"debugfs:  ls -p /upper/cnt_env/deep/nested/pa\r\x1bM\x1b[C\x1b[C\x1b[Kth\n" +
		"/20/100644/10295/10295/tool/512/\n" +
		"debugfs:  quit\n"

	entries, err := parseListing(out, dirs)
	if err != nil {
		t.Fatalf("parseListing: %v", err)
	}
	got := paths(entries)
	want := []string{"cnt_env/.wh..wh..opq", "cnt_env/deep/nested/path/tool"}
	if len(got) != len(want) {
		t.Fatalf("entries = %v, want %v", got, want)
	}
	for i := range want {
		if got[i] != want[i] {
			t.Errorf("entry %d = %q, want %q", i, got[i], want[i])
		}
	}
}

// Output that cannot be attributed to every directory asked for is an error, not
// a partial answer: entries after the gap would be filed under the wrong parent.
func TestParseListingRefusesUnattributableOutput(t *testing.T) {
	_, err := parseListing("debugfs:  ls -p /upper/a\n/15/040755/0/0/f/1/\n",
		[]string{"a", "b"})
	if err == nil {
		t.Fatal("a short listing was accepted")
	}
}
