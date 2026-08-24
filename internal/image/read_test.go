package image

import (
	"errors"
	"os"
	"os/exec"
	"path/filepath"
	"testing"

	"github.com/Justype/condatainer/internal/image/tool"
)

// packSqf builds a real SquashFS archive holding one directory, and returns it
// under the given filename. The name is the point of these tests, so it is the
// caller's to choose.
func packSqf(t *testing.T, filename string) string {
	t.Helper()
	if _, err := exec.LookPath("mksquashfs"); err != nil {
		t.Skip("mksquashfs not available")
	}
	if _, err := exec.LookPath("unsquashfs"); err != nil {
		t.Skip("unsquashfs not available")
	}

	root := t.TempDir()
	inner := filepath.Join(root, ".cnt")
	if err := os.MkdirAll(inner, 0o775); err != nil {
		t.Fatal(err)
	}
	if err := os.WriteFile(filepath.Join(inner, "manifest.json"), []byte(`{"schema_version":1}`), 0o664); err != nil {
		t.Fatal(err)
	}

	out := filepath.Join(t.TempDir(), filename)
	cmd := exec.Command("mksquashfs", root, out, "-no-progress", "-noappend", "-quiet", "-no-xattrs")
	if output, err := cmd.CombinedOutput(); err != nil {
		t.Fatalf("mksquashfs: %v\n%s", err, output)
	}
	return out
}

// The extension is the fast path, never the authority. A producer's in-progress
// output is `<name>.sqf.<tag>.part`, so store.Commit has to verify a file whose
// name says nothing — every commit fails if dispatch stops at the extension.
func TestExtractDirReadsAnImageWhoseExtensionNamesNoFormat(t *testing.T) {
	image := packSqf(t, "artifact.sqf.local-node-1234.part")

	dest := t.TempDir()
	if err := ExtractDir(image, "/.cnt", dest); err != nil {
		t.Fatalf("ExtractDir on a .part: %v", err)
	}
	if _, err := os.Stat(filepath.Join(dest, ".cnt", "manifest.json")); err != nil {
		t.Fatalf("extracted tree is missing the file: %v", err)
	}
}

// The ordinary case must keep working, and by the extension: sniffing every
// image would open and read each one before dispatching.
func TestExtractDirStillReadsAPlainSqf(t *testing.T) {
	image := packSqf(t, "artifact.sqf")

	dest := t.TempDir()
	if err := ExtractDir(image, "/.cnt", dest); err != nil {
		t.Fatalf("ExtractDir: %v", err)
	}
	if _, err := os.Stat(filepath.Join(dest, ".cnt", "manifest.json")); err != nil {
		t.Fatalf("extracted tree is missing the file: %v", err)
	}
}

// Sniffing must not turn "unreadable" into a confusing tool failure: a file
// that is not an image at all is corrupt, and says so.
func TestExtractDirRejectsBytesThatAreNoImage(t *testing.T) {
	path := filepath.Join(t.TempDir(), "notanimage.part")
	if err := os.WriteFile(path, []byte("this is not an image, it is a note"), 0o664); err != nil {
		t.Fatal(err)
	}
	err := ExtractDir(path, "/.cnt", t.TempDir())
	if !errors.Is(err, tool.ErrCorrupt) {
		t.Fatalf("want ErrCorrupt, got %v", err)
	}
}

// A file shorter than the header cannot be sniffed. Reading it must fail as
// corrupt rather than panic on a short slice.
func TestExtractDirRejectsAFileTooShortToSniff(t *testing.T) {
	path := filepath.Join(t.TempDir(), "truncated.part")
	if err := os.WriteFile(path, []byte("hsq"), 0o664); err != nil {
		t.Fatal(err)
	}
	err := ExtractDir(path, "/.cnt", t.TempDir())
	if !errors.Is(err, tool.ErrCorrupt) {
		t.Fatalf("want ErrCorrupt, got %v", err)
	}
}
