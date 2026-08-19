package sif

import (
	"os"
	"path/filepath"
	"testing"
)

// realSIFs are apptainer-produced base images to validate the parser against.
//
// A fixture this package builds itself proves only that buildSIF and the parser
// agree — if a layout constant is wrong, both are wrong the same way and the
// test still passes. That is not hypothetical: dataPartition was 4 until a real
// SIF showed the data type enumeration is 0x4000-based, so the on-disk value is
// 16388. Only a file this package did not write can catch that class of error.
func realSIFs(t *testing.T) []string {
	t.Helper()
	root, err := filepath.Abs(filepath.Join("..", "..", "..", "images"))
	if err != nil {
		t.Skip("cannot resolve the images directory")
	}
	matches, err := filepath.Glob(filepath.Join(root, "*.sif"))
	if err != nil || len(matches) == 0 {
		t.Skipf("no .sif images in %s to validate against", root)
	}
	return matches
}

// TestAgainstRealSIF is the parser's ground truth: the partition it reports has
// to start on an actual SquashFS superblock, and lie inside the file.
func TestAgainstRealSIF(t *testing.T) {
	for _, path := range realSIFs(t) {
		t.Run(filepath.Base(path), func(t *testing.T) {
			fi, err := os.Stat(path)
			if err != nil {
				t.Skipf("stat: %v", err)
			}

			part, err := PrimarySystemPartition(path)
			if err != nil {
				t.Fatalf("PrimarySystemPartition: %v", err)
			}
			if part.Offset <= 0 || part.Size <= 0 {
				t.Fatalf("partition = %+v", part)
			}
			if part.Offset+part.Size > fi.Size() {
				t.Fatalf("partition %+v runs past the end of a %d byte file", part, fi.Size())
			}

			f, err := os.Open(path)
			if err != nil {
				t.Fatalf("open: %v", err)
			}
			defer f.Close()

			// The one check that cannot be faked by agreeing with ourselves.
			magicAt := make([]byte, 4)
			if _, err := f.ReadAt(magicAt, part.Offset); err != nil {
				t.Fatalf("read at %d: %v", part.Offset, err)
			}
			if string(magicAt) != "hsqs" {
				t.Fatalf("offset %d holds %q, not the SquashFS magic", part.Offset, magicAt)
			}
			t.Logf("%s: partition at %d, %d bytes, superblock verified",
				filepath.Base(path), part.Offset, part.Size)
		})
	}
}
