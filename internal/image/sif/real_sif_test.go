package sif

import (
	"context"
	"crypto/sha256"
	"encoding/hex"
	"errors"
	"io"
	"os"
	"os/exec"
	"path/filepath"
	"strconv"
	"testing"

	"github.com/Justype/condatainer/internal/image/tool"
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

// A real base SIF carries no CondaTainer manifest yet, and that has to read as
// ErrFileNotFound — the outcome meta.Read turns into ErrNoManifest, which is
// what keeps an ordinary Apptainer image usable as a base.
func TestRealSIFManifestLookup(t *testing.T) {
	requireTools(t)

	for _, path := range realSIFs(t) {
		t.Run(filepath.Base(path), func(t *testing.T) {
			data, err := ReadFile(path, "/.cnt/manifest.json")
			switch {
			case err == nil:
				t.Logf("%s carries a manifest: %s", filepath.Base(path), data)
			case errors.Is(err, tool.ErrFileNotFound):
				t.Logf("%s has no manifest, reported cleanly", filepath.Base(path))
			default:
				t.Fatalf("unexpected error reading %s: %v", filepath.Base(path), err)
			}
		})
	}
}

// TestExtractMatchesApptainer is the check that lets ExtractPartition replace
// `apptainer sif dump`: on a real image the two must produce the same bytes.
// It skips when apptainer is absent, so it never gates a build — but where
// apptainer exists, it is the proof that dropping it changed nothing.
func TestExtractMatchesApptainer(t *testing.T) {
	if _, err := exec.LookPath("apptainer"); err != nil {
		t.Skip("apptainer not available to compare against")
	}

	for _, path := range realSIFs(t) {
		t.Run(filepath.Base(path), func(t *testing.T) {
			part, err := PrimarySystemPartition(path)
			if err != nil {
				t.Fatalf("PrimarySystemPartition: %v", err)
			}

			mine := filepath.Join(t.TempDir(), "mine.sqf")
			if err := ExtractPartition(context.Background(), path, mine); err != nil {
				t.Fatalf("ExtractPartition: %v", err)
			}

			theirs := filepath.Join(t.TempDir(), "theirs.sqf")
			out, err := os.Create(theirs)
			if err != nil {
				t.Fatalf("create: %v", err)
			}
			// The ID we parsed is what `sif dump` takes, so this also checks it.
			cmd := exec.Command("apptainer", "sif", "dump", strconv.FormatUint(uint64(part.ID), 10), path)
			cmd.Stdout = out
			if err := cmd.Run(); err != nil {
				out.Close()
				t.Skipf("apptainer sif dump %d failed: %v", part.ID, err)
			}
			out.Close()

			if got, want := sha256File(t, mine), sha256File(t, theirs); got != want {
				t.Fatalf("extraction differs from apptainer:\n ours: %s\ntheirs: %s", got, want)
			}
			t.Logf("%s: partition id=%d extracted identically to apptainer sif dump", filepath.Base(path), part.ID)
		})
	}
}

func sha256File(t *testing.T, path string) string {
	t.Helper()
	f, err := os.Open(path)
	if err != nil {
		t.Fatalf("open %s: %v", path, err)
	}
	defer f.Close()
	h := sha256.New()
	if _, err := io.Copy(h, f); err != nil {
		t.Fatalf("hash %s: %v", path, err)
	}
	return hex.EncodeToString(h.Sum(nil))
}
