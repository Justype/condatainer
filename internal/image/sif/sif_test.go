package sif

import (
	"encoding/binary"
	"errors"
	"os"
	"os/exec"
	"path/filepath"
	"testing"

	"github.com/Justype/condatainer/internal/image/tool"
)

// partitionSpec is one descriptor a fixture should contain.
type partitionSpec struct {
	used     bool
	dataType int32
	fsType   int32
	partType int32
}

// primarySquashfs is the descriptor a real CondaTainer base image has.
var primarySquashfs = partitionSpec{used: true, dataType: dataPartition, fsType: fsSquash, partType: partPrimSys}

// buildSIF writes a SIF whose descriptor table holds specs and whose payload is
// a real SquashFS archive of dir, placed after the table. The layout mirrors
// what apptainer produces closely enough to exercise the offset arithmetic:
// the partition is read in place, so if the offset is wrong, unsquashfs fails.
func buildSIF(t *testing.T, dir string, specs []partitionSpec) string {
	t.Helper()

	sqf := filepath.Join(t.TempDir(), "payload.sqf")
	cmd := exec.Command("mksquashfs", dir, sqf, "-no-progress", "-noappend", "-quiet")
	if out, err := cmd.CombinedOutput(); err != nil {
		t.Fatalf("mksquashfs: %v\n%s", err, out)
	}
	payload, err := os.ReadFile(sqf)
	if err != nil {
		t.Fatalf("read payload: %v", err)
	}

	tableOffset := int64(headerSize)
	tableBytes := int64(len(specs)) * descriptorSize
	// Round the data start up to 4 KiB, as a real SIF does.
	dataOffset := (tableOffset + tableBytes + 4095) / 4096 * 4096

	out := make([]byte, dataOffset+int64(len(payload)))
	copy(out[magicOffset:], magic)
	copy(out[42:], []byte("01\x00"))
	binary.LittleEndian.PutUint64(out[descrTotalField:], uint64(len(specs)))
	binary.LittleEndian.PutUint64(out[descrOffsetField:], uint64(tableOffset))
	binary.LittleEndian.PutUint64(out[descrSizeField:], uint64(tableBytes))
	binary.LittleEndian.PutUint64(out[112:], uint64(dataOffset))
	binary.LittleEndian.PutUint64(out[120:], uint64(len(payload)))

	for i, spec := range specs {
		d := out[tableOffset+int64(i)*descriptorSize:]
		binary.LittleEndian.PutUint32(d[dDataType:], uint32(spec.dataType))
		if spec.used {
			d[dUsed] = 1
		}
		binary.LittleEndian.PutUint64(d[dOffset:], uint64(dataOffset))
		binary.LittleEndian.PutUint64(d[dSize:], uint64(len(payload)))
		binary.LittleEndian.PutUint32(d[pFsType:], uint32(spec.fsType))
		binary.LittleEndian.PutUint32(d[pPartType:], uint32(spec.partType))
	}
	copy(out[dataOffset:], payload)

	path := filepath.Join(t.TempDir(), "image.sif")
	if err := os.WriteFile(path, out, 0o644); err != nil {
		t.Fatalf("write sif: %v", err)
	}
	return path
}

// payloadDir returns a directory holding one file at innerPath.
func payloadDir(t *testing.T, innerPath, content string) string {
	t.Helper()
	root := t.TempDir()
	full := filepath.Join(root, innerPath)
	if err := os.MkdirAll(filepath.Dir(full), 0o755); err != nil {
		t.Fatalf("mkdir: %v", err)
	}
	if err := os.WriteFile(full, []byte(content), 0o644); err != nil {
		t.Fatalf("write: %v", err)
	}
	return root
}

func requireTools(t *testing.T) {
	t.Helper()
	for _, bin := range []string{"mksquashfs", "unsquashfs"} {
		if _, err := exec.LookPath(bin); err != nil {
			t.Skipf("%s not available", bin)
		}
	}
}

func TestPrimarySystemPartition(t *testing.T) {
	requireTools(t)

	path := buildSIF(t, payloadDir(t, ".cnt/manifest.json", "{}"), []partitionSpec{primarySquashfs})
	part, err := PrimarySystemPartition(path)
	if err != nil {
		t.Fatalf("PrimarySystemPartition: %v", err)
	}
	if part.Offset <= 0 {
		t.Errorf("offset = %d, want a positive byte offset", part.Offset)
	}
	if part.Size <= 0 {
		t.Errorf("size = %d", part.Size)
	}

	// The offset has to land on the SquashFS superblock, which starts "hsqs".
	f, err := os.Open(path)
	if err != nil {
		t.Fatalf("open: %v", err)
	}
	defer f.Close()
	magicAt := make([]byte, 4)
	if _, err := f.ReadAt(magicAt, part.Offset); err != nil {
		t.Fatalf("read at partition offset: %v", err)
	}
	if string(magicAt) != "hsqs" {
		t.Errorf("bytes at offset = %q, want the SquashFS magic", magicAt)
	}
}

// The reader must pick the primary system partition specifically, not the first
// partition descriptor it happens to see.
func TestPrimarySystemPartitionSkipsOthers(t *testing.T) {
	requireTools(t)

	specs := []partitionSpec{
		{used: false, dataType: dataPartition, fsType: fsSquash, partType: partPrimSys}, // unused
		{used: true, dataType: 0x4001, fsType: fsSquash, partType: partPrimSys},         // a deffile, not a partition
		{used: true, dataType: dataPartition, fsType: 2, partType: partPrimSys},         // ext3, not squashfs
		{used: true, dataType: dataPartition, fsType: fsSquash, partType: 3},            // a data partition
		primarySquashfs,
	}
	path := buildSIF(t, payloadDir(t, ".cnt/manifest.json", "{}"), specs)

	if _, err := PrimarySystemPartition(path); err != nil {
		t.Fatalf("PrimarySystemPartition: %v", err)
	}
}

func TestPrimarySystemPartitionMissing(t *testing.T) {
	requireTools(t)

	specs := []partitionSpec{{used: true, dataType: dataPartition, fsType: 2, partType: partPrimSys}}
	path := buildSIF(t, payloadDir(t, "file", "x"), specs)

	_, err := PrimarySystemPartition(path)
	if !errors.Is(err, tool.ErrCorrupt) {
		t.Fatalf("err = %v, want ErrCorrupt", err)
	}
}

// "This is not a SIF" is a different answer from "this SIF has no metadata",
// and must not be reported as the latter.
func TestPrimarySystemPartitionNotASIF(t *testing.T) {
	path := filepath.Join(t.TempDir(), "plain.sif")
	if err := os.WriteFile(path, make([]byte, headerSize*2), 0o644); err != nil {
		t.Fatalf("write: %v", err)
	}
	if _, err := PrimarySystemPartition(path); !errors.Is(err, tool.ErrCorrupt) {
		t.Fatalf("err = %v, want ErrCorrupt", err)
	}

	// A file too short to hold a header is equally not a SIF.
	short := filepath.Join(t.TempDir(), "short.sif")
	if err := os.WriteFile(short, []byte("tiny"), 0o644); err != nil {
		t.Fatalf("write: %v", err)
	}
	if _, err := PrimarySystemPartition(short); !errors.Is(err, tool.ErrCorrupt) {
		t.Fatalf("err = %v, want ErrCorrupt", err)
	}

	if _, err := PrimarySystemPartition(filepath.Join(t.TempDir(), "absent.sif")); !errors.Is(err, tool.ErrUnreadable) {
		t.Error("a missing file should be ErrUnreadable")
	}
}

func TestReadFile(t *testing.T) {
	requireTools(t)

	const manifest = `{"schema_version":1,"name":"ubuntu24/base"}`
	path := buildSIF(t, payloadDir(t, ".cnt/manifest.json", manifest), []partitionSpec{primarySquashfs})

	data, err := ReadFile(path, "/.cnt/manifest.json")
	if err != nil {
		t.Fatalf("ReadFile: %v", err)
	}
	if string(data) != manifest {
		t.Errorf("read %q, want %q", data, manifest)
	}
}

// A SIF without the file is ErrFileNotFound, which is what lets a caller report
// a plain Apptainer image as "no manifest" rather than as a fault.
func TestReadFileMissing(t *testing.T) {
	requireTools(t)

	path := buildSIF(t, payloadDir(t, "usr/bin/sh", "x"), []partitionSpec{primarySquashfs})

	_, err := ReadFile(path, "/.cnt/manifest.json")
	if !errors.Is(err, tool.ErrFileNotFound) {
		t.Fatalf("err = %v, want ErrFileNotFound", err)
	}
}
