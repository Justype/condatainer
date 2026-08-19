package sif

import (
	"bytes"
	"encoding/binary"
	"errors"
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

// payloadSlot is how much room each descriptor's payload is given.
const payloadSlot = 512

// partitionOffset is where the nth descriptor's payload lands in sifBytes.
func partitionOffset(specs []partitionSpec, n int) int64 {
	tableBytes := int64(len(specs)) * descriptorSize
	dataOffset := (int64(headerSize) + tableBytes + 4095) / 4096 * 4096
	return dataOffset + int64(n)*payloadSlot
}

// sifBytes returns a SIF image in memory whose descriptor table holds specs.
//
// The payload is filler that only has to start with the SquashFS magic: this
// package locates a partition, it never decompresses one, so packing a real
// archive would cost an mksquashfs run per test to exercise the same arithmetic.
//
// Each descriptor gets its own offset, one payload slot apart, so a test can
// tell which one the parser chose. Giving them all the same offset would let a
// reader that picked the wrong descriptor still return a valid-looking result.
func sifBytes(t *testing.T, specs []partitionSpec) []byte {
	t.Helper()

	payload := append([]byte("hsqs"), make([]byte, payloadSlot-4)...)
	tableOffset := int64(headerSize)
	tableBytes := int64(len(specs)) * descriptorSize
	// Round the data start up to 4 KiB, as a real SIF does.
	dataOffset := (tableOffset + tableBytes + 4095) / 4096 * 4096

	out := make([]byte, dataOffset+int64(len(specs))*payloadSlot)
	copy(out[magicOffset:], magic)
	copy(out[versionOffset:], version)
	binary.LittleEndian.PutUint64(out[descrTotalField:], uint64(len(specs)))
	binary.LittleEndian.PutUint64(out[descrOffsetField:], uint64(tableOffset))
	binary.LittleEndian.PutUint64(out[descrSizeField:], uint64(tableBytes))

	for i, spec := range specs {
		d := out[tableOffset+int64(i)*descriptorSize:]
		at := partitionOffset(specs, i)
		binary.LittleEndian.PutUint32(d[dDataType:], uint32(spec.dataType))
		if spec.used {
			d[dUsed] = 1
		}
		binary.LittleEndian.PutUint32(d[dID:], uint32(i+1))
		binary.LittleEndian.PutUint64(d[dOffset:], uint64(at))
		binary.LittleEndian.PutUint64(d[dSize:], uint64(len(payload)))
		binary.LittleEndian.PutUint32(d[pFsType:], uint32(spec.fsType))
		binary.LittleEndian.PutUint32(d[pPartType:], uint32(spec.partType))
		copy(out[at:], payload)
	}
	return out
}

// parseSIF runs the parser over an in-memory image.
func parseSIF(t *testing.T, image []byte) (Partition, error) {
	t.Helper()
	return parsePrimarySystemPartition(bytes.NewReader(image), int64(len(image)), "image.sif")
}

func TestPrimarySystemPartition(t *testing.T) {
	specs := []partitionSpec{primarySquashfs}
	image := sifBytes(t, specs)

	part, err := parseSIF(t, image)
	if err != nil {
		t.Fatalf("parse: %v", err)
	}
	if want := partitionOffset(specs, 0); part.Offset != want {
		t.Errorf("offset = %d, want %d", part.Offset, want)
	}
	if part.Size <= 0 {
		t.Errorf("size = %d", part.Size)
	}
	// The offset has to land on the SquashFS superblock, which starts "hsqs".
	if got := string(image[part.Offset : part.Offset+4]); got != "hsqs" {
		t.Errorf("bytes at offset = %q, want the SquashFS magic", got)
	}
}

// The reader must pick the primary system partition specifically, not the first
// partition descriptor it happens to see. Each descriptor points at its own
// payload, so the returned offset says which one was chosen.
func TestPrimarySystemPartitionSkipsOthers(t *testing.T) {
	specs := []partitionSpec{
		{used: false, dataType: dataPartition, fsType: fsSquash, partType: partPrimSys}, // unused
		{used: true, dataType: 0x4001, fsType: fsSquash, partType: partPrimSys},         // a deffile, not a partition
		{used: true, dataType: dataPartition, fsType: 2, partType: partPrimSys},         // ext3, not squashfs
		{used: true, dataType: dataPartition, fsType: fsSquash, partType: 3},            // a data partition
		primarySquashfs,
	}

	part, err := parseSIF(t, sifBytes(t, specs))
	if err != nil {
		t.Fatalf("parse: %v", err)
	}
	last := len(specs) - 1
	if want := partitionOffset(specs, last); part.Offset != want {
		t.Errorf("offset = %d, want the primary partition at %d", part.Offset, want)
	}
	if want := uint32(last + 1); part.ID != want {
		t.Errorf("ID = %d, want the primary partition's %d", part.ID, want)
	}
}

func TestPrimarySystemPartitionMissing(t *testing.T) {
	specs := []partitionSpec{{used: true, dataType: dataPartition, fsType: 2, partType: partPrimSys}}

	_, err := parseSIF(t, sifBytes(t, specs))
	if !errors.Is(err, tool.ErrCorrupt) {
		t.Fatalf("err = %v, want ErrCorrupt", err)
	}
}

// A descriptor may claim anything, so a partition that does not fit inside the
// file has to be refused rather than handed out as an offset to read at.
func TestPrimarySystemPartitionOutOfBounds(t *testing.T) {
	specs := []partitionSpec{primarySquashfs}
	image := sifBytes(t, specs)

	d := image[headerSize:]
	binary.LittleEndian.PutUint64(d[dSize:], uint64(len(image)))

	if _, err := parseSIF(t, image); !errors.Is(err, tool.ErrCorrupt) {
		t.Fatalf("err = %v, want ErrCorrupt", err)
	}
}

// "This is not a SIF" is a different answer from "this SIF has no metadata",
// and must not be reported as the latter.
func TestPrimarySystemPartitionNotASIF(t *testing.T) {
	if _, err := parseSIF(t, make([]byte, headerSize*2)); !errors.Is(err, tool.ErrCorrupt) {
		t.Fatalf("err = %v, want ErrCorrupt", err)
	}

	// A file too short to hold a header is equally not a SIF.
	if _, err := parseSIF(t, []byte("tiny")); !errors.Is(err, tool.ErrCorrupt) {
		t.Fatalf("err = %v, want ErrCorrupt", err)
	}

	// A SIF this parser does not know the layout of is refused, not misread.
	wrongVersion := sifBytes(t, []partitionSpec{primarySquashfs})
	copy(wrongVersion[versionOffset:], []byte("99\x00"))
	if _, err := parseSIF(t, wrongVersion); !errors.Is(err, tool.ErrCorrupt) {
		t.Fatalf("err = %v, want ErrCorrupt", err)
	}

	// A file that is not there at all is the wrapper's answer, not the parser's.
	if _, err := PrimarySystemPartition("/nonexistent/absent.sif"); !errors.Is(err, tool.ErrUnreadable) {
		t.Error("a missing file should be ErrUnreadable")
	}
}
