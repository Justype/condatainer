// Package sif reads the Singularity Image Format container layout: enough of a
// SIF's global header and descriptor table to find where its SquashFS payload
// starts, so a single file can be read out of it in place.
//
// This is deliberately not a dependency on github.com/apptainer/sif. The layout
// needed here is small and fixed, and this module's dependency set is kept
// minimal on purpose.
package sif

import (
	"bytes"
	"encoding/binary"
	"fmt"
	"io"
	"os"

	"github.com/Justype/condatainer/internal/image/tool"
)

// Global header layout, all integers little-endian. The header is the first 128
// bytes of the file: a 32-byte launch script, the magic, version and arch, a
// UUID, and then eight int64 fields.
const (
	headerSize = 128

	magicOffset = 32
	magicLen    = 10

	versionOffset = 42 // 3 bytes, "01\0"
	versionLen    = 3

	archOffset = 45 // 3 bytes, "02\0" for amd64 — see archNames
	archLen    = 3

	descrOffsetField = 96  // int64: byte offset of the descriptor table
	descrTotalField  = 88  // int64: number of descriptor slots
	descrSizeField   = 104 // int64: byte length of the descriptor table
)

// magic is what marks a file as a SIF. The trailing NUL is part of it.
var magic = []byte("SIF_MAGIC\x00")

// HasMagic reports whether a header begins a SIF, for a caller identifying a
// file whose extension does not name its format.
func HasMagic(header []byte) bool {
	if len(header) < magicOffset+magicLen {
		return false
	}
	return string(header[magicOffset:magicOffset+magicLen]) == string(magic)
}

// version is the SIF layout this parser knows. Every offset in this file is
// specific to it, so a different version is refused rather than misread: a
// wrong offset would silently produce a plausible number pointing at nothing.
var version = []byte("01\x00")

// Descriptor layout. Entries are fixed size, so the table is indexed rather
// than parsed sequentially.
const (
	descriptorSize = 585

	dDataType = 0   // int32
	dUsed     = 4   // bool, one byte
	dID       = 5   // uint32: the ID column `apptainer sif list` prints
	dOffset   = 17  // int64: where this object's bytes start in the file
	dSize     = 25  // int64: how many bytes it occupies
	dExtra    = 201 // 384 bytes of type-specific data

	// A partition descriptor's extra begins with these two int32 fields.
	pFsType   = dExtra + 0
	pPartType = dExtra + 4
)

// The descriptor and partition constants this reader matches against; the full
// enumerations are larger.
//
// Data types are 0x4000-based on disk, so a partition descriptor reads as 16388,
// not 4. A real base SIF holds 0x4001 (definition), 0x4006 (inspect metadata)
// and 0x4004 (the SquashFS partition).
const (
	dataPartition int32 = 0x4004 // descriptor holds a filesystem partition

	fsSquash    int32 = 1 // that partition is SquashFS
	partPrimSys int32 = 2 // and it is the primary system partition
)

// Partition is where a filesystem lives inside a SIF. Offset and Size bound a
// complete SquashFS archive, byte-identical to a standalone .sqf. ID is what
// `apptainer sif list` prints, carried only so a diagnostic can be compared
// against that output.
type Partition struct {
	ID     uint32
	Offset int64
	Size   int64
}

// PrimarySystemPartition returns the SquashFS partition a SIF boots from —
// what `apptainer sif list` shows as FS (Squashfs/*System).
//
// A file whose magic does not match is reported as ErrCorrupt rather than as a
// missing manifest, because "this is not a SIF" and "this SIF has no metadata"
// are different answers to different questions.
func PrimarySystemPartition(path string) (Partition, error) {
	f, err := os.Open(path)
	if err != nil {
		return Partition{}, fmt.Errorf("%w: %s: %w", tool.ErrUnreadable, path, err)
	}
	defer f.Close()

	fi, err := f.Stat()
	if err != nil {
		return Partition{}, fmt.Errorf("%w: %s: %w", tool.ErrUnreadable, path, err)
	}
	return parsePrimarySystemPartition(f, fi.Size(), path)
}

// parsePrimarySystemPartition is PrimarySystemPartition over bytes already
// opened. Everything the format decides is here, so a test can drive it from an
// in-memory image instead of packing a real archive to read four fields back
// out. path only names the source in errors.
func parsePrimarySystemPartition(r io.ReaderAt, size int64, path string) (Partition, error) {
	header := make([]byte, headerSize)
	if _, err := r.ReadAt(header, 0); err != nil {
		return Partition{}, fmt.Errorf("%w: %s: short header: %w", tool.ErrCorrupt, path, err)
	}
	if string(header[magicOffset:magicOffset+magicLen]) != string(magic) {
		return Partition{}, fmt.Errorf("%w: %s is not a SIF", tool.ErrCorrupt, path)
	}
	if got := header[versionOffset : versionOffset+versionLen]; string(got) != string(version) {
		return Partition{}, fmt.Errorf("%w: %s is SIF version %q, this build reads %q",
			tool.ErrCorrupt, path, trimNUL(got), trimNUL(version))
	}

	tableOffset := int64(binary.LittleEndian.Uint64(header[descrOffsetField:]))
	tableCount := int64(binary.LittleEndian.Uint64(header[descrTotalField:]))
	// The header is attacker-controllable in the sense that a corrupt file can
	// claim anything, so the table has to fit in the file before it is walked.
	if tableOffset < headerSize || tableCount <= 0 ||
		tableOffset+tableCount*descriptorSize > size {
		return Partition{}, fmt.Errorf("%w: %s: descriptor table does not fit (offset %d, count %d, file %d bytes)",
			tool.ErrCorrupt, path, tableOffset, tableCount, size)
	}

	entry := make([]byte, descriptorSize)
	for i := range tableCount {
		if _, err := r.ReadAt(entry, tableOffset+i*descriptorSize); err != nil {
			return Partition{}, fmt.Errorf("%w: %s: truncated descriptor table: %w", tool.ErrCorrupt, path, err)
		}
		if entry[dUsed] == 0 {
			continue
		}
		if int32(binary.LittleEndian.Uint32(entry[dDataType:])) != dataPartition {
			continue
		}
		if int32(binary.LittleEndian.Uint32(entry[pFsType:])) != fsSquash {
			continue
		}
		if int32(binary.LittleEndian.Uint32(entry[pPartType:])) != partPrimSys {
			continue
		}
		part := Partition{
			ID:     binary.LittleEndian.Uint32(entry[dID:]),
			Offset: int64(binary.LittleEndian.Uint64(entry[dOffset:])),
			Size:   int64(binary.LittleEndian.Uint64(entry[dSize:])),
		}
		if part.Offset < headerSize || part.Size <= 0 || part.Offset+part.Size > size {
			return Partition{}, fmt.Errorf("%w: %s: partition %+v does not fit in a %d byte file",
				tool.ErrCorrupt, path, part, size)
		}
		return part, nil
	}
	return Partition{}, fmt.Errorf("%w: %s has no primary SquashFS system partition", tool.ErrCorrupt, path)
}

// trimNUL renders a NUL-padded fixed-width header field for a message.
func trimNUL(b []byte) string {
	return string(bytes.TrimRight(b, "\x00"))
}

// archNames maps a SIF header's two-digit architecture code to a GOARCH
// string, values confirmed against github.com/apptainer/sif's pkg/sif/arch.go.
// Only the architectures a CondaTainer host runs on are listed; Arch returns
// any other code as-is.
var archNames = map[string]string{
	"02": "amd64",
	"04": "arm64",
}

// Arch returns the architecture recorded in a SIF's global header, as a
// GOARCH string ("amd64", "arm64") or the raw code if archNames does not
// name it.
func Arch(path string) (string, error) {
	f, err := os.Open(path)
	if err != nil {
		return "", fmt.Errorf("%w: %s: %w", tool.ErrUnreadable, path, err)
	}
	defer f.Close()

	header := make([]byte, headerSize)
	if _, err := f.ReadAt(header, 0); err != nil {
		return "", fmt.Errorf("%w: %s: short header: %w", tool.ErrCorrupt, path, err)
	}
	if string(header[magicOffset:magicOffset+magicLen]) != string(magic) {
		return "", fmt.Errorf("%w: %s is not a SIF", tool.ErrCorrupt, path)
	}
	if got := header[versionOffset : versionOffset+versionLen]; string(got) != string(version) {
		return "", fmt.Errorf("%w: %s is SIF version %q, this build reads %q",
			tool.ErrCorrupt, path, trimNUL(got), trimNUL(version))
	}

	code := trimNUL(header[archOffset : archOffset+archLen])
	if name, ok := archNames[code]; ok {
		return name, nil
	}
	return code, nil
}
