package freeze

import (
	"bytes"
	"context"
	"fmt"
	"os"
	"strconv"
	"strings"

	"github.com/Justype/condatainer/internal/artifact/key"
	"github.com/Justype/condatainer/internal/artifact/meta"
	execpkg "github.com/Justype/condatainer/internal/runtime/exec"
	"github.com/Justype/condatainer/internal/utils"
)

// treeScript walks the mounted artifact at mnt and prints one NUL-terminated
// field per value: type, mode, path and link target for every entry, a sha256
// for every regular file, and major:minor for every device — which find
// cannot print, so stat supplies it.
//
// NUL rather than newlines because a path may contain one, and the paths come
// from whatever the user installed. `xargs` batches the hashing into a handful of
// processes; one sha256sum per file is 190x slower and is the single mistake that
// makes this look unaffordable *(measured)*.
func treeScript(mnt string) string {
	return fmt.Sprintf(`cd %s || exit 1
find . -mindepth 1 -printf 'E\0%%y\0%%m\0%%p\0%%l\0'
find . -mindepth 1 -type f -print0 | xargs -0 -r sha256sum -z |
	while IFS= read -r -d '' line; do printf 'H\0%%s\0' "$line"; done
find . -mindepth 1 \( -type c -o -type b \) -print0 |
	xargs -0 -r stat --printf='D\0%%n\0%%t\0%%T\0'
`, shellQuote(mnt))
}

// TreeIdentity hashes the payload of a packed artifact into its identity.
//
// It reads the finished .sqf rather than the overlay it came from, for two
// reasons. The identity then describes what the archive actually contains, so a
// pack that dropped or renamed something is reflected rather than papered over;
// and reading a .sqf through squashfuse costs a fraction of reading the .img
// through fuse2fs — 2.4x the cost of a local read against 45x *(measured)*.
func TreeIdentity(ctx context.Context, artifact string) (meta.KeyRef, error) {
	squashfuse, err := findSquashfuse()
	if err != nil {
		return meta.KeyRef{}, err
	}

	scratch := utils.GetTmpDir()
	if err := os.MkdirAll(scratch, 0o755); err != nil {
		return meta.KeyRef{}, fmt.Errorf("stage identity mount: %w", err)
	}
	mnt, err := os.MkdirTemp(scratch, "cnt-ident-")
	if err != nil {
		return meta.KeyRef{}, fmt.Errorf("stage identity mount: %w", err)
	}
	defer os.RemoveAll(mnt)

	var stdout, stderr bytes.Buffer
	err = MountedRun(ctx, squashfuse, []string{artifact}, mnt, treeScript(mnt),
		execpkg.IO{Stdout: &stdout, Stderr: &stderr})
	if err != nil {
		return meta.KeyRef{}, fmt.Errorf("read %s to identify it: %w: %s", artifact, err, stderr.String())
	}

	records, err := parseTree(stdout.String())
	if err != nil {
		return meta.KeyRef{}, fmt.Errorf("read %s to identify it: %w", artifact, err)
	}
	return key.SnapshotIdentity(records)
}

// parseTree turns the walk's output into records, dropping the metadata
// directory: the identity covers the payload, so that the manifest carrying it
// can live inside the archive it describes.
//
// Two passes arrive interleaved and are joined by path — 'E' records carry every
// entry, 'H' records carry a hash for the regular ones. A hash for a path no 'E'
// record claimed means the two walks disagreed, which is a broken read rather
// than an entry to invent.
func parseTree(out string) ([]key.TreeRecord, error) {
	fields := strings.Split(out, "\x00")
	byPath := map[string]int{}
	var records []key.TreeRecord

	for i := 0; i < len(fields); {
		switch fields[i] {
		case "E":
			if i+5 > len(fields) {
				return nil, fmt.Errorf("truncated entry record at field %d", i)
			}
			typ, mode, p, link := fields[i+1], fields[i+2], fields[i+3], fields[i+4]
			i += 5
			clean := archivePath(p)
			if clean == "" || underMetaDir(clean) {
				continue
			}
			bits, err := strconv.ParseUint(mode, 8, 32)
			if err != nil {
				return nil, fmt.Errorf("%s: unreadable mode %q", clean, mode)
			}
			r := key.TreeRecord{Type: typ[0], Mode: uint32(bits), Path: clean}
			if r.Type == 'l' {
				r.ID = link
			}
			byPath[clean] = len(records)
			records = append(records, r)
		case "H":
			if i+2 > len(fields) {
				return nil, fmt.Errorf("truncated hash record at field %d", i)
			}
			line := fields[i+1]
			i += 2
			sum, p, ok := strings.Cut(line, "  ")
			if !ok {
				return nil, fmt.Errorf("unreadable sha256sum line %q", line)
			}
			clean := archivePath(p)
			if underMetaDir(clean) {
				continue
			}
			at, ok := byPath[clean]
			if !ok {
				return nil, fmt.Errorf("hashed %s, which the entry walk did not list", clean)
			}
			records[at].ID = sum
		case "D":
			if i+4 > len(fields) {
				return nil, fmt.Errorf("truncated device record at field %d", i)
			}
			p, major, minor := fields[i+1], fields[i+2], fields[i+3]
			i += 4
			clean := archivePath(p)
			if underMetaDir(clean) {
				continue
			}
			at, ok := byPath[clean]
			if !ok {
				return nil, fmt.Errorf("stat'd %s, which the entry walk did not list", clean)
			}
			// stat prints them in hex; a whiteout is 0:0 and every deletion in the
			// archive is one, so this is the field that keeps them distinct from a
			// device node someone really created.
			maj, err := strconv.ParseUint(major, 16, 32)
			if err != nil {
				return nil, fmt.Errorf("%s: unreadable device major %q", clean, major)
			}
			min, err := strconv.ParseUint(minor, 16, 32)
			if err != nil {
				return nil, fmt.Errorf("%s: unreadable device minor %q", clean, minor)
			}
			records[at].ID = fmt.Sprintf("%d:%d", maj, min)
		default:
			// Trailing empty field from the final NUL, and nothing else: a real
			// unknown tag means the streams have desynchronised.
			if strings.TrimSpace(fields[i]) != "" {
				return nil, fmt.Errorf("unexpected record tag %q at field %d", fields[i], i)
			}
			i++
		}
	}

	for _, r := range records {
		switch {
		case r.Type == 'f' && r.ID == "":
			return nil, fmt.Errorf("%s is a regular file with no hash", r.Path)
		case (r.Type == 'c' || r.Type == 'b') && r.ID == "":
			return nil, fmt.Errorf("%s is a device with no major:minor", r.Path)
		}
	}
	return records, nil
}

// archivePath turns find's "./x/y" into the "/x/y" a mount exposes, which is
// what the record has to name: the identity is about what the artifact provides,
// not about where it was read from.
func archivePath(p string) string {
	p = strings.TrimPrefix(p, ".")
	if p == "" || p == "/" {
		return ""
	}
	if !strings.HasPrefix(p, "/") {
		p = "/" + p
	}
	return p
}

// underMetaDir reports a path inside /.cnt, which the identity excludes.
func underMetaDir(p string) bool {
	dir := "/" + meta.DirName
	return p == dir || strings.HasPrefix(p, dir+"/")
}
