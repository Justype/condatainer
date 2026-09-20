package key

import (
	"context"
	"crypto/sha256"
	"encoding/hex"
	"errors"
	"fmt"
	"io"
	"io/fs"
	"os"
	"path"
	"path/filepath"
	"sort"
	"strings"
	"sync"
	"syscall"

	"github.com/Justype/condatainer/internal/artifact/meta"
)

// ErrUnkeyed reports an artifact whose keys cannot be regenerated from vendored
// sources. Only a snapshot is, and callers that can recompute from the payload
// branch on this rather than treating it as a missing key.
var ErrUnkeyed = errors.New("artifact carries no regenerable keys")

// IsSnapshot reports whether a manifest describes a frozen writable overlay.
func IsSnapshot(m meta.Manifest) bool { return m.BuildType == meta.BuildTypeSnapshot }

// PayloadTreeV1 names how an artifact's payload is keyed: a sha256 over one
// record per archive entry, sorted by path. It is a frozen environment's identity
// and every built artifact's payload key.
//
// It is the artifact's payload that is hashed, never the packed file. A repack
// changes the file — mksquashfs stamps its own creation time into the superblock
// — while a payload that has not changed must keep its key.
const PayloadTreeV1 Scheme = "payload-tree-v1"

// TreeRecord is one archive entry's contribution to a payload key's preimage.
//
// Only what changes behaviour is carried. Mode is here because the executable
// bit decides whether a binary runs; mtime, uid and gid are not, because
// touching a file does not make a different environment and -all-root flattens
// ownership to a constant. Size is absent as the content hash subsumes it.
type TreeRecord struct {
	// Type is the entry's kind, spelled as find does: f d l c b p s.
	Type byte
	// Mode is the permission bits.
	Mode uint32
	// Path is the entry's full path inside the archive.
	Path string
	// ID identifies the entry's content: the sha256 of the bytes for a regular
	// file, the target for a symlink, "major:minor" for a device, empty for
	// anything with no content of its own.
	ID string
}

// PayloadKey hashes a payload's records into its key.
//
// Records are sorted here rather than trusted in the order they arrive, so the
// ordering is a property of the scheme and not of whatever walked the archive.
// Byte order, because a locale-aware comparison would give one tree different
// identities on two machines.
func PayloadKey(records []TreeRecord) (meta.KeyRef, error) {
	if len(records) == 0 {
		return meta.KeyRef{}, errors.New("payload key: no entries to hash")
	}
	sum := sha256.Sum256([]byte(PayloadPreimage(records)))
	return meta.KeyRef{Scheme: string(PayloadTreeV1), SHA256: hex.EncodeToString(sum[:])}, nil
}

// preimage renders one record. The NUL before the content id is what keeps a
// path holding a space or a newline from shifting the framing — paths come from
// whatever the user installed, so they are untrusted input.
func (r TreeRecord) preimage() string {
	return fmt.Sprintf("%c %04o %s\x00%s\n", r.Type, r.Mode&0o7777, r.Path, r.ID)
}

// PayloadPreimage renders the records PayloadKey hashes — the one place
// the ordering and the framing are decided, so a diff of two of these is exactly
// what the two digests were taken over. Exported because that is how a mismatch
// is read: by comparing these, not by staring at two digests.
func PayloadPreimage(records []TreeRecord) string {
	sorted := make([]TreeRecord, len(records))
	copy(sorted, records)
	sort.Slice(sorted, func(i, j int) bool { return sorted[i].Path < sorted[j].Path })
	var b strings.Builder
	for _, r := range sorted {
		b.WriteString(r.preimage())
	}
	return b.String()
}

// ArtifactDigest is the sha256 of the packed file, in OCI form.
//
// Not an identity: it changes on every repack, so nothing may pin it. It exists
// because a registry descriptor speaks this and nothing else.
func ArtifactDigest(path string) (string, error) {
	sum, err := sumFile(path)
	if err != nil {
		return "", err
	}
	return "sha256:" + sum, nil
}

// sumFile streams a sha256 and returns it as bare hex, which is what meta.KeyRef
// holds and what validSHA256 accepts.
func sumFile(path string) (string, error) {
	f, err := os.Open(path)
	if err != nil {
		return "", fmt.Errorf("digest %s: %w", path, err)
	}
	defer f.Close()

	sum := sha256.New()
	if _, err := io.Copy(sum, f); err != nil {
		return "", fmt.Errorf("digest %s: %w", path, err)
	}
	return hex.EncodeToString(sum.Sum(nil)), nil
}

// TreeOf records the payload under dir as the archive will hold it: the entries
// are named base plus their path below dir, and dir itself is the entry base.
// An empty base is an archive whose root is dir, so dir is not an entry.
//
// Contents are hashed by up to workers goroutines, each holding one buffer. A
// file is never split, so one large file costs one core.
func TreeOf(ctx context.Context, dir, base string, workers int) ([]TreeRecord, error) {
	if workers < 1 {
		workers = 1
	}
	var records []TreeRecord
	var files []int
	err := filepath.WalkDir(dir, func(p string, d fs.DirEntry, err error) error {
		if err != nil {
			return err
		}
		if err := ctx.Err(); err != nil {
			return err
		}
		rel, err := filepath.Rel(dir, p)
		if err != nil {
			return err
		}
		name := base
		if rel != "." {
			name = base + "/" + filepath.ToSlash(rel)
		} else if base == "" {
			return nil
		}
		// The metadata directory carries the manifest that holds this key.
		if name == "/"+meta.DirName {
			return fs.SkipDir
		}
		info, err := os.Lstat(p)
		if err != nil {
			return err
		}
		r, err := recordOf(p, path.Clean("/"+name), info)
		if err != nil {
			return err
		}
		if r.Type == 'f' {
			files = append(files, len(records))
		}
		records = append(records, r)
		return nil
	})
	if err != nil {
		return nil, fmt.Errorf("read payload %s: %w", dir, err)
	}
	if err := hashFiles(ctx, dir, base, records, files, workers); err != nil {
		return nil, err
	}
	return records, nil
}

// recordOf describes one entry apart from a regular file's content.
func recordOf(p, name string, info fs.FileInfo) (TreeRecord, error) {
	r := TreeRecord{Path: name, Mode: uint32(info.Mode().Perm())}
	if info.Mode()&fs.ModeSetuid != 0 {
		r.Mode |= 0o4000
	}
	if info.Mode()&fs.ModeSetgid != 0 {
		r.Mode |= 0o2000
	}
	if info.Mode()&fs.ModeSticky != 0 {
		r.Mode |= 0o1000
	}
	switch mode := info.Mode(); {
	case mode.IsRegular():
		r.Type = 'f'
	case mode.IsDir():
		r.Type = 'd'
	case mode&fs.ModeSymlink != 0:
		target, err := os.Readlink(p)
		if err != nil {
			return r, err
		}
		r.Type, r.ID = 'l', target
	case mode&fs.ModeNamedPipe != 0:
		r.Type = 'p'
	case mode&fs.ModeSocket != 0:
		r.Type = 's'
	case mode&fs.ModeDevice != 0:
		st, ok := info.Sys().(*syscall.Stat_t)
		if !ok {
			return r, fmt.Errorf("%s: no device numbers", p)
		}
		r.Type = 'b'
		if mode&fs.ModeCharDevice != 0 {
			r.Type = 'c'
		}
		r.ID = fmt.Sprintf("%d:%d", unixMajor(uint64(st.Rdev)), unixMinor(uint64(st.Rdev)))
	default:
		return r, fmt.Errorf("%s: unsupported file type %s", p, mode.Type())
	}
	return r, nil
}

// unixMajor and unixMinor split a Linux device number.
func unixMajor(dev uint64) uint32 {
	return uint32((dev>>8)&0xfff) | uint32((dev>>32)&^0xfff)
}

func unixMinor(dev uint64) uint32 {
	return uint32(dev&0xff) | uint32((dev>>12)&^0xff)
}

// hashFiles fills in the content id of the regular files at records[files[i]].
func hashFiles(ctx context.Context, dir, base string, records []TreeRecord, files []int, workers int) error {
	if workers > len(files) {
		workers = len(files)
	}
	jobs := make(chan int)
	var wg sync.WaitGroup
	var once sync.Once
	var firstErr error
	fail := func(err error) { once.Do(func() { firstErr = err }) }

	for w := 0; w < workers; w++ {
		wg.Add(1)
		go func() {
			defer wg.Done()
			buf := make([]byte, 1<<20)
			for i := range jobs {
				rel := strings.TrimPrefix(strings.TrimPrefix(records[i].Path, base), "/")
				id, err := hashFile(filepath.Join(dir, filepath.FromSlash(rel)), buf)
				if err != nil {
					fail(err)
					continue
				}
				records[i].ID = id
			}
		}()
	}
	for _, i := range files {
		if ctx.Err() != nil {
			break
		}
		jobs <- i
	}
	close(jobs)
	wg.Wait()
	if firstErr != nil {
		return firstErr
	}
	return ctx.Err()
}

func hashFile(p string, buf []byte) (string, error) {
	f, err := os.Open(p)
	if err != nil {
		return "", err
	}
	defer f.Close()
	sum := sha256.New()
	if _, err := io.CopyBuffer(sum, f, buf); err != nil {
		return "", fmt.Errorf("hash %s: %w", p, err)
	}
	return hex.EncodeToString(sum.Sum(nil)), nil
}
