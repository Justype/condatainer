package store

import (
	"errors"
	"fmt"
	"os"
	"path/filepath"
	"slices"
	"sort"
	"strings"
	"syscall"
	"time"

	"github.com/Justype/condatainer/internal/artifact/meta"
	"github.com/Justype/condatainer/internal/config"
	"github.com/Justype/condatainer/internal/image"
	"github.com/Justype/condatainer/internal/image/producer"
	"github.com/Justype/condatainer/internal/utils"
)

// ErrUnscoped reports a destructive GC that named no store to act on.
var ErrUnscoped = errors.New("gc --apply requires a store to act on")

// StagingGrace is how long an abandoned .part file is left before GC reports
// it. Fixed rather than configurable: it is crash recovery, not retention
// policy, and the only question it answers is whether a producer could still be
// writing. A day is far past any real build's publication step.
const StagingGrace = 24 * time.Hour

// Basis names the timestamp that decided an entry's age.
type Basis string

const (
	BasisAtime Basis = "atime"
	BasisMtime Basis = "mtime"
	BasisCtime Basis = "ctime"
)

// Collectable is one store entry old enough and free enough to delete.
type Collectable struct {
	Path     string      `json:"path"`
	Root     string      `json:"root"`
	Name     string      `json:"name"`
	Identity meta.KeyRef `json:"identity"`
	Size     int64       `json:"size"`
	// Age is measured from Basis, the newest of the three timestamps.
	Age   time.Duration `json:"age"`
	Basis Basis         `json:"basis"`
	// Staging marks an abandoned .part rather than a published entry.
	Staging bool `json:"staging,omitempty"`
	// Removed is set by an applied run.
	Removed bool `json:"removed,omitempty"`
}

// Retained is one entry GC will not touch, and why.
type Retained struct {
	Path   string `json:"path"`
	Reason string `json:"reason"`
}

// GCReport is one pass over the selected stores.
type GCReport struct {
	// Dirs are the image roots examined, in the order they were scanned.
	Dirs        []string      `json:"dirs"`
	Grace       time.Duration `json:"grace"`
	Applied     bool          `json:"applied"`
	Collectable []Collectable `json:"collectable,omitempty"`
	Retained    []Retained    `json:"retained,omitempty"`
	// Reclaimable is the total size of Collectable, whether or not it was freed.
	Reclaimable int64 `json:"reclaimable"`
}

// GCOptions selects the stores and the retention rule.
type GCOptions struct {
	// Dirs are image roots to examine. Empty means every writable configured
	// root, which is what an unscoped report walks.
	Dirs []string
	// Grace is the age below which nothing is collectable. Zero uses the
	// configured store_gc_grace.
	Grace time.Duration
	// Apply deletes rather than reporting. It refuses an empty Dirs.
	Apply bool
	// now is injected for tests.
	now time.Time
}

// GC reports what could be reclaimed from the selected stores, and with Apply
// set, reclaims it.
//
// Report and apply are one traversal rather than two phases, which is what
// closes the gap between deciding and deleting: the exclusive lock taken to
// judge an entry is still held when it is unlinked, so nothing can start
// reading it in between. A separate apply run repeats the whole judgement for
// the same reason — the report it was shown may be hours old.
//
// Every uncertainty retains. An entry that cannot be opened, locked, stat'd or
// validated stays, because a file that cannot be shown to be garbage is not
// garbage, and deleting on a read error turns corruption into data loss.
func GC(opts GCOptions) (*GCReport, error) {
	if opts.Apply && len(opts.Dirs) == 0 {
		return nil, fmt.Errorf("%w: name one with --dir or --layer", ErrUnscoped)
	}
	grace := opts.Grace
	if grace <= 0 {
		grace = config.Global.StoreGCGrace
	}
	now := opts.now
	if now.IsZero() {
		now = time.Now()
	}
	dirs := opts.Dirs
	if len(dirs) == 0 {
		dirs = writableRoots()
	}

	report := &GCReport{Dirs: dirs, Grace: grace, Applied: opts.Apply}
	verified := verifiedByPath(dirs)
	for _, root := range dirs {
		storeDir := filepath.Join(root, DirName)
		entries, err := os.ReadDir(storeDir)
		if err != nil {
			if !os.IsNotExist(err) {
				report.Retained = append(report.Retained, Retained{Path: storeDir, Reason: err.Error()})
			}
			continue
		}
		for _, entry := range entries {
			if entry.IsDir() {
				continue
			}
			path := filepath.Join(storeDir, entry.Name())
			collectable, retained := judge(path, root, verified, grace, now)
			if retained != nil {
				report.Retained = append(report.Retained, *retained)
				continue
			}
			if opts.Apply {
				if err := removeLocked(path); err != nil {
					report.Retained = append(report.Retained, Retained{Path: path, Reason: err.Error()})
					continue
				}
				collectable.Removed = true
			}
			report.Reclaimable += collectable.Size
			report.Collectable = append(report.Collectable, *collectable)
		}
	}
	sort.SliceStable(report.Collectable, func(i, j int) bool {
		return report.Collectable[i].Size > report.Collectable[j].Size
	})
	sort.Slice(report.Retained, func(i, j int) bool { return report.Retained[i].Path < report.Retained[j].Path })
	return report, nil
}

// judge decides one file's fate, returning exactly one of the two.
//
// The order is deliberate: age is cheapest and rejects most entries, then the
// lock, which is the only check that can say "a container is using this" or
// "someone pinned this". Validation comes last because it is the one that
// already ran during the scan.
func judge(path, root string, verified map[string]Candidate, grace time.Duration, now time.Time) (*Collectable, *Retained) {
	info, err := os.Lstat(path)
	if err != nil {
		return nil, &Retained{Path: path, Reason: err.Error()}
	}
	if !info.Mode().IsRegular() {
		return nil, &Retained{Path: path, Reason: "not a regular file"}
	}
	staging := strings.HasSuffix(path, producer.PreparedSuffix)
	window := grace
	if staging {
		window = StagingGrace
	}
	age, basis := ageOf(info, now)
	if age < window {
		return nil, &Retained{Path: path,
			Reason: fmt.Sprintf("%s is %s old, within the %s grace", basis, utils.FormatDuration(age), utils.FormatDuration(window))}
	}

	candidate, known := verified[path]
	if !staging && !known {
		// The scan already regenerated keys for everything readable here, so an
		// absence means the entry did not validate. Something unreadable under a
		// store filename is a repair job, not garbage.
		return nil, &Retained{Path: path, Reason: "did not validate; `condatainer store validate` explains why"}
	}

	lock, err := image.AcquireLock(path, true)
	if err != nil {
		return nil, &Retained{Path: path, Reason: err.Error()}
	}
	lock.Close() //nolint:errcheck

	return &Collectable{
		Path: path, Root: root, Name: candidate.Name, Identity: candidate.Identity,
		Size: info.Size(), Age: age, Basis: basis, Staging: staging,
	}, nil
}

// ageOf is how long ago the file last mattered, from the newest of its three
// timestamps.
//
// atime is the one that tracks use, and it is what HPC scratch purging is built
// on. Taking the newest is what makes every degradation safe: a spuriously
// refreshed atime — a backup, an indexer, a tree walk — makes an entry look
// newer and retains it, and a noatime mount freezes atime at creation so the
// other two answer, which is no worse than reading ctime alone. Neither
// deletes. Reporting the basis is what tells those two situations apart without
// having to go read mount options.
func ageOf(info os.FileInfo, now time.Time) (time.Duration, Basis) {
	newest, basis := info.ModTime(), BasisMtime
	if raw, ok := info.Sys().(*syscall.Stat_t); ok {
		if atime := time.Unix(raw.Atim.Sec, raw.Atim.Nsec); atime.After(newest) {
			newest, basis = atime, BasisAtime
		}
		if ctime := time.Unix(raw.Ctim.Sec, raw.Ctim.Nsec); ctime.After(newest) {
			newest, basis = ctime, BasisCtime
		}
	}
	age := now.Sub(newest)
	if age < 0 {
		// A clock skew across an NFS server is not evidence of anything; treat
		// the entry as brand new, which retains it.
		age = 0
	}
	return age, basis
}

// verifiedByPath indexes the verified store entries of dirs by absolute path,
// so one scan answers for every file GC walks.
func verifiedByPath(dirs []string) map[string]Candidate {
	report := Scan(ScanOptions{Dirs: dirs, Stored: true})
	out := make(map[string]Candidate, len(report.Candidates))
	for _, candidate := range report.Candidates {
		out[candidate.Path] = candidate
	}
	return out
}

// writableRoots is the configured image roots GC may delete from. A read-only
// root is not a candidate, and a root with no store is skipped rather than
// given one.
func writableRoots() []string {
	var out []string
	for _, root := range config.GetImageSearchPaths() {
		if !utils.CanWriteToDir(root) {
			continue
		}
		if !slices.Contains(out, root) {
			out = append(out, root)
		}
	}
	return out
}
