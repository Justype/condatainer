package store

import (
	"errors"
	"os"
	"path/filepath"
	"strings"
	"testing"
	"time"

	"github.com/Justype/condatainer/internal/image"
	"github.com/Justype/condatainer/internal/image/producer"
)

// storeEntry writes a file into root's store/.
//
// It does not backdate. os.Chtimes moves atime and mtime but sets ctime to now,
// and ctime cannot be written from userspace at all, so a file on a real
// filesystem can never be made to look old to a rule that takes the newest of
// the three. Tests move the clock forward instead, which is what GCOptions.now
// is for.
func storeEntry(t *testing.T, root, filename string, size int) string {
	t.Helper()
	dir := filepath.Join(root, DirName)
	if err := os.MkdirAll(dir, 0o775); err != nil {
		t.Fatal(err)
	}
	path := filepath.Join(dir, filename)
	if err := os.WriteFile(path, make([]byte, size), 0o664); err != nil {
		t.Fatal(err)
	}
	return path
}

// later is a clock far enough ahead that everything written now reads as old.
func later() time.Time { return time.Now().Add(90 * 24 * time.Hour) }

// collected indexes a report's collectable entries by path.
func collected(report *GCReport) map[string]Collectable {
	out := map[string]Collectable{}
	for _, entry := range report.Collectable {
		out[entry.Path] = entry
	}
	return out
}

// retained indexes a report's retained entries by path.
func retained(report *GCReport) map[string]string {
	out := map[string]string{}
	for _, entry := range report.Retained {
		out[entry.Path] = entry.Reason
	}
	return out
}

// A published entry that does not validate is a repair job, not garbage, so
// even an ancient unreadable one is kept.
func TestGCKeepsAnEntryThatDoesNotValidate(t *testing.T) {
	root := t.TempDir()
	path := storeEntry(t, root, "star--2.7@abc.sqf", 1024)

	report, err := GC(GCOptions{Dirs: []string{root}, Grace: 30 * 24 * time.Hour, now: later()})
	if err != nil {
		t.Fatal(err)
	}
	if len(report.Collectable) != 0 {
		t.Fatalf("collectable = %v, want none", report.Collectable)
	}
	if reason := retained(report)[path]; !strings.Contains(reason, "validate") {
		t.Errorf("retained for %q, want a validation reason", reason)
	}
}

// Nothing inside the grace is collectable, whatever else is true of it.
func TestGCKeepsAnEntryInsideTheGrace(t *testing.T) {
	root := t.TempDir()
	path := storeEntry(t, root, "star--2.7@abc.sqf", 1024)

	report, err := GC(GCOptions{Dirs: []string{root}, Grace: 30 * 24 * time.Hour})
	if err != nil {
		t.Fatal(err)
	}
	if reason := retained(report)[path]; !strings.Contains(reason, "grace") {
		t.Errorf("retained for %q, want the grace reason", reason)
	}
}

// An abandoned .part is crash debris rather than an artifact, so it is reported
// without having to validate and on its own short fixed window.
func TestGCReportsAbandonedStagingWithoutValidating(t *testing.T) {
	root := t.TempDir()
	path := storeEntry(t, root, "star--2.7@abc.sqf.slurm-1"+producer.PreparedSuffix, 2048)

	report, err := GC(GCOptions{Dirs: []string{root}, Grace: 365 * 24 * time.Hour, now: later()})
	if err != nil {
		t.Fatal(err)
	}
	entry, ok := collected(report)[path]
	if !ok {
		t.Fatalf("staging file not reported: %+v", report)
	}
	if !entry.Staging {
		t.Error("staging file was not marked as such")
	}
	if report.Reclaimable != 2048 {
		t.Errorf("reclaimable = %d, want the staging file's size", report.Reclaimable)
	}
}

// A young .part may still be being written by a live producer.
func TestGCKeepsRecentStaging(t *testing.T) {
	root := t.TempDir()
	path := storeEntry(t, root, "star--2.7@abc.sqf.local-1"+producer.PreparedSuffix, 2048)

	report, err := GC(GCOptions{Dirs: []string{root}, Grace: 0})
	if err != nil {
		t.Fatal(err)
	}
	if _, ok := collected(report)[path]; ok {
		t.Error("a staging file younger than the fixed window was reported collectable")
	}
}

// A container reading an entry holds a shared lock, and the reason has to say
// so rather than blaming age or validity.
func TestGCKeepsAnEntryInUse(t *testing.T) {
	root := t.TempDir()
	path := storeEntry(t, root, "star--2.7@abc.sqf.local-1"+producer.PreparedSuffix, 1024)
	lock, err := image.AcquireLock(path, false)
	if err != nil {
		t.Fatal(err)
	}
	defer lock.Close() //nolint:errcheck

	report, err := GC(GCOptions{Dirs: []string{root}, now: later()})
	if err != nil {
		t.Fatal(err)
	}
	if _, ok := collected(report)[path]; ok {
		t.Fatal("an entry under a shared lock was reported collectable")
	}
	if reason := retained(report)[path]; !strings.Contains(reason, "in use") {
		t.Errorf("retained for %q, want the in-use reason", reason)
	}
}

// Clearing the write bit is how an identity is pinned, and it survives GC.
func TestGCKeepsAProtectedEntry(t *testing.T) {
	root := t.TempDir()
	path := storeEntry(t, root, "star--2.7@abc.sqf.local-1"+producer.PreparedSuffix, 1024)
	if err := os.Chmod(path, 0o444); err != nil {
		t.Fatal(err)
	}
	t.Cleanup(func() { os.Chmod(path, 0o664) }) //nolint:errcheck

	report, err := GC(GCOptions{Dirs: []string{root}, now: later()})
	if err != nil {
		t.Fatal(err)
	}
	if _, ok := collected(report)[path]; ok {
		t.Fatal("a write-protected entry was reported collectable")
	}
	if reason := retained(report)[path]; !strings.Contains(reason, "protected") {
		t.Errorf("retained for %q, want the protection reason", reason)
	}
}

// Reporting is the default, so a plain run must leave everything on disk.
func TestGCReportsWithoutDeleting(t *testing.T) {
	root := t.TempDir()
	path := storeEntry(t, root, "star--2.7@abc.sqf.local-1"+producer.PreparedSuffix, 4096)

	report, err := GC(GCOptions{Dirs: []string{root}, now: later()})
	if err != nil {
		t.Fatal(err)
	}
	if len(report.Collectable) != 1 || report.Applied {
		t.Fatalf("report = %+v, want one collectable and no apply", report)
	}
	if _, err := os.Stat(path); err != nil {
		t.Fatalf("a report deleted something: %v", err)
	}
}

func TestGCApplyDeletesAndReportsWhatItFreed(t *testing.T) {
	root := t.TempDir()
	path := storeEntry(t, root, "star--2.7@abc.sqf.local-1"+producer.PreparedSuffix, 4096)

	report, err := GC(GCOptions{Dirs: []string{root}, Apply: true, now: later()})
	if err != nil {
		t.Fatal(err)
	}
	if len(report.Collectable) != 1 || !report.Collectable[0].Removed {
		t.Fatalf("report = %+v, want one removed entry", report)
	}
	if report.Reclaimable != 4096 {
		t.Errorf("reclaimable = %d, want 4096", report.Reclaimable)
	}
	if _, err := os.Stat(path); !os.IsNotExist(err) {
		t.Errorf("entry survived --apply: %v", err)
	}
}

// Deleting across every writable tier by omission is the one thing --apply may
// not do, and it is refused before anything is opened.
func TestGCApplyRefusesAnUnscopedRun(t *testing.T) {
	report, err := GC(GCOptions{Apply: true})
	if !errors.Is(err, ErrUnscoped) {
		t.Fatalf("err = %v, want ErrUnscoped", err)
	}
	if report != nil {
		t.Errorf("report = %+v, want none", report)
	}
}

// Largest first, because the point of the report is deciding what is worth
// reclaiming.
func TestGCSortsBySizeDescending(t *testing.T) {
	root := t.TempDir()
	storeEntry(t, root, "small@a.sqf.local-1"+producer.PreparedSuffix, 1024)
	storeEntry(t, root, "large@b.sqf.local-1"+producer.PreparedSuffix, 8192)

	report, err := GC(GCOptions{Dirs: []string{root}, now: later()})
	if err != nil {
		t.Fatal(err)
	}
	if len(report.Collectable) != 2 {
		t.Fatalf("collectable = %+v, want two", report.Collectable)
	}
	if report.Collectable[0].Size < report.Collectable[1].Size {
		t.Errorf("sorted ascending: %d then %d", report.Collectable[0].Size, report.Collectable[1].Size)
	}
}

// A root with no store/ is skipped, not given one.
func TestGCCreatesNoStore(t *testing.T) {
	root := t.TempDir()
	if _, err := GC(GCOptions{Dirs: []string{root}, now: later()}); err != nil {
		t.Fatal(err)
	}
	if _, err := os.Stat(filepath.Join(root, DirName)); !os.IsNotExist(err) {
		t.Errorf("gc created a store directory: %v", err)
	}
}

// A flat artifact answers to a bare name and belongs to `remove`; taking one
// out of service through identity housekeeping would be a surprise.
func TestRemoveRefusesAFlatArtifact(t *testing.T) {
	root := t.TempDir()
	identity := strings.Repeat("a", 64)
	flat := filepath.Join(root, "star--2.7.sqf")
	if err := os.WriteFile(flat, []byte("payload"), 0o664); err != nil {
		t.Fatal(err)
	}
	previous := resolveIdentityFor
	resolveIdentityFor = func(string, IdentityQuery, []string) (Candidate, Report, error) {
		return Candidate{Name: "star/2.7", Path: flat, Root: root, Layout: LayoutFlat}, Report{}, nil
	}
	t.Cleanup(func() { resolveIdentityFor = previous })

	_, err := Remove("star/2.7", IdentityQuery{Scheme: "identity-v1", SHA256: identity}, []string{root})
	if !errors.Is(err, ErrNotStored) {
		t.Fatalf("err = %v, want ErrNotStored", err)
	}
	if _, statErr := os.Stat(flat); statErr != nil {
		t.Errorf("the flat artifact was removed anyway: %v", statErr)
	}
}

func TestRemoveDeletesAStoreEntry(t *testing.T) {
	root := t.TempDir()
	path := storeEntry(t, root, "star--2.7@abc.sqf", 1024)
	previous := resolveIdentityFor
	resolveIdentityFor = func(string, IdentityQuery, []string) (Candidate, Report, error) {
		return Candidate{Name: "star/2.7", Path: path, Root: root, Layout: LayoutStored}, Report{}, nil
	}
	t.Cleanup(func() { resolveIdentityFor = previous })

	if _, err := Remove("star/2.7", IdentityQuery{}, []string{root}); err != nil {
		t.Fatal(err)
	}
	if _, err := os.Stat(path); !os.IsNotExist(err) {
		t.Errorf("entry survived removal: %v", err)
	}
}

// Clearing the write bit pins an identity, and `store rm` honours it for the
// same reason GC does.
func TestRemoveRefusesAProtectedEntry(t *testing.T) {
	root := t.TempDir()
	path := storeEntry(t, root, "star--2.7@abc.sqf", 1024)
	if err := os.Chmod(path, 0o444); err != nil {
		t.Fatal(err)
	}
	t.Cleanup(func() { os.Chmod(path, 0o664) }) //nolint:errcheck
	previous := resolveIdentityFor
	resolveIdentityFor = func(string, IdentityQuery, []string) (Candidate, Report, error) {
		return Candidate{Name: "star/2.7", Path: path, Root: root, Layout: LayoutStored}, Report{}, nil
	}
	t.Cleanup(func() { resolveIdentityFor = previous })

	_, err := Remove("star/2.7", IdentityQuery{}, []string{root})
	if !errors.Is(err, image.ErrProtected) {
		t.Fatalf("err = %v, want ErrProtected", err)
	}
	if _, statErr := os.Stat(path); statErr != nil {
		t.Errorf("a protected entry was removed: %v", statErr)
	}
}

// timestamped writes a file with the given atime and mtime.
//
// ctime is always now and cannot be set from userspace, so it is the timestamp
// the other two have to beat — which is why the cases below put atime or mtime
// in the *future* rather than backdating. That is not contrived: it is the only
// way to drive the rule on a real filesystem, and it is exactly what a
// spuriously refreshed atime looks like.
func timestamped(t *testing.T, atime, mtime time.Time) os.FileInfo {
	t.Helper()
	path := filepath.Join(t.TempDir(), "entry.sqf")
	if err := os.WriteFile(path, []byte("payload"), 0o664); err != nil {
		t.Fatal(err)
	}
	if err := os.Chtimes(path, atime, mtime); err != nil {
		t.Fatal(err)
	}
	info, err := os.Lstat(path)
	if err != nil {
		t.Fatal(err)
	}
	return info
}

// atime is the signal that tracks use, so when it is the newest it decides.
func TestAgeOfPrefersAtimeWhenItIsNewest(t *testing.T) {
	now := time.Now().Add(24 * time.Hour)
	info := timestamped(t, now.Add(-time.Hour), now.Add(-10*time.Hour))

	age, basis := ageOf(info, now)
	if basis != BasisAtime {
		t.Fatalf("basis = %q, want atime", basis)
	}
	if age < 50*time.Minute || age > 70*time.Minute {
		t.Errorf("age = %s, want about an hour", age)
	}
}

func TestAgeOfPrefersMtimeWhenItIsNewest(t *testing.T) {
	now := time.Now().Add(24 * time.Hour)
	info := timestamped(t, now.Add(-10*time.Hour), now.Add(-time.Hour))

	if _, basis := ageOf(info, now); basis != BasisMtime {
		t.Fatalf("basis = %q, want mtime", basis)
	}
}

// The noatime case: atime frozen at creation and mtime older still, so ctime
// answers. That is no worse than reading ctime alone, which is the floor the
// rule is designed to have.
func TestAgeOfFallsBackToCtimeWhenAtimeIsFrozen(t *testing.T) {
	old := time.Now().Add(-90 * 24 * time.Hour)
	info := timestamped(t, old, old)

	age, basis := ageOf(info, time.Now())
	if basis != BasisCtime {
		t.Fatalf("basis = %q, want ctime", basis)
	}
	if age > time.Minute {
		t.Errorf("age = %s, want ctime's own age rather than the backdated stamps", age)
	}
}

// An NFS server whose clock runs ahead is not evidence that a file is from the
// future. Treating it as brand new retains it, which is the safe direction.
func TestAgeOfClampsClockSkewToZero(t *testing.T) {
	now := time.Now()
	info := timestamped(t, now.Add(time.Hour), now.Add(time.Hour))

	if age, _ := ageOf(info, now); age != 0 {
		t.Fatalf("age = %s, want 0 for a timestamp ahead of the clock", age)
	}
}

// The safety property the whole rule exists for: a backup or an indexer that
// reads everything makes entries look newer, and newer means kept.
func TestGCRetainsAnEntryWhoseAtimeWasRefreshed(t *testing.T) {
	root := t.TempDir()
	path := storeEntry(t, root, "star--2.7@abc.sqf.local-1"+producer.PreparedSuffix, 1024)
	clock := later()

	if report, err := GC(GCOptions{Dirs: []string{root}, now: clock}); err != nil {
		t.Fatal(err)
	} else if _, ok := collected(report)[path]; !ok {
		t.Fatal("the entry was not collectable before its atime moved")
	}

	// Something read it — a backup, an indexer, a tree walk. A read moves atime
	// and nothing else, which a zero mtime asks os.Chtimes to leave alone.
	if err := os.Chtimes(path, clock.Add(-time.Minute), time.Time{}); err != nil {
		t.Fatal(err)
	}
	report, err := GC(GCOptions{Dirs: []string{root}, now: clock})
	if err != nil {
		t.Fatal(err)
	}
	if _, ok := collected(report)[path]; ok {
		t.Fatal("a refreshed atime did not retain the entry")
	}
	if reason := retained(report)[path]; !strings.Contains(reason, "atime") {
		t.Errorf("retained for %q, want atime named as the deciding stamp", reason)
	}
}
