package container

import (
	"os"
	"path/filepath"
	"strings"

	"github.com/Justype/condatainer/catalog"
	"github.com/Justype/condatainer/internal/artifact/meta"
	"github.com/Justype/condatainer/internal/utils"
)

// SnapshotLookup is what LookupSnapshot found beside a writable .img.
type SnapshotLookup struct {
	// Path is the env-typed snapshot to pair with the .img, or "" if none.
	Path string
	// Blocked is the candidate path when the first slot that exists on disk
	// is not an env-typed artifact — occupied by something unrelated, or a
	// mistake. Path is empty whenever Blocked is set.
	Blocked string
}

// LookupSnapshot finds the env-typed .sqf that pairs with a writable .img,
// checking a personal line before the shared one. It is the one place this
// pairing is computed — reused in both directions: autoloading a snapshot at
// mount time (container.Setup), and resolving a bare `overlay freeze`'s
// replacement target, which asks "which snapshot did this .img actually
// autoload against."
//
// Naming is derived from the .img's own basename, not hardcoded to "env":
// stripping a trailing "-<user>" suffix (if present) gives a stem, and the
// candidates are, in order:
//
//  1. <stem>-<user>.sqf — a personal snapshot line, distinct from the shared one.
//  2. <stem>.sqf — the shared snapshot.
//
// The first candidate that exists on disk decides the outcome: if its
// runtime.json reports Type == catalog.TypeEnv, it is the pair. Otherwise the
// slot is Blocked — something else is occupying the derived path, and this
// function does not fall through to the next candidate, since guessing that a
// different slot was meant is exactly what autoload must not do. A caller on
// the mount path treats Blocked the same as "nothing found" (the .img just
// mounts alone); `overlay freeze` treats it as a hard refusal.
func LookupSnapshot(imgPath string) SnapshotLookup {
	for _, candidate := range SnapshotCandidates(imgPath) {
		if !utils.FileExists(candidate) {
			continue
		}
		rt, err := meta.ReadRuntime(candidate)
		if err != nil || rt.Type != catalog.TypeEnv {
			return SnapshotLookup{Blocked: candidate}
		}
		return SnapshotLookup{Path: candidate}
	}
	return SnapshotLookup{}
}

// SnapshotCandidates returns the ordered paths LookupSnapshot considers for
// imgPath: the personal line (if $USER is set), then the shared line.
func SnapshotCandidates(imgPath string) []string {
	dir := filepath.Dir(imgPath)
	stem := snapshotStem(imgPath)
	var candidates []string
	if user := os.Getenv("USER"); user != "" {
		candidates = append(candidates, filepath.Join(dir, stem+"-"+user+".sqf"))
	}
	return append(candidates, filepath.Join(dir, stem+".sqf"))
}

// snapshotStem is an .img's basename with its extension and, if present, its
// trailing "-<user>" suffix removed — so "rnaseq-alice.img" and "rnaseq.img"
// both resolve against the same two candidates.
func snapshotStem(imgPath string) string {
	base := strings.TrimSuffix(filepath.Base(imgPath), filepath.Ext(imgPath))
	if user := os.Getenv("USER"); user != "" && strings.HasSuffix(base, "-"+user) {
		return strings.TrimSuffix(base, "-"+user)
	}
	return base
}

// PairedSize returns the combined size, in bytes, of path and — when path is
// a writable .img that autoloads one (LookupSnapshot) — its paired snapshot,
// plus the snapshot path itself ("" when there is none). Either file may be
// absent; a missing one simply contributes zero, which is what lets this
// double as the third-state probe: no .img on disk yet, but a snapshot's own
// size and path reported anyway.
func PairedSize(path string) (sizeBytes int64, snapshot string) {
	if info, err := os.Stat(path); err == nil {
		sizeBytes = info.Size()
	}
	if utils.IsImg(path) {
		if lookup := LookupSnapshot(path); lookup.Path != "" {
			snapshot = lookup.Path
			if info, err := os.Stat(snapshot); err == nil {
				sizeBytes += info.Size()
			}
		}
	}
	return sizeBytes, snapshot
}

// isEnvSnapshotSqf reports whether path is a .sqf recording Type == TypeEnv.
func isEnvSnapshotSqf(path string) bool {
	if !utils.IsSqf(path) {
		return false
	}
	rt, err := meta.ReadRuntime(path)
	return err == nil && rt.Type == catalog.TypeEnv
}
