package project

import (
	"path"
	"path/filepath"
	"sort"
	"strings"

	"github.com/Justype/condatainer/internal/artifact/key"
	"github.com/Justype/condatainer/internal/helperhistory"
	"github.com/Justype/condatainer/internal/project/lock"
	"github.com/Justype/condatainer/internal/runtime/container"
	"github.com/Justype/condatainer/internal/utils"
)

// ProjectStatus is a plain, read-only view of a project — the CLI's
// `project status` and the dashboard's status panel both read this, so the
// two never describe a project's state two different ways.
type ProjectStatus struct {
	HasProject bool   `json:"has_project"`
	Root       string `json:"root,omitempty"`
	PinCount   int    `json:"pin_count"`
	// Unpinned lists #DEP:-derived declarations with no pin yet. "Unpinned"
	// to match this package's own vocabulary (README, "Manual pins"), not
	// "drift".
	Unpinned []string `json:"unpinned,omitempty"`
	// UnpinnedHelperOverlays and ManualPinUsage are the same two
	// usageIndex-derived signals project/orchestrate.Lock returns.
	UnpinnedHelperOverlays []string            `json:"unpinned_helper_overlays,omitempty"`
	ManualPinUsage         map[string][]string `json:"manual_pin_usage,omitempty"`
}

// Status reads a project's current state without writing anything, walking
// up from cwd the way `Standing` always does. HasProject is false, and every
// other field zero, when no ancestor of cwd has one — the caller's cue that
// there is nothing here to show.
//
// For an ambient caller with no `--project` concept of its own — the
// dashboard's raw cwd field is the only one today — walking up is exactly
// right, the same reason a helper's #REQUIRED_OVERLAYS: resolves this way.
// The CLI's `project status` is not that caller: it sits in the `project`
// subcommand family, which resolves `--project`/ambient cwd through the
// strict, non-walking `lock.RootAt` (see README, "Acting in a project") —
// it calls StatusAt on its own already-resolved root instead.
func Status(cwd string) (*ProjectStatus, error) {
	standing, err := StandingAt(cwd)
	if err != nil {
		return nil, err
	}
	if standing == nil {
		return &ProjectStatus{}, nil
	}
	return statusAt(standing.Root, standing.Lock)
}

// StatusAt reads status for the project rooted exactly at root — no walking
// up to an ancestor. root must already contain cnt-lock/; ErrNoProject
// otherwise (matching lock.RootAt).
func StatusAt(root string) (*ProjectStatus, error) {
	if _, err := lock.RootAt(root); err != nil {
		return nil, err
	}
	l, err := lock.Load(root)
	if err != nil {
		return nil, err
	}
	return statusAt(root, l)
}

// statusAt is Status/StatusAt's shared body once a root and its loaded lock
// are already in hand — the only difference between the two is how that
// root was found.
func statusAt(root string, l *lock.Lock) (*ProjectStatus, error) {
	scanned, err := lock.Scan(root, lock.ScanOptions{})
	if err != nil {
		return nil, err
	}
	unpinnedOverlays, err := UnpinnedHelperOverlays(root, l)
	if err != nil {
		return nil, err
	}
	usage, err := ManualPinUsage(root, l)
	if err != nil {
		return nil, err
	}
	return &ProjectStatus{
		HasProject:             true,
		Root:                   root,
		PinCount:               len(l.Pins),
		Unpinned:               unpinnedRequests(l, scanned),
		UnpinnedHelperOverlays: unpinnedOverlays,
		ManualPinUsage:         usage,
	}, nil
}

// unpinnedRequests lists pinnable declarations with no pin yet, without
// mutating the lock the way lock.Reconcile does — Status only ever reads.
func unpinnedRequests(l *lock.Lock, scanned *lock.ScanResult) []string {
	var out []string
	for _, request := range scanned.Requests {
		if !request.Kind.Pinnable() {
			continue
		}
		if _, _, ok := lock.MatchPin(l, request); !ok {
			out = append(out, request.Key)
		}
	}
	sort.Strings(out)
	return out
}

// UnpublishedFrozenEnv lists pinned frozen-environment artifacts with no
// recorded remote — what `project lock`'s old noteUnpublished checked for,
// minus the printing: each caller (CLI text, SSE stream) decides how to
// show it.
func UnpublishedFrozenEnv(l *lock.Lock, pinned []*lock.Pinned) []string {
	var out []string
	for _, p := range pinned {
		if p.Identity.Scheme == string(key.SnapshotEnvV1) && len(l.Remotes[p.Artifact]) == 0 {
			out = append(out, p.Artifact)
		}
	}
	return out
}

// usageIndex walks helperhistory.ListAll(root) and maps each classified
// overlay key to the sorted, deduplicated helper names that have some
// recorded combination containing it. Built once and read two ways —
// UnpinnedHelperOverlays and ManualPinUsage are its only two callers.
//
// Classification is the same three-way split lock/scan.go's Kind makes,
// since only two of its four kinds are ever pinnable — an entry that can't
// be pinned contributes no key at all, in either direction.
func usageIndex(root string) (map[string][]string, error) {
	all, err := helperhistory.ListAll(root)
	if err != nil {
		return nil, err
	}

	sets := map[string]map[string]bool{} // pin key -> set of helper names
	add := func(pinKey, helperName string) {
		set := sets[pinKey]
		if set == nil {
			set = map[string]bool{}
			sets[pinKey] = set
		}
		set[helperName] = true
	}

	for helperName, byLocation := range all {
		for location, combos := range byLocation {
			for _, combo := range combos {
				for _, overlay := range combo.Overlays {
					if pinKey, ok := classifyOverlay(overlay, location); ok {
						add(pinKey, helperName)
					}
				}
			}
			// The frozen env.sqf autoloaded at this location, if any, is
			// folded in for free — it isn't a fourth, special case, just
			// another KindPath candidate once resolved to a path.
			if pinKey, ok := frozenEnvKey(root, location); ok {
				add(pinKey, helperName)
			}
		}
	}

	out := make(map[string][]string, len(sets))
	for pinKey, set := range sets {
		names := make([]string, 0, len(set))
		for name := range set {
			names = append(names, name)
		}
		sort.Strings(names)
		out[pinKey] = names
	}
	return out, nil
}

// classifyOverlay computes overlay's pin-key form, given the
// (root-relative) location it was recorded from, or reports it contributes
// no key at all when nothing could ever pin it.
//
//   - A catalog name/version (normalizeOverlayForHistory's internal form —
//     no overlay extension) → the key is the overlay itself.
//   - A project-relative path that stays under the root once joined with
//     location → a KindPath key.
//   - An absolute path, a writable .img, or a joined path that still
//     escapes the root → KindExternal/KindWritable — neither is ever
//     pinnable, so neither contributes a key.
func classifyOverlay(overlay, location string) (string, bool) {
	if !utils.IsOverlay(overlay) {
		return overlay, true
	}
	if utils.IsImg(overlay) || filepath.IsAbs(overlay) {
		return "", false
	}
	candidate := path.Clean(path.Join(location, filepath.ToSlash(overlay)))
	if candidate == ".." || strings.HasPrefix(candidate, "../") {
		return "", false
	}
	return lock.PathPrefix + candidate, true
}

// frozenEnvKey resolves location's shared, frozen env.sqf — never the
// writable .img, which has no identity to pin, and never a personal
// env-$USER.sqf, which has no audience beyond whoever created it — to the
// same KindPath key classifyOverlay would give an ordinary recorded overlay
// at that path. It isn't a fourth, special case, once it's the frozen form.
func frozenEnvKey(root, location string) (string, bool) {
	dir := filepath.Join(root, filepath.FromSlash(location))
	resolved := container.ResolveEnvOverlay("", dir)
	if resolved == "" || utils.IsImg(resolved) || filepath.Base(resolved) != "env.sqf" {
		return "", false
	}
	rel, err := filepath.Rel(root, resolved)
	if err != nil {
		return "", false
	}
	rel = filepath.ToSlash(rel)
	if rel == ".." || strings.HasPrefix(rel, "../") {
		return "", false
	}
	return lock.PathPrefix + path.Clean(rel), true
}

// UnpinnedHelperOverlays is usageIndex's keys with no matching entry in
// l.Pins — a suggestion, never a fallback. It never writes a pin and never
// becomes a scan Finding: usage is still never evidence of a project
// dependency.
func UnpinnedHelperOverlays(root string, l *lock.Lock) ([]string, error) {
	index, err := usageIndex(root)
	if err != nil {
		return nil, err
	}
	var out []string
	for pinKey := range index {
		if _, pinned := l.Pins[pinKey]; !pinned {
			out = append(out, pinKey)
		}
	}
	sort.Strings(out)
	return out, nil
}

// ManualPinUsage reports, for every manual pin, which helpers have a
// recorded combination using it — a fact, never a verdict. An absent (nil)
// entry never means "unused," only "no recorded helper usage": a cold-start
// or infrequently run helper looks the same either way, and only a person
// weighing which helpers actually matter can tell those apart.
func ManualPinUsage(root string, l *lock.Lock) (map[string][]string, error) {
	index, err := usageIndex(root)
	if err != nil {
		return nil, err
	}
	out := make(map[string][]string, len(l.Pins))
	for pinKey, entry := range l.Pins {
		if !entry.Manual {
			continue
		}
		out[pinKey] = index[pinKey]
	}
	return out, nil
}
