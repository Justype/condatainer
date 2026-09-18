// Package orchestrate composes internal/project and internal/project/lock
// steps that a single caller — `cmd/project.go`'s `project lock`, and the
// dashboard's lock action — needs run in a fixed sequence. It exists as its
// own package because internal/project/publish already imports
// internal/project (for LookupAt/LookupLocal), so internal/project itself
// can never import internal/project/publish back without a cycle; this
// package sits above both instead.
package orchestrate

import (
	"context"

	"github.com/Justype/condatainer/internal/config"
	"github.com/Justype/condatainer/internal/logging"
	"github.com/Justype/condatainer/internal/project"
	"github.com/Justype/condatainer/internal/project/lock"
	"github.com/Justype/condatainer/internal/project/publish"
)

// LockOptions tunes Lock. Empty today — a placeholder so this function's
// signature does not need to change if it later grows an option.
type LockOptions struct{}

// LockResult is what one `project lock` run produced, for the caller to
// report however it likes — CLI text, JSON, or an SSE stream.
type LockResult struct {
	Root    string
	Current *lock.Lock
	Pinned  []*lock.Pinned
	Failed  []error
	Scanned *lock.ScanResult
	// Upstream maps a vendored artifact to the repository recordUpstream
	// found already serving its exact identity, for display.
	Upstream map[string]string
	// UnpublishedFrozenEnv lists pinned frozen-environment artifacts with no
	// recorded remote yet — noteUnpublished's old condition, minus the
	// printing.
	UnpublishedFrozenEnv []string
	// UnpinnedHelperOverlays and ManualPinUsage are project's two
	// usageIndex-derived signals: overlays a helper has used but nobody has
	// pinned, and which helpers are recorded using each manual pin.
	UnpinnedHelperOverlays []string
	ManualPinUsage         map[string][]string
}

// Lock creates or updates root's lock from its #DEP: declarations: scan,
// reconcile, pin, derive the base, and report. This is the orchestration
// `cmd/project.go`'s `project lock` used to inline directly inside its
// RunE — extracted so the dashboard's lock action is a second caller of
// exactly this sequence, never a reimplementation of it.
func Lock(ctx context.Context, root string, opts LockOptions) (*LockResult, error) {
	current, err := lock.Load(root)
	if err != nil {
		return nil, err
	}
	scanned, err := lock.Scan(root, lock.ScanOptions{})
	if err != nil {
		return nil, err
	}
	unpinned := lock.Reconcile(root, current, scanned)
	if err := lock.Publish(root, current); err != nil {
		return nil, err
	}

	// Pinning is what makes this a lock rather than a scan: a declaration
	// names what is needed, and the lock has to say which exact build
	// answers it.
	pinned, failed := lock.PinAll(root, current, unpinned, lock.PinOptions{})
	upstream := map[string]string{}
	for _, p := range pinned {
		for artifact, remote := range RecordUpstream(ctx, root, current, p) {
			upstream[artifact] = remote
		}
	}
	// Every project's root is pinned, unconditionally — there is no closure
	// to walk deciding whether one is needed. A manual `project
	// select-distro` override is left alone.
	if basePinned, err := lock.DeriveBase(root, current, config.ResolvedDefaultDistro(), lock.PinOptions{}); err != nil {
		failed = append(failed, err)
	} else if basePinned != nil {
		for artifact, remote := range RecordUpstream(ctx, root, current, basePinned) {
			upstream[artifact] = remote
		}
		pinned = append(pinned, basePinned)
	}
	unpublished := project.UnpublishedFrozenEnv(current, pinned)
	if len(pinned) > 0 {
		if err := lock.Publish(root, current); err != nil {
			return nil, err
		}
	}

	unpinnedOverlays, err := project.UnpinnedHelperOverlays(root, current)
	if err != nil {
		return nil, err
	}
	usage, err := project.ManualPinUsage(root, current)
	if err != nil {
		return nil, err
	}

	return &LockResult{
		Root: root, Current: current, Pinned: pinned, Failed: failed, Scanned: scanned,
		Upstream: upstream, UnpublishedFrozenEnv: unpublished,
		UnpinnedHelperOverlays: unpinnedOverlays, ManualPinUsage: usage,
	}, nil
}

// RecordUpstream adds a fetch location for every artifact this pin
// vendored that its own recipe collection already publishes at the exact
// same identity, and reports which, for display.
//
// Moved here from cmd/project.go (it used to be unexported and local): it
// mutates the lock (l.AddRemote) in a way Publish persists, so it is part
// of what locking produces, not CLI decoration — the dashboard's lock
// action has to write the same cnt-lock/lock.json content the CLI would for
// the same run. Exported because `project pin` and `project select-distro`
// call it too, each pinning outside a full Lock run.
//
// It never fails the pin. No network, no configured source, no declared
// endpoint and no match all record nothing — locking has to work offline,
// and an absent remote costs a rebuild rather than an error.
func RecordUpstream(ctx context.Context, root string, l *lock.Lock, pinned *lock.Pinned) map[string]string {
	cat, err := config.OpenCatalog(ctx)
	if err != nil {
		logging.FromContext(ctx).Debug("no catalog, so no upstream locations were recorded", "err", err)
		return nil
	}
	found := publish.Upstream(ctx, root, pinned.Vendored, cat)
	shown := make(map[string]string, len(found))
	for artifact, remotes := range found {
		for _, remote := range remotes {
			if err := l.AddRemote(artifact, remote); err != nil {
				logging.FromContext(ctx).Debug("could not record an upstream location", "artifact", artifact, "err", err)
				continue
			}
			shown[artifact] = remote.Repository
		}
	}
	return shown
}
