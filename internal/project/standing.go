package project

import (
	"context"
	"errors"
	"fmt"
	"strings"

	"github.com/Justype/condatainer/internal/config"
	"github.com/Justype/condatainer/internal/project/lock"
)

// Standing is the project a caller is standing in, loaded once from its
// working directory. Every overlay- and base-resolution entry point in this
// package hangs off it, so `exec`/`run`'s `-o` handling, a helper's
// `#REQUIRED_OVERLAYS:`, and root resolution all read one lock the same way
// instead of three packages repeating "find the root, load the lock, handle
// no-project" on their own.
type Standing struct {
	Root string
	Lock *lock.Lock
}

// StandingAt resolves the project rooted at or above cwd — walking up through
// ancestors, so standing anywhere inside a project's tree is recognized, not
// only at its root — or (nil, nil) when no ancestor has one, the caller's cue
// to fall back to its own ordinary resolution unchanged.
func StandingAt(cwd string) (*Standing, error) {
	root, err := lock.RootAbove(cwd)
	if errors.Is(err, lock.ErrNoProject) {
		return nil, nil
	}
	if err != nil {
		return nil, err
	}
	current, err := lock.Load(root)
	if err != nil {
		return nil, err
	}
	return &Standing{Root: root, Lock: current}, nil
}

// ResolveComplete is Resolve against the standing project, refusing with an
// error naming `project restore` when anything comes back unresolved.
//
// Every caller reaching for this resolves rather than acquires: an absent
// artifact is an error naming the remedy, never a fetch or a build. Whether
// an unpinned name may fall back to whatever currently answers to it is
// opts.LiveResolve's call, made by the caller — see Resolve's own comment.
func (s *Standing) ResolveComplete(ctx context.Context, requests []lock.Request, opts ResolveOptions) (*Resolution, error) {
	resolution, err := Resolve(ctx, s.Root, s.Lock, requests, opts)
	if err != nil {
		return nil, err
	}
	if !resolution.Complete() {
		return nil, unresolvedError(s.Root, resolution)
	}
	return resolution, nil
}

// DefaultDistro is the distro a bare name expands under here: the project's
// selected root, else the configured default_distro.
func (s *Standing) DefaultDistro() string {
	if distro := s.SelectedDistro(); distro != "" {
		return distro
	}
	return config.ResolvedDefaultDistro()
}

// ResolveNames resolves each of names — parsed with the same grammar a
// `#DEP:` uses — the way `exec -o` does: a pinned name through its pin, any
// other name live against what is installed. It returns one Mount per name, in
// order; a name nothing installed answers has an empty Path, for the caller to
// acquire. A pinned artifact absent here, and an unpinned path, refuse naming
// `project restore`.
func (s *Standing) ResolveNames(ctx context.Context, names []string) ([]Mount, error) {
	if len(names) == 0 {
		return nil, nil
	}
	requests := make([]lock.Request, len(names))
	for i, name := range names {
		request, reason := lock.ParseDeclaration(name)
		if reason != "" {
			return nil, fmt.Errorf("%s: %s", name, reason)
		}
		requests[i] = request
	}
	resolution, err := Resolve(ctx, s.Root, s.Lock, requests, ResolveOptions{LiveResolve: true, Distro: s.DefaultDistro()})
	if err != nil {
		return nil, err
	}
	answered := make(map[string]Mount, len(resolution.Mounts))
	for _, mount := range resolution.Mounts {
		answered[mount.Request] = mount
	}
	var refused Resolution
	mounts := make([]Mount, len(requests))
	for i, request := range requests {
		if mount, ok := answered[request.Key]; ok {
			mounts[i] = mount
			continue
		}
		if _, _, pinned := lock.MatchPin(s.Lock, request); pinned || request.Kind != lock.KindName {
			refused.Unresolved = append(refused.Unresolved, unresolvedFor(resolution, request.Key))
		}
		mounts[i] = Mount{Request: request.Key}
	}
	if len(refused.Unresolved) > 0 {
		return nil, unresolvedError(s.Root, &refused)
	}
	return mounts, nil
}

// unresolvedFor is the entry of resolution for key.
func unresolvedFor(resolution *Resolution, key string) Unresolved {
	for _, unresolved := range resolution.Unresolved {
		if unresolved.Request == key {
			return unresolved
		}
	}
	return Unresolved{Request: key}
}

// Base resolves the standing project's reserved root pin (lock.BaseKey) to a
// local path. Empty and no error when the project has no base pin at all —
// which `project lock` no longer leaves possible, but an older or hand-built
// lock might — so a caller falls through to its own ordinary default.
//
// Strict, unlike ResolveNames: an unresolved pin refuses naming
// `project restore` rather than silently falling back to this machine's
// configured default_distro. See README, "The project's root".
func (s *Standing) Base(ctx context.Context) (string, error) {
	if _, ok := s.Lock.Pins[lock.BaseKey]; !ok {
		return "", nil
	}
	resolution, err := s.ResolveComplete(ctx,
		[]lock.Request{{Key: lock.BaseKey, Kind: lock.KindName}}, ResolveOptions{})
	if err != nil {
		return "", err
	}
	return resolution.Mounts[0].Path, nil
}

// SelectedDistro reports the distro named by the standing project's base pin,
// or "" when there is no base pin or its lookup fails for any reason — the
// safe default for a completion or display path, which falls back to a
// configured default rather than erroring.
//
// Reads the vendored manifest only, never the local images roots: unlike
// Base this never has to resolve a local copy, so it answers the same in a
// fresh clone that has restored nothing.
func (s *Standing) SelectedDistro() string {
	pin, ok := s.Lock.Pins[lock.BaseKey]
	if !ok {
		return ""
	}
	verified, _ := lock.Verify(s.Root, s.Lock)
	entry, ok := verified.Entries[pin.Artifact]
	if !ok {
		return ""
	}
	// select-distro and DeriveBase only ever compose "<distro>/base", so this
	// always has one slash to cut on.
	distro, _, found := strings.Cut(entry.Manifest.Name, "/")
	if !found {
		return ""
	}
	return distro
}

// unresolvedError explains why a resolution came back incomplete and names
// the one remedy. Shared wording for every caller standing in a project: an
// unresolved declaration always means the same thing, whether it came from a
// script scan, a manual -o, a helper's required overlays, or the reserved
// base pin.
func unresolvedError(root string, resolution *Resolution) error {
	var out strings.Builder
	fmt.Fprintf(&out, "this project cannot supply what was asked of it:")
	for _, unresolved := range resolution.Unresolved {
		fmt.Fprintf(&out, "\n  %s: %s", unresolved.Request, unresolved.Reason)
	}
	fmt.Fprintf(&out, "\n\nrun `condatainer project restore --project %s` to make them available", root)
	return errors.New(out.String())
}
