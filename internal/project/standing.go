package project

import (
	"errors"
	"fmt"
	"strings"

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
// artifact is an error naming the remedy, never a fetch, a build, or a
// fallback to whatever currently answers to the name — the failure a lock
// exists to prevent.
func (s *Standing) ResolveComplete(requests []lock.Request, opts ResolveOptions) (*Resolution, error) {
	resolution, err := Resolve(s.Root, s.Lock, requests, opts)
	if err != nil {
		return nil, err
	}
	if !resolution.Complete() {
		return nil, unresolvedError(s.Root, resolution)
	}
	return resolution, nil
}

// ResolveNames resolves each of names — parsed with the same grammar a
// `#DEP:` uses — against the standing project, in order.
//
// Every name must already be pinned: a caller reaching for this declares its
// own fixed requirements outside any `#DEP:`, so they reach the lock only as
// manual pins (see README, "Manual pins"). An unresolved name is therefore an
// error naming `project restore`, never a fall back to whatever currently
// answers to the name.
func (s *Standing) ResolveNames(names []string) ([]string, error) {
	if len(names) == 0 {
		return nil, nil
	}
	requests := make([]lock.Request, 0, len(names))
	for _, name := range names {
		request, reason := lock.ParseDeclaration(name)
		if reason != "" {
			return nil, fmt.Errorf("%s: %s", name, reason)
		}
		requests = append(requests, request)
	}
	resolution, err := s.ResolveComplete(requests, ResolveOptions{})
	if err != nil {
		return nil, err
	}
	paths := make([]string, len(resolution.Mounts))
	for i, mount := range resolution.Mounts {
		paths[i] = mount.Path
	}
	return paths, nil
}

// Base resolves the standing project's reserved root pin (lock.BaseKey) to a
// local path. Empty and no error when the project has no base pin at all —
// which `project lock` no longer leaves possible, but an older or hand-built
// lock might — so a caller falls through to its own ordinary default.
//
// Strict otherwise, like ResolveNames: an unresolved pin refuses naming
// `project restore` rather than silently falling back to this machine's
// configured default_distro. See README, "The project's root".
func (s *Standing) Base() (string, error) {
	if _, ok := s.Lock.Pins[lock.BaseKey]; !ok {
		return "", nil
	}
	resolution, err := s.ResolveComplete(
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
