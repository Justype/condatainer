package cmd

import (
	"errors"
	"fmt"
	"os"
	"strings"

	"github.com/Justype/condatainer/internal/project"
	"github.com/Justype/condatainer/internal/project/lock"
	"github.com/Justype/condatainer/internal/utils"
)

// projectOverlays resolves `-o` arguments through the lock of the project the
// caller is standing in, returning them unchanged when there is no project.
//
// `exec` and `e` share it, so both act *in* whatever project the caller is
// standing in. There is no flag either way: `--project DIR` belongs to the
// commands that act *on* a project, and standing somewhere else is the opt-out.
//
// It resolves rather than acquires — an artifact that is absent is an error
// naming `project restore`, never a fetch or a build, and never a fallback to
// whatever currently answers to the name. That fallback is the failure a lock
// exists to prevent.
//
// This deliberately does not live in container.ResolveOverlayPaths, though that
// is the one place a name becomes a path. Five callers there resolve a build's
// own `#DEP:` or a base image, and none of them may pick up the lock of whatever
// directory the user happened to be standing in.
func projectOverlays(overlays []string) ([]string, error) {
	if len(overlays) == 0 {
		return overlays, nil
	}
	cwd, err := os.Getwd()
	if err != nil {
		return nil, err
	}
	root, err := lock.RootAt(cwd)
	if errors.Is(err, lock.ErrNoProject) {
		return overlays, nil
	}
	if err != nil {
		return nil, err
	}
	current, err := lock.Load(root)
	if err != nil {
		return nil, err
	}

	// The suffix is a mount mode, not part of the name, and it has to survive
	// resolution to reach the runtime.
	requests := make([]lock.Request, 0, len(overlays))
	suffixes := make([]string, 0, len(overlays))
	for _, overlay := range overlays {
		value, suffix := splitOverlayMode(overlay)
		request, reason := lock.ParseDeclaration(value)
		if reason != "" {
			return nil, fmt.Errorf("-o %s: %s", overlay, reason)
		}
		requests = append(requests, request)
		suffixes = append(suffixes, suffix)
	}

	resolution, err := project.Resolve(root, current, requests, project.ResolveOptions{})
	if err != nil {
		return nil, err
	}
	if !resolution.Complete() {
		return nil, unresolvedError(root, resolution)
	}
	// One mount per request, in request order, so the suffixes line up.
	if len(resolution.Mounts) != len(suffixes) {
		return nil, fmt.Errorf("resolved %d of %d overlays", len(resolution.Mounts), len(suffixes))
	}

	utils.PrintMessage("Project: %s", utils.StylePath(root))
	resolved := make([]string, 0, len(resolution.Mounts))
	for i, mount := range resolution.Mounts {
		resolved = append(resolved, mount.Path+suffixes[i])
		if mount.Found != "" {
			utils.PrintNote("%s is mounted from an equivalent artifact, not %s",
				utils.StyleName(mount.Name), short(mount.Identity))
		}
	}
	return resolved, nil
}

// splitOverlayMode separates a `:ro`/`:rw` mount mode from what it applies to.
func splitOverlayMode(overlay string) (value, suffix string) {
	for _, mode := range []string{":ro", ":rw"} {
		if trimmed, ok := strings.CutSuffix(overlay, mode); ok {
			return trimmed, mode
		}
	}
	return overlay, ""
}
