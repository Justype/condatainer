package cmd

import (
	"fmt"
	"os"
	"strings"

	"github.com/Justype/condatainer/internal/config"
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
func projectOverlays(overlays []string) ([]string, error) {
	if len(overlays) == 0 || noProjectRequested {
		return overlays, nil
	}
	cwd, err := os.Getwd()
	if err != nil {
		return nil, err
	}
	standing, err := project.StandingAt(cwd)
	if err != nil || standing == nil {
		return overlays, err
	}
	announceProject(standing.Root)

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

	resolution, err := standing.ResolveComplete(requests, project.ResolveOptions{})
	if err != nil {
		return nil, err
	}
	// One mount per request, in request order, so the suffixes line up.
	if len(resolution.Mounts) != len(suffixes) {
		return nil, fmt.Errorf("resolved %d of %d overlays", len(resolution.Mounts), len(suffixes))
	}

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

// projectDefaultDistro reports the distro every bare-name alias expands
// against: the standing project's selected root if there is one, else the
// configured default_distro. Every site that expands "<name>" to
// "<distro>/<name>" — completion, display, and expandBareName's build target
// — reads this instead of config.ResolvedDefaultDistro() directly, so
// `project select-distro` changes what all of them mean at once.
func projectDefaultDistro() string {
	if distro := projectSelectedDistro(); distro != "" {
		return distro
	}
	return config.ResolvedDefaultDistro()
}

// projectSelectedDistro reports the distro named by the standing project's
// base pin, or "" when there is no project, no base pin, or its lookup fails
// for any reason — the safe default for a completion or display path, which
// falls back to the configured default_distro rather than erroring.
func projectSelectedDistro() string {
	cwd, err := os.Getwd()
	if err != nil {
		return ""
	}
	standing, err := project.StandingAt(cwd)
	if err != nil || standing == nil {
		return ""
	}
	return standing.SelectedDistro()
}

// projectBaseImage resolves the project's locked root to a local path, for a
// caller standing in a project that did not itself request a root — "" and no
// error otherwise, so ensureRootBaseImage falls through to the ordinary
// configured default.
//
// Strict like projectOverlays: an unresolved pin is a refusal naming
// `project restore`, never a silent fall back to this machine's
// default_distro. See internal/project/README.md, "The project's root".
func projectBaseImage() (string, error) {
	if noProjectRequested {
		return "", nil
	}
	cwd, err := os.Getwd()
	if err != nil {
		return "", err
	}
	standing, err := project.StandingAt(cwd)
	if err != nil || standing == nil {
		return "", err
	}
	announceProject(standing.Root)
	return standing.Base()
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
