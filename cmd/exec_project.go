package cmd

import (
	"context"
	"os"
	"strings"

	"github.com/Justype/condatainer/internal/config"
	"github.com/Justype/condatainer/internal/project"
)

// standingProject reports the project the caller is standing in, or nil when
// there is none or --no-project was given.
func standingProject() (*project.Standing, error) {
	if noProjectRequested {
		return nil, nil
	}
	cwd, err := os.Getwd()
	if err != nil {
		return nil, err
	}
	standing, err := project.StandingAt(cwd)
	if err != nil || standing == nil {
		return nil, err
	}
	announceProject(standing.Root)
	return standing, nil
}

// projectOverlays resolves `-o` arguments: through the lock of the project the
// caller is standing in, else against what is installed. `exec` and `e` share
// it. There is no flag either way: `--project DIR` belongs to the commands
// that act *on* a project, and standing somewhere else is the opt-out.
func projectOverlays(ctx context.Context, overlays []string) ([]string, error) {
	if len(overlays) == 0 {
		return overlays, nil
	}
	standing, err := standingProject()
	if err != nil {
		return nil, err
	}
	return resolveOverlayValues(ctx, overlays, standing, false)
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
// Strict, unlike projectOverlays: the root is never a bare or partial name,
// so there is nothing for LiveResolve to apply to. An unresolved pin is a
// refusal naming `project restore`, never a silent fall back to this
// machine's default_distro. See internal/project/README.md, "The project's root".
func projectBaseImage(ctx context.Context) (string, error) {
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
	return standing.Base(ctx)
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
