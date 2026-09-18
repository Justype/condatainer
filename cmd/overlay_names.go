package cmd

import (
	"context"
	"fmt"

	"github.com/Justype/condatainer/catalog"
	"github.com/Justype/condatainer/internal/image"
	"github.com/Justype/condatainer/internal/project"
	"github.com/Justype/condatainer/internal/project/lock"
	"github.com/Justype/condatainer/internal/utils"
)

// resolveOverlayValues turns each overlay value (a name, possibly bare or
// partial, or a file, with an optional :ro/:rw suffix) into what the runtime
// mounts. Standing in a project (standing non-nil) every value answers to its
// lock, and a name no pin answers resolves live; with no project a name
// resolves against what is installed and a file is its own answer. Either way a
// version constraint is refused, since it is a build-recipe feature — unless
// allowConstraint, where such a value is left as it was.
//
// A name nothing installed answers, and any file, is returned unchanged
// outside a project for container.ResolveOverlayPaths to report; inside one it
// is an error naming `project restore`.
func resolveOverlayValues(ctx context.Context, values []string, standing *project.Standing, allowConstraint bool) ([]string, error) {
	requests := make([]lock.Request, 0, len(values))
	suffixes := make([]string, 0, len(values))
	kept := make([]int, 0, len(values))
	out := append([]string(nil), values...)
	for i, value := range values {
		name, suffix := splitOverlayMode(value)
		if name == "" {
			continue
		}
		request, reason := lock.ParseDeclaration(name)
		if reason != "" {
			if dep, err := catalog.ParseDep(name); allowConstraint && err == nil && dep.Op != "" {
				continue
			}
			return nil, fmt.Errorf("%s: %s", name, reason)
		}
		requests = append(requests, request)
		suffixes = append(suffixes, suffix)
		kept = append(kept, i)
	}
	if len(requests) == 0 {
		return out, nil
	}

	opts := project.ResolveOptions{LiveResolve: true, Distro: projectDefaultDistro()}
	var mounts []project.Mount
	if standing != nil {
		resolution, err := standing.ResolveComplete(ctx, requests, opts)
		if err != nil {
			return nil, err
		}
		mounts = resolution.Mounts
	} else {
		var err error
		if mounts, err = project.ResolveUnlocked(ctx, requests, opts); err != nil {
			return nil, err
		}
	}
	if len(mounts) != len(requests) {
		return nil, fmt.Errorf("resolved %d of %d overlays", len(mounts), len(requests))
	}

	for n, mount := range mounts {
		// Outside a project a file or an unanswered name is the resolver's to report.
		if standing == nil && (mount.Unpinned || mount.Path == "") {
			continue
		}
		out[kept[n]] = mount.Path + suffixes[n]
		switch {
		case mount.Found != "":
			utils.PrintNote("%s is mounted from an equivalent artifact, not %s",
				utils.StyleName(mount.Name), short(mount.Identity))
		case mount.Live && standing != nil:
			utils.PrintNote("%s resolved to %s, not pinned — run `condatainer project pin %s` to lock it",
				utils.StyleName(mount.Request), utils.StyleName(mount.Name), mount.Name)
		case mount.Live && mount.Name != mount.Request:
			utils.PrintNote("Expanding '%s' to '%s'", mount.Request, mount.Name)
		}
	}
	return out, nil
}

// installedOverlayPath resolves a name — bare, partial or exact — to the
// installed overlay it means, for the commands that act on one overlay.
func installedOverlayPath(name string) (string, bool, error) {
	scan, err := image.ScanOverlays(image.ScanOptions{})
	if err != nil {
		utils.PrintWarning("%v", err)
	}
	_, path, found, err := catalog.SolveInstalled(context.Background(), scan, projectDefaultDistro(), name)
	return path, found, err
}
