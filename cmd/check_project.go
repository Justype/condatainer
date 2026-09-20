package cmd

import (
	"context"
	"fmt"
	"os"
	"path/filepath"

	"github.com/Justype/condatainer/internal/project"
	"github.com/Justype/condatainer/internal/project/lock"
	"github.com/Justype/condatainer/internal/utils"
)

// projectCheck reports whether one script can run right now against the project
// in the current directory, and returns false when there is no project here.
//
// The project commands answer for the whole lock: `project validate --installed`
// and `project restore --dry-run` ask whether every pin is here. This answers
// "can I run *this script*", through the same resolution and the same anchor
// `run` uses —
// including its live resolve of an unpinned #DEP: name — so the two cannot
// disagree about whether a script is runnable.
func projectCheck(ctx context.Context, scriptPaths []string, metaDeps []string) (handled bool, err error) {
	// A `<name>` addresses a recipe and belongs to no project, so checking one
	// stays ordinary wherever it is typed.
	if len(scriptPaths) == 0 || noProjectRequested {
		return false, nil
	}
	cwd, err := os.Getwd()
	if err != nil {
		return false, err
	}
	standing, err := project.StandingAt(cwd)
	if err != nil {
		return true, err
	}
	if standing == nil {
		return false, nil
	}
	root := standing.Root
	announceProject(root)
	if err := oneScript(root, scriptPaths, metaDeps); err != nil {
		return true, err
	}
	script := scriptPaths[0]
	if checkAutoInstall {
		// -a installs by name into a shared images directory, which is one
		// checkout's dependency silently replacing the artifact every other user
		// of that directory sees. This holds even when nothing is pinned yet:
		// the answer is to lock them, not to install them by name.
		return true, fmt.Errorf(
			"-a installs by name and cannot run inside a project, which pins exact identities\n"+
				"run `condatainer project restore --project %s` instead", root)
	}

	scanned, err := lock.ScanScript(root, script)
	if err != nil {
		return true, err
	}
	resolution, err := project.Resolve(ctx, root, standing.Lock, scanned.Requests,
		project.ResolveOptions{LiveResolve: true, Distro: projectDefaultDistro()})
	if err != nil {
		return true, err
	}

	for _, mount := range resolution.Mounts {
		switch {
		case mount.Unpinned:
			utils.PrintMessage("  %s %s", utils.StyleWarning("cannot be pinned"), mount.Request)
		default:
			utils.PrintMessage("  %s %s → %s", utils.StyleSuccess("✓"), mount.Request, mount.Path)
		}
	}
	for _, unresolved := range resolution.Unresolved {
		utils.PrintError("%s: %s", utils.StyleName(unresolved.Request), unresolved.Reason)
	}
	if resolution.Complete() {
		utils.PrintSuccess("%s can run: %d declaration(s) resolved.", filepath.Base(script), len(resolution.Mounts))
		return true, nil
	}
	utils.PrintHint("Run %s to make them available.",
		utils.StyleAction("condatainer project restore --project "+root))
	return true, fmt.Errorf("%d declaration(s) cannot be resolved", len(resolution.Unresolved))
}

// oneScript refuses anything this command cannot answer per script.
//
// check merges every argument's declarations into one set and answers once,
// which is right when the answer is "install these by name". In a project the
// answer is per script — `run a.sh` mounts a.sh's pins and `run b.sh`
// mounts b.sh's — so a merged answer would describe a set nobody runs. The
// merged question is `project restore --dry-run`, which already reports per
// artifact.
func oneScript(root string, scriptPaths, metaDeps []string) error {
	if len(metaDeps) > 0 {
		return fmt.Errorf("cannot mix a project script with a `<name>`, which addresses a recipe and belongs to no project\n"+
			"check the project with `condatainer project restore --project %s --dry-run`", root)
	}
	if len(scriptPaths) > 1 {
		return fmt.Errorf("in a project, check answers for one script at a time, because each script mounts its own pins\n"+
			"for the whole project run `condatainer project restore --project %s --dry-run`", root)
	}
	return nil
}
