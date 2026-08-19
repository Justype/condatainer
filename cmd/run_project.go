package cmd

import (
	"errors"
	"fmt"
	"os"
	"path/filepath"
	"strings"

	"github.com/Justype/condatainer/internal/project"
	"github.com/Justype/condatainer/internal/project/lock"
	"github.com/Justype/condatainer/internal/scheduler"
	"github.com/Justype/condatainer/internal/utils"
)

// projectContext is a resolved project run: what to mount, and where the job
// must run. Nil when the command was not run from a project root.
type projectContext struct {
	Root string
	// Overlays are absolute paths, in declaration order.
	Overlays []string
}

// projectRunContext resolves a script's declarations through its project lock,
// or reports that there is no project to resolve against.
//
// The project is the current directory, and only when cnt-lock/ is directly in
// it. The script may live anywhere — `run scripts/align.sh` from the root is an
// ordinary project run — because the root, not the script's directory, is what
// relative declarations and the job's working directory resolve against.
func projectRunContext(contentScript string, specs *scheduler.ScriptSpecs) (*projectContext, error) {
	cwd, err := os.Getwd()
	if err != nil {
		return nil, err
	}
	root, err := lock.RootAt(cwd)
	if errors.Is(err, lock.ErrNoProject) {
		return nil, nil
	}
	if err != nil {
		return nil, err
	}
	script, err := filepath.Abs(contentScript)
	if err != nil {
		return nil, err
	}

	// A submitted job runs where the scheduler puts it, which is the submission
	// directory under SLURM and LSF and $HOME under PBS. Inside a project that
	// is the root, and a directive naming anywhere else is refused rather than
	// overridden — see §13.2.
	workDir, err := project.WorkDir(root, specs.Control.WorkDir)
	if err != nil {
		return nil, err
	}
	specs.Control.WorkDir = workDir

	current, err := lock.Load(root)
	if err != nil {
		return nil, err
	}
	scanned, err := lock.ScanScript(root, script)
	if err != nil {
		return nil, err
	}
	resolution, err := project.Resolve(root, current, scanned.Requests, project.ResolveOptions{})
	if err != nil {
		return nil, err
	}
	if !resolution.Complete() {
		return nil, unresolvedError(root, resolution)
	}

	context := &projectContext{Root: root}
	for _, mount := range resolution.Mounts {
		context.Overlays = append(context.Overlays, mount.Path)
		if mount.Found != "" {
			utils.PrintNote("%s is mounted from an equivalent artifact, not %s",
				utils.StyleName(mount.Name), short(mount.Identity))
		}
	}
	return context, nil
}

// unresolvedError explains why a project cannot run and names the one remedy.
//
// Inside a project this is always fatal. `run` never falls back to whatever
// currently answers to the name, and never restores on its own: a compute node
// routinely lacks the network, credentials, build tools and writable images
// directory that acquiring would need.
func unresolvedError(root string, resolution *project.Resolution) error {
	var out strings.Builder
	fmt.Fprintf(&out, "this project cannot run as locked:")
	for _, unresolved := range resolution.Unresolved {
		fmt.Fprintf(&out, "\n  %s: %s", unresolved.Request, unresolved.Reason)
	}
	fmt.Fprintf(&out, "\n\nrun `condatainer project restore --project %s` to make them available", root)
	return errors.New(out.String())
}
