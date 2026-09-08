package cmd

import (
	"os"
	"path/filepath"

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
	script, err := filepath.Abs(contentScript)
	if err != nil {
		return nil, err
	}

	// A submitted job runs where the scheduler puts it, which is the submission
	// directory under SLURM and LSF and $HOME under PBS. Inside a project that
	// is the root, and a directive naming anywhere else is refused rather than
	// overridden.
	workDir, err := project.WorkDir(standing.Root, specs.Control.WorkDir)
	if err != nil {
		return nil, err
	}
	specs.Control.WorkDir = workDir

	scanned, err := lock.ScanScript(standing.Root, script)
	if err != nil {
		return nil, err
	}
	resolution, err := standing.ResolveComplete(scanned.Requests, project.ResolveOptions{})
	if err != nil {
		return nil, err
	}

	context := &projectContext{Root: standing.Root}
	for _, mount := range resolution.Mounts {
		context.Overlays = append(context.Overlays, mount.Path)
		if mount.Found != "" {
			utils.PrintNote("%s is mounted from an equivalent artifact, not %s",
				utils.StyleName(mount.Name), short(mount.Identity))
		}
	}
	return context, nil
}
