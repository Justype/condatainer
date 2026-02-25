package exec

import (
	"context"
	"fmt"
	"os"
	"sort"
	"strings"

	"github.com/Justype/condatainer/internal/apptainer"
	"github.com/Justype/condatainer/internal/container"
	"github.com/Justype/condatainer/internal/overlay"
	"github.com/Justype/condatainer/internal/utils"
)

// Run executes a command inside a configured Apptainer container.
func Run(ctx context.Context, options Options) error {
	options = options.ensureDefaults()

	if err := apptainer.SetBin(options.ApptainerBin); err != nil {
		return err
	}

	// Use shared container setup logic
	setupResult, err := container.Setup(container.SetupConfig{
		Overlays:       options.Overlays,
		WritableImg:    options.WritableImg,
		EnvSettings:    options.EnvSettings,
		BindPaths:      options.BindPaths,
		Fakeroot:       options.Fakeroot,
		ApptainerFlags: options.ApptainerFlags,
	})
	if err != nil {
		return err
	}

	// Auto-enable fakeroot if needed for writable .img
	fakeroot := container.AutoEnableFakeroot(setupResult.LastImg, options.WritableImg, setupResult.Fakeroot)

	if utils.DebugMode {
		utils.PrintDebug("[EXEC]Exec overlays: %v", setupResult.Overlays)
		utils.PrintDebug("[EXEC]Bind paths: %v", setupResult.BindPaths)
		utils.PrintDebug("[EXEC]Overlay mounts: %v", setupResult.OverlayArgs)
		utils.PrintDebug("[EXEC]Env list: %v", setupResult.EnvList)
		utils.PrintDebug("[EXEC]Command: %s", strings.Join(options.Command, " "))
		if len(options.EnvSettings) > 0 {
			utils.PrintDebug("[EXEC]Env overrides: %v", options.EnvSettings)
		}
	}

	// Print environment notes if interactive
	if utils.IsInteractiveShell() && !options.HideOutput && !options.HidePrompt {
		if len(setupResult.EnvNotes) > 0 {
			utils.PrintMessage("Overlay envs:")
			sortedNotes := make([]string, 0, len(setupResult.EnvNotes))
			for key := range setupResult.EnvNotes {
				sortedNotes = append(sortedNotes, key)
			}
			sort.Strings(sortedNotes)
			for _, key := range sortedNotes {
				fmt.Printf("  %s: %s\n", utils.StyleName(key), utils.StyleInfo(setupResult.EnvNotes[key]))
			}
			fmt.Println("")
		} else if setupResult.LastImg != "" && options.WritableImg {
			utils.PrintMessage("Overlay env:")
			fmt.Printf("  %s: %s\n", utils.StyleName("CNT_CONDA_PREFIX"), utils.StylePath("/ext3/env"))
			fmt.Println("")
		}
	}

	// Acquire read locks on all overlay files for the duration of exec.
	// This prevents concurrent remove/update from deleting files in use.
	var execLocks []*overlay.Lock
	releaseLocks := func() {
		for _, l := range execLocks {
			l.Close()
		}
	}
	for _, ol := range setupResult.Overlays {
		if !utils.FileExists(ol) || utils.DirExists(ol) {
			continue
		}
		// Skip .img files: Apptainer flocks them itself during execution, so
		// acquiring our own lock would conflict. They are pre-checked via
		// container.Setup() → CheckAvailable() before exec starts.
		if utils.IsImg(ol) {
			continue
		}
		lock, err := overlay.AcquireLock(ol, false) // .sqf: always read-only
		if err != nil {
			releaseLocks()
			return err
		}
		execLocks = append(execLocks, lock)
	}
	if utils.FileExists(options.BaseImage) {
		lock, err := overlay.AcquireLock(options.BaseImage, false)
		if err != nil {
			releaseLocks()
			return err
		}
		execLocks = append(execLocks, lock)
	}
	defer releaseLocks()

	opts := &apptainer.ExecOptions{
		Bind:       setupResult.BindPaths,
		Overlay:    setupResult.OverlayArgs,
		Env:        setupResult.EnvList,
		Fakeroot:   fakeroot,
		HideOutput: options.HideOutput,
		Additional: setupResult.ApptainerFlags,
	}

	// Pass through stdin if requested (for interactive build scripts)
	if options.PassThruStdin {
		if options.Stdin != nil {
			opts.Stdin = options.Stdin
		} else {
			opts.Stdin = os.Stdin
		}
	}

	if err := apptainer.Exec(ctx, options.BaseImage, options.Command, opts); err != nil {
		return err
	}
	return nil
}
