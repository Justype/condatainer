package cmd

import (
	"errors"
	"fmt"
	"os"
	"path/filepath"
	"strings"

	"github.com/Justype/condatainer/internal/config"
	"github.com/Justype/condatainer/internal/image"
	"github.com/Justype/condatainer/internal/image/freeze"
	"github.com/Justype/condatainer/internal/utils"
	"github.com/spf13/cobra"
)

var (
	freezeDescription string
	freezeBlockSize   string
	freezeUseTmp      bool
	freezeCompFlags   map[string]*bool
)

const overlayFreezeHelp = `Pack a writable overlay into an immutable .sqf artifact.

  - The environment is kept exactly as it is, not rebuilt.
  - It has no recipe behind it, so the .sqf is the only copy: back it up, or push it.
  - 'overlay unfreeze' turns one back into a writable .img.`

var overlayFreezeCmd = &cobra.Command{
	Use:   "freeze <overlay.img> [artifact.sqf]",
	Short: "Pack a writable overlay into an immutable artifact",
	Long:  overlayFreezeHelp,
	Args:  cobra.RangeArgs(1, 2),
	Example: `  condatainer overlay freeze env.img                     # → env.sqf beside it
  condatainer overlay freeze env.img ./overlays/env.sqf  # → that path
  condatainer overlay freeze env.img --zstd-high         # → pack harder than the default`,
	Run: func(cmd *cobra.Command, args []string) {
		ctx := cmd.Context()
		source := args[0]

		if !utils.IsImg(source) {
			ExitWithError("%s is not a writable overlay; freeze takes the .img an environment lives in.", source)
		}
		dest := ""
		if len(args) > 1 {
			dest = args[1]
		}
		// The base is only used to translate a directory-level deletion into
		// individual whiteouts against what that directory held; without one,
		// freeze still runs, just without that translation.
		base, err := config.GetBaseImage()
		if err != nil {
			utils.PrintWarning("no base image found (%v); a directory deleted wholesale (not file by file) will not be recorded as deleted, and its original contents will come back when this artifact is used", err)
		}
		// Held across the whole freeze rather than probed: a payload being written
		// to has no defined content to pack, and the pack takes minutes, so a
		// check that released would only prove the overlay was idle when the
		// command started. Shared is the whole requirement — it conflicts with
		// the exclusive lock a writer takes — and it opens read-only, so a
		// pinned overlay stays freezable.
		held, err := image.AcquireLock(source, false)
		if err != nil {
			ExitWithError("%v", err)
		}
		defer held.Close() //nolint:errcheck

		// Compression and block size come from the build settings, so a freeze
		// packs the way a build does; the flags override them for this pack only,
		// exactly as create's do.
		compressArgs := config.Global.Build.CompressArgs
		if args, err := compressArgsFromFlags(freezeCompFlags); err != nil {
			ExitWithError("%v", err)
		} else if args != "" {
			compressArgs = args
		}
		blockSize := config.Global.Build.BlockSize
		if freezeBlockSize != "" {
			if !config.IsValidBlockSize(freezeBlockSize) {
				ExitWithError("Invalid --block-size %q: must be a power of two between 4096 and 1M (e.g. 64k, 128k, 512k, 1m)", freezeBlockSize)
			}
			blockSize = freezeBlockSize
		}

		target, err := resolveFreezeTarget(source, dest)
		if err != nil {
			ExitWithError("%v", err)
		}
		if utils.FileExists(target) {
			ExitWithError("%s already exists.", target)
		}

		res, freezeErr := freeze.Freeze(ctx, freeze.Options{
			Image:        source,
			Target:       target,
			Description:  freezeDescription,
			Base:         base,
			CompressArgs: compressArgs,
			BlockSize:    blockSize,
			UseTmp:       freezeUseTmp,
			Processors:   config.Global.Build.Defaults.CpusPerTask,
		})
		if freezeErr != nil {
			if errors.Is(freezeErr, freeze.ErrEmptyOverlay) {
				ExitWithError("%s has no payload to freeze.", source)
			}
			ExitWithError("%v", freezeErr)
		}
		reportFreeze(res, source)
	},
}

// resolveFreezeTarget decides where the artifact is written — only where, since
// every frozen environment is called meta.EnvName. With no destination it lands
// beside the source with the extension changed.
func resolveFreezeTarget(source, dest string) (string, error) {
	named := strings.TrimSuffix(filepath.Base(source), filepath.Ext(source)) + ".sqf"
	switch {
	case dest == "":
		dest = filepath.Join(filepath.Dir(source), named)
	case strings.HasSuffix(dest, string(filepath.Separator)):
		// A trailing separator says directory whether or not one is there yet.
		// Without this a missing one falls through to the extension branch and
		// becomes a file named ".sqf".
		dest = filepath.Join(dest, named)
	default:
		if info, err := os.Stat(dest); err == nil && info.IsDir() {
			// A directory is a place to put it, not the artifact itself.
			dest = filepath.Join(dest, named)
		} else if !strings.HasSuffix(dest, ".sqf") {
			dest += ".sqf"
		}
	}

	abs, err := filepath.Abs(dest)
	if err != nil {
		return "", err
	}
	// The pack binds the artifact's directory into the container, so a missing one
	// surfaces there as a mount failure, after the walk and the translation have
	// already run.
	if info, err := os.Stat(filepath.Dir(abs)); err != nil || !info.IsDir() {
		return "", fmt.Errorf("%s does not exist; create it or choose another destination", filepath.Dir(abs))
	}
	if dir := imagesDirContaining(abs); dir != "" {
		return "", fmt.Errorf("%s is inside the images directory %s; a frozen environment belongs to its project and is addressed by path, so it is never filed by name", abs, dir)
	}
	return abs, nil
}

// imagesDirContaining reports which images directory holds path, if any. An
// images directory addresses artifacts by filename, and every frozen environment
// is called env, so two would collide there.
func imagesDirContaining(path string) string {
	for _, dir := range config.GetImageSearchPaths() {
		abs, err := filepath.Abs(dir)
		if err != nil {
			continue
		}
		if rel, err := filepath.Rel(abs, path); err == nil &&
			rel != ".." && !strings.HasPrefix(rel, ".."+string(filepath.Separator)) {
			return abs
		}
	}
	return ""
}

// reportFreeze prints what was produced. What the artifact contains is what
// `info` is for, read back from the file itself.
func reportFreeze(res freeze.Result, source string) {
	utils.PrintSuccess("frozen %s → %s", utils.StylePath(source), utils.StylePath(res.Path))
	if n := res.Translation.Deletions(); n > 0 {
		utils.PrintMessage("  carries %d deletion(s), translated from %s markers",
			n, string(res.Translation.Convention))
	}
}

func init() {
	overlayCmd.AddCommand(overlayFreezeCmd)
	f := overlayFreezeCmd.Flags()
	f.StringVarP(&freezeDescription, "description", "d", "", "Description recorded in the artifact")
	f.StringVar(&freezeBlockSize, "block-size", "", "SquashFS block size (default: build.block_size)")
	f.BoolVar(&freezeUseTmp, "use-tmp", false, "Copy to temp dir and pack (faster; needs payload-sized space)")
	freezeCompFlags = make(map[string]*bool, len(config.CompressOptions))
	for _, opt := range config.CompressOptions {
		freezeCompFlags[opt.Name] = f.Bool(opt.Name, false, opt.Description)
	}
	overlayFreezeCmd.RegisterFlagCompletionFunc("block-size", //nolint:errcheck
		func(cmd *cobra.Command, args []string, toComplete string) ([]string, cobra.ShellCompDirective) {
			return config.BlockSizeCompletions, cobra.ShellCompDirectiveNoFileComp
		})
}
