package cmd

import (
	"context"
	"errors"
	"os"
	"path/filepath"
	"strings"

	"github.com/Justype/condatainer/internal/overlay"
	"github.com/Justype/condatainer/internal/utils"
	"github.com/spf13/cobra"
)

// oCmd is a shortcut for "overlay create"
var oCmd = &cobra.Command{
	Use:   "o [image_path]",
	Short: "Shortcut for 'overlay create'",
	Long: `Creates an ext3 overlay image optimized for specific workloads.
This is a shortcut for the 'overlay create' command.

If no image path is provided, defaults to 'env.img'.`,
	Example: `  condatainer o # 10G with default inode ratio
  condatainer o my_data.img -s 50G -t data
  condatainer o --fakeroot --sparse
  condatainer o -f environment.yml # Initialize with conda env file`,

	Args: cobra.RangeArgs(0, 1),

	Run: func(cmd *cobra.Command, args []string) {
		// 1. Handle Positional Argument (Image Path)
		path := "env.img"
		if len(args) > 0 {
			path = args[0]
		}

		// Auto-append .img extension if not present and path doesn't have an extension
		if !utils.IsImg(path) {
			if strings.Contains(filepath.Base(path), ".") {
				ExitWithUsageError("Overlay image must have a .img extension.")
			}
			path += ".img"
		}

		// Convert to absolute path
		absPath, err := filepath.Abs(path)
		if err != nil {
			ExitWithError("Failed to resolve path: %v", err)
		}
		path = absPath

		if utils.FileExists(path) || utils.DirExists(path) {
			ExitWithUsageError("Path %s already exists.", utils.StylePath(path))
		}

		// 2. Parse Flags
		sizeStr, _ := cmd.Flags().GetString("size")
		fakeroot, _ := cmd.Flags().GetBool("fakeroot")
		sparse, _ := cmd.Flags().GetBool("sparse")
		typeFlag, _ := cmd.Flags().GetString("type")
		envFile, _ := cmd.Flags().GetString("file")
		noTmp, _ := cmd.Flags().GetBool("no-tmp")

		// 3. Parse size
		sizeMB, err := utils.ParseSizeToMB(sizeStr)
		if err != nil {
			ExitWithError("Invalid size: %v", err)
		}

		uid, gid := os.Getuid(), os.Getgid()
		if fakeroot {
			uid, gid = 0, 0
		}
		opts := &overlay.CreateOptions{
			Path:           path,
			SizeMB:         sizeMB,
			UID:            uid,
			GID:            gid,
			Profile:        overlay.GetProfile(typeFlag),
			Sparse:         sparse,
			FilesystemType: "ext3",
		}

		if noTmp {
			// 4a. Create directly at target path (no tmp), then conda init there.
			if err := overlay.CreateDirectly(cmd.Context(), opts); err != nil {
				if errors.Is(err, context.Canceled) || errors.Is(cmd.Context().Err(), context.Canceled) {
					utils.PrintWarning("Overlay creation cancelled.")
					return
				}
				ExitWithError("%v", err)
			}
			if err := initializeOverlayWithConda(cmd.Context(), path, envFile, fakeroot); err != nil {
				os.Remove(path)
				if errors.Is(err, context.Canceled) || errors.Is(cmd.Context().Err(), context.Canceled) {
					utils.PrintWarning("Overlay initialization cancelled.")
					return
				}
				ExitWithError("Failed to initialize overlay with conda environment: %v", err)
			}
		} else {
			// 4b. Create sparse at local tmp (fast I/O), conda init there, then move + allocate.
			tmpPath, err := overlay.CreateInTmp(cmd.Context(), opts)
			if err != nil {
				if errors.Is(err, context.Canceled) || errors.Is(cmd.Context().Err(), context.Canceled) {
					utils.PrintWarning("Overlay creation cancelled.")
					return
				}
				ExitWithError("%v", err)
			}

			if err := initializeOverlayWithConda(cmd.Context(), tmpPath, envFile, fakeroot); err != nil {
				os.Remove(tmpPath)
				if errors.Is(err, context.Canceled) || errors.Is(cmd.Context().Err(), context.Canceled) {
					utils.PrintWarning("Overlay initialization cancelled.")
					return
				}
				ExitWithError("Failed to initialize overlay with conda environment: %v", err)
			}

			utils.PrintMessage("Moving overlay to %s", utils.StylePath(path))
			copied, err := overlay.MoveOverlayCopied(tmpPath, path, sparse)
			if err != nil {
				os.Remove(tmpPath)
				ExitWithError("Failed to move overlay to destination: %v", err)
			}

			// Skip AllocateOverlay when io.Copy was used: zeros already written physically.
			if !sparse && !copied {
				overlay.AllocateOverlay(cmd.Context(), path, sizeMB)
			}
		}

		utils.PrintSuccess("Created overlay %s", utils.StylePath(path))
	},
}

func init() {
	// 1. Attach to root command
	rootCmd.AddCommand(oCmd)

	// 2. Define flags (same as overlay create)
	oCmd.Flags().StringP("size", "s", "10G", "Set overlay size (e.g., 500M, 10G)")
	oCmd.Flags().StringP("type", "t", "balanced", "Overlay profile: small/balanced/large files")
	oCmd.Flags().Bool("fakeroot", false, "Create a fakeroot-compatible overlay (owned by root)")
	oCmd.Flags().BoolP("sparse", "S", false, "Create a sparse overlay image (no pre-allocation)")
	oCmd.Flags().StringP("file", "f", "", "Initialize with Conda environment file (.yml or .yaml)")
	oCmd.Flags().Bool("no-tmp", false, "Create directly at target path instead of local tmp (slower on network filesystems)")

	oCmd.RegisterFlagCompletionFunc("type", func(cmd *cobra.Command, args []string, toComplete string) ([]string, cobra.ShellCompDirective) {
		opts := []string{"small", "balanced", "large"}
		res := make([]string, 0, len(opts))
		for _, o := range opts {
			if toComplete == "" || strings.HasPrefix(o, toComplete) {
				res = append(res, o)
			}
		}
		return res, cobra.ShellCompDirectiveNoFileComp
	})

	oCmd.RegisterFlagCompletionFunc("file", func(cmd *cobra.Command, args []string, toComplete string) ([]string, cobra.ShellCompDirective) {
		return nil, cobra.ShellCompDirectiveFilterFileExt
	})
}
