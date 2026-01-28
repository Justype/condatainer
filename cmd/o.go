package cmd

import (
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
				utils.PrintError("Overlay image must have a .img extension.")
				os.Exit(1)
			}
			path += ".img"
		}

		// Convert to absolute path
		absPath, err := filepath.Abs(path)
		if err != nil {
			utils.PrintError("Failed to resolve path: %v", err)
			os.Exit(1)
		}
		path = absPath

		if utils.FileExists(path) || utils.DirExists(path) {
			utils.PrintError("Path %s already exists.", utils.StylePath(path))
			os.Exit(1)
		}

		// 2. Parse Flags
		sizeStr, _ := cmd.Flags().GetString("size")
		fakeroot, _ := cmd.Flags().GetBool("fakeroot")
		sparse, _ := cmd.Flags().GetBool("sparse")
		typeFlag, _ := cmd.Flags().GetString("type")
		fsType, _ := cmd.Flags().GetString("fs")
		envFile, _ := cmd.Flags().GetString("file")

		// 3. Parse size
		sizeMB, err := utils.ParseSizeToMB(sizeStr)
		if err != nil {
			utils.PrintError("Invalid size: %v", err)
			os.Exit(1)
		}

		// 4. Create the Overlay
		if fakeroot {
			err = overlay.CreateForRoot(path, sizeMB, typeFlag, sparse, fsType)
		} else {
			err = overlay.CreateForCurrentUser(path, sizeMB, typeFlag, sparse, fsType)
		}

		if err != nil {
			utils.PrintError("%v", err)
			os.Exit(1)
		}

		// 5. Initialize with Conda environment if file specified
		if envFile != "" {
			if err := initializeOverlayWithConda(path, envFile, fakeroot); err != nil {
				utils.PrintError("Failed to initialize overlay with conda environment: %v", err)
				os.Exit(1)
			}
		} else {
			// Initialize with minimal conda environment (zlib)
			if err := initializeOverlayWithConda(path, "", fakeroot); err != nil {
				utils.PrintError("Failed to initialize overlay with conda environment: %v", err)
				os.Exit(1)
			}
		}
	},
}

func init() {
	// 1. Attach to root command
	rootCmd.AddCommand(oCmd)

	// 2. Define flags (same as overlay create)
	oCmd.Flags().StringP("size", "s", "10G", "Set overlay size (e.g., 500M, 10G)")
	oCmd.Flags().StringP("type", "t", "balanced", "Overlay profile: small/balanced/large files")
	oCmd.Flags().String("fs", "ext3", "Filesystem type: ext3 or ext4 (some system does not support ext4)")
	oCmd.Flags().Bool("fakeroot", false, "Create a fakeroot-compatible overlay (owned by root)")
	oCmd.Flags().Bool("sparse", false, "Create a sparse overlay image")
	oCmd.Flags().StringP("file", "f", "", "Initialize with Conda environment file (.yml or .yaml)")
}
