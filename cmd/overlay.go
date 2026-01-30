package cmd

import (
	"fmt"
	"os"
	"path/filepath"
	"strings"

	"github.com/Justype/condatainer/internal/apptainer"
	"github.com/Justype/condatainer/internal/exec"
	"github.com/Justype/condatainer/internal/overlay"
	"github.com/Justype/condatainer/internal/utils"
	"github.com/spf13/cobra"
)

// overlayCmd is the parent command for overlay management
var overlayCmd = &cobra.Command{
	Use:   "overlay",
	Short: "Manage persistent overlay images (create, resize, check, info)",
	Long: `Utilities to manage ext3 overlay images for Apptainer.
Allows creating images, resizing existing images,
and verifying filesystem integrity.`,
}

// ---------------------------------------------------------
// 1. Create Command
// ---------------------------------------------------------

var overlayCreateCmd = &cobra.Command{
	Use:   "create [image_path]",
	Short: "Create a new sparse overlay image",
	Long: `Creates an ext3 overlay image optimized for specific workloads.
	
If no image path is provided, defaults to 'env.img'.`,
	Example: `  condatainer overlay create # 10G with default inode ratio
  condatainer overlay create my_data.img -s 50G -t data
  condatainer overlay create --fakeroot --sparse`,

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

		sizeMB, err := utils.ParseSizeToMB(sizeStr)
		if err != nil {
			utils.PrintError("Invalid size format '%s': %v", sizeStr, err)
			os.Exit(1)
		}

		// 3. Create the Overlay
		if fakeroot {
			err = overlay.CreateForRoot(path, sizeMB, typeFlag, sparse, fsType, false)
		} else {
			err = overlay.CreateForCurrentUser(path, sizeMB, typeFlag, sparse, fsType, false)
		}

		if err != nil {
			utils.PrintError("%v", err)
			os.Exit(1)
		}

		// 4. Initialize with Conda environment if file specified
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

// ---------------------------------------------------------
// 2. Resize Command
// ---------------------------------------------------------

var resizeCmd = &cobra.Command{
	Use:   "resize [image_path]",
	Short: "Expand or shrink an existing overlay image",
	Example: `  condatainer overlay resize env.img -s 20G
  condatainer overlay resize data.img --size 512M`,
	Args: cobra.ExactArgs(1),

	// Enable Smart Tab Completion for .img files
	ValidArgsFunction: completeImages,

	Run: func(cmd *cobra.Command, args []string) {
		path := args[0]

		// 1. Get Size Flag
		sizeStr, _ := cmd.Flags().GetString("size")
		if sizeStr == "" {
			utils.PrintError("Size flag is required (e.g., -s 20G)")
			_ = cmd.Usage()
			os.Exit(1)
		}

		// 2. Parse Size
		sizeMB, err := utils.ParseSizeToMB(sizeStr)
		if err != nil {
			utils.PrintError("Invalid size format '%s': %v", sizeStr, err)
			os.Exit(1)
		}

		// 3. Execute
		err = overlay.Resize(path, sizeMB)
		if err != nil {
			utils.PrintError("%v", err)
			os.Exit(1)
		}
	},
}

// ---------------------------------------------------------
// 3. Info Command
// ---------------------------------------------------------

var infoCmd = &cobra.Command{
	Use:   "info [path]",
	Short: "Display disk usage and filesystem stats",
	Args:  cobra.ExactArgs(1),

	// Enable Smart Tab Completion for .img files
	ValidArgsFunction: completeImages,

	Run: func(cmd *cobra.Command, args []string) {
		path := args[0]

		stats, err := overlay.GetStats(path)
		if err != nil {
			utils.PrintError("%v", err)
			os.Exit(1)
		}

		usedBytes, diskPct := stats.Usage()
		totalBytes := stats.TotalBlocks * stats.BlockSize
		freeBytes := stats.FreeBlocks * stats.BlockSize
		reservedBytes := stats.ReservedBlocks * stats.BlockSize
		usedInodes := stats.TotalInodes - stats.FreeInodes
		inodePct := stats.InodeUsage()

		// File info
		fmt.Println(utils.StyleTitle("File"))
		fmt.Printf("  %-15s %s\n", "Path:", utils.StylePath(path))
		fmt.Printf("  %-15s %s\n", "Size:", utils.FormatBytes(stats.FileSizeBytes))
		if stats.IsSparse {
			fmt.Printf("  %-15s %s (%s on disk)\n", "Type:", utils.StyleInfo("Sparse"), utils.FormatBytes(stats.FileBlocksUsed))
		} else {
			fmt.Printf("  %-15s %s\n", "Type:", utils.StyleInfo("Allocated"))
		}

		// Filesystem info
		fmt.Println(utils.StyleTitle("Filesystem"))
		fmt.Printf("  %-15s %s\n", "Format:", utils.StyleInfo(stats.FilesystemType))
		fmt.Printf("  %-15s %s\n", "State:", stats.FilesystemState)
		fmt.Printf("  %-15s %d bytes\n", "Block Size:", stats.BlockSize)
		fmt.Printf("  %-15s %s\n", "UUID:", stats.FilesystemUUID)
		fmt.Printf("  %-15s %s\n", "Created:", stats.CreatedTime)
		if stats.LastMounted != "" && stats.LastMounted != "<not available>" {
			fmt.Printf("  %-15s %s\n", "Last Mounted:", stats.LastMounted)
		}

		// Ownership
		fmt.Println(utils.StyleTitle("Ownership"))
		if stats.UpperUID == 0 && stats.UpperGID == 0 {
			fmt.Printf("  %-15s %s (use with --fakeroot)\n", "Owner:", utils.StyleInfo("root"))
		} else if stats.UpperUID >= 0 {
			fmt.Printf("  %-15s UID=%d GID=%d\n", "Owner:", stats.UpperUID, stats.UpperGID)
		} else {
			fmt.Printf("  %-15s %s\n", "Owner:", "unknown")
		}

		// Disk Usage
		fmt.Println(utils.StyleTitle("Disk Usage"))
		fmt.Printf("  %-15s %s / %s (%.2f%%)\n", "Used:",
			utils.FormatBytes(usedBytes),
			utils.FormatBytes(totalBytes),
			diskPct)
		fmt.Printf("  %-15s %s\n", "Free:", utils.FormatBytes(freeBytes))
		if stats.ReservedBlocks > 0 {
			reservedPct := float64(stats.ReservedBlocks) / float64(stats.TotalBlocks) * 100
			fmt.Printf("  %-15s %s (%.1f%% for root)\n", "Reserved:", utils.FormatBytes(reservedBytes), reservedPct)
		}

		// Inode Usage
		fmt.Println(utils.StyleTitle("Inode Usage"))
		fmt.Printf("  %-15s %d / %d (%.2f%%)\n", "Used:",
			usedInodes,
			stats.TotalInodes,
			inodePct)
		fmt.Printf("  %-15s %d\n", "Free:", stats.FreeInodes)
	},
}

// ---------------------------------------------------------
// 4. Check Command
// ---------------------------------------------------------

var checkCmd = &cobra.Command{
	Use:   "check [path]",
	Short: "Verify filesystem integrity (e2fsck)",
	Args:  cobra.ExactArgs(1),

	// Enable Smart Tab Completion for .img files
	ValidArgsFunction: completeImages,

	Run: func(cmd *cobra.Command, args []string) {
		path := args[0]
		force, _ := cmd.Flags().GetBool("force")

		err := overlay.CheckIntegrity(path, force)
		if err != nil {
			utils.PrintError("%v", err)
			os.Exit(1)
		}
	},
}

// ---------------------------------------------------------
// 5. Chown Command
// ---------------------------------------------------------

var chownCmd = &cobra.Command{
	Use:   "chown [image_path]",
	Short: "Recursively set internal files to specific UID/GID",
	Long: `Walks the overlay filesystem (without mounting) and updates the UID/GID.

Defaults to the current user and the '/ext3' directory inside the image.
Use --root to force ownership to 0:0.`,
	Example: `  condatainer overlay chown env.img                  # Set /ext3 to current user
  condatainer overlay chown env.img --root -p /      # Set entire image to root
  condatainer overlay chown env.img -u 1001 -g 1001  # Set to specific ID`,
	Args: cobra.ExactArgs(1),

	// Enable Smart Tab Completion for .img files
	ValidArgsFunction: completeImages,

	Run: func(cmd *cobra.Command, args []string) {
		path := args[0]

		// Parse Flags
		isRoot, _ := cmd.Flags().GetBool("root")
		uidFlag, _ := cmd.Flags().GetInt("uid")
		gidFlag, _ := cmd.Flags().GetInt("gid")
		internalPath, _ := cmd.Flags().GetString("path")

		var targetUID, targetGID int
		currentUser := os.Getuid()
		currentGroup := os.Getgid()

		// Logic: Root Flag > Specific Flag > Current User (Default)
		if isRoot {
			targetUID = 0
			targetGID = 0
			utils.PrintDebug("Mode: Root (0:0)")
		} else {
			if uidFlag == -1 {
				targetUID = currentUser
			} else {
				targetUID = uidFlag
			}

			if gidFlag == -1 {
				targetGID = currentGroup
			} else {
				targetGID = gidFlag
			}
		}

		// Feedback to user
		utils.PrintDebug("Chown Target: %s inside %s", internalPath, path)
		utils.PrintDebug("New Owner:    UID=%d GID=%d", targetUID, targetGID)

		// Perform the recursive chown
		err := overlay.ChownRecursively(path, targetUID, targetGID, internalPath)
		if err != nil {
			utils.PrintError("%v", err)
			os.Exit(1)
		}
	},
}

// ---------------------------------------------------------
// Helper: Shell Completion
// ---------------------------------------------------------

// completeImages tells the shell to only suggest file extensions ending in "img".
func completeImages(cmd *cobra.Command, args []string, toComplete string) ([]string, cobra.ShellCompDirective) {
	return []string{"img"}, cobra.ShellCompDirectiveFilterFileExt
}

// ---------------------------------------------------------
// Initialization
// ---------------------------------------------------------

func init() {
	// 1. Attach to root command
	rootCmd.AddCommand(overlayCmd)

	// 2. Attach subcommands
	overlayCmd.AddCommand(overlayCreateCmd)
	overlayCmd.AddCommand(resizeCmd)
	overlayCmd.AddCommand(infoCmd)
	overlayCmd.AddCommand(checkCmd)
	overlayCmd.AddCommand(chownCmd)

	// 3. Define flags

	// --- Create ---
	overlayCreateCmd.Flags().StringP("size", "s", "10G", "Set overlay size (e.g., 500M, 10G)")
	overlayCreateCmd.Flags().StringP("type", "t", "balanced", "Overlay profile: small/balanced/large files")
	overlayCreateCmd.Flags().String("fs", "ext3", "Filesystem type: ext3 or ext4 (some system does not support ext4)")
	overlayCreateCmd.Flags().Bool("fakeroot", false, "Create a fakeroot-compatible overlay (owned by root)")
	overlayCreateCmd.Flags().Bool("sparse", false, "Create a sparse overlay image")
	overlayCreateCmd.Flags().StringP("file", "f", "", "Initialize with Conda environment file (.yml or .yaml)")

	// --- Flag completion helpers ---
	overlayCreateCmd.RegisterFlagCompletionFunc("fs", func(cmd *cobra.Command, args []string, toComplete string) ([]string, cobra.ShellCompDirective) {
		opts := []string{"ext3", "ext4"}
		res := make([]string, 0, len(opts))
		for _, o := range opts {
			if toComplete == "" || strings.HasPrefix(o, toComplete) {
				res = append(res, o)
			}
		}
		return res, cobra.ShellCompDirectiveNoFileComp
	})

	overlayCreateCmd.RegisterFlagCompletionFunc("type", func(cmd *cobra.Command, args []string, toComplete string) ([]string, cobra.ShellCompDirective) {
		opts := []string{"small", "balanced", "large"}
		res := make([]string, 0, len(opts))
		for _, o := range opts {
			if toComplete == "" || strings.HasPrefix(o, toComplete) {
				res = append(res, o)
			}
		}
		return res, cobra.ShellCompDirectiveNoFileComp
	})

	// Let the shell complete files for --file (yaml/yml)
	overlayCreateCmd.RegisterFlagCompletionFunc("file", func(cmd *cobra.Command, args []string, toComplete string) ([]string, cobra.ShellCompDirective) {
		return nil, cobra.ShellCompDirectiveFilterFileExt
	})

	// --- Resize ---
	resizeCmd.Flags().StringP("size", "s", "", "New size (e.g., 20G, 2048M)")
	_ = resizeCmd.MarkFlagRequired("size")

	// --- Check ---
	checkCmd.Flags().BoolP("force", "f", false, "Force check even if filesystem appears clean")

	// --- Chown ---
	uidDefault := fmt.Sprintf("current user %d", os.Getuid())
	gidDefault := fmt.Sprintf("current group %d", os.Getgid())

	chownCmd.Flags().IntP("uid", "u", -1, "User ID to set (default: "+uidDefault+")")
	chownCmd.Flags().IntP("gid", "g", -1, "Group ID to set (default: "+gidDefault+")")
	chownCmd.Flags().Bool("root", false, "Set UID and GID to 0 (root); will override -u and -g options")
	chownCmd.Flags().StringP("path", "p", "/ext3", "Path inside the overlay to change")
}

// initializeOverlayWithConda installs a conda environment from a YAML file into an overlay
// If envFile is empty, it creates a minimal environment with zlib
func initializeOverlayWithConda(overlayPath, envFile string, fakeroot bool) error {
	// Validate environment file if provided
	if envFile != "" {
		absEnvFile, err := filepath.Abs(envFile)
		if err != nil {
			return fmt.Errorf("failed to get absolute path of environment file: %w", err)
		}
		envFile = absEnvFile

		if !utils.FileExists(envFile) {
			return fmt.Errorf("environment file %s not found", envFile)
		}

		if !strings.HasSuffix(envFile, ".yml") && !strings.HasSuffix(envFile, ".yaml") {
			return fmt.Errorf("environment file must be .yml or .yaml")
		}
	}

	// Ensure base image exists
	if err := apptainer.EnsureBaseImage(); err != nil {
		return err
	}

	absOverlayPath, err := filepath.Abs(overlayPath)
	if err != nil {
		return fmt.Errorf("failed to get absolute path of overlay: %w", err)
	}

	var mmCreateCmd []string
	if envFile != "" {
		utils.PrintMessage("Initializing conda environment using %s...", utils.StylePath(envFile))
		mmCreateCmd = []string{"mm-create", "-f", envFile, "-y", "-q"}
	} else {
		utils.PrintMessage("Initializing minimal conda environment with small package (zlib)...")
		mmCreateCmd = []string{"mm-create", "zlib", "-y", "-q"}
	}

	// Run mm-create command
	opts := exec.Options{
		Overlays:    []string{absOverlayPath},
		Command:     mmCreateCmd,
		WritableImg: true,
		Fakeroot:    fakeroot,
		HideOutput:  true, // Suppress mm-create verbose output
	}

	if err := exec.Run(opts); err != nil {
		os.Remove(absOverlayPath)
		return fmt.Errorf("failed to run mm-create: %w", err)
	}

	// Clean up micromamba cache
	utils.PrintMessage("Cleaning up micromamba cache...")
	cleanOpts := exec.Options{
		Overlays:    []string{absOverlayPath},
		Command:     []string{"mm-clean", "-a", "-y", "-q"},
		WritableImg: true,
		Fakeroot:    fakeroot,
		HideOutput:  true, // Suppress mm-clean verbose output
	}

	if err := exec.Run(cleanOpts); err != nil {
		utils.PrintWarning("Failed to clean micromamba cache: %v", err)
		// Not a critical error, continue
	}

	utils.PrintSuccess("Conda env is created inside %s.", utils.StylePath(overlayPath))
	return nil
}
