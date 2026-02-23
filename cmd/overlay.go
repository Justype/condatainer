package cmd

import (
	"context"
	"errors"
	"fmt"
	"os"
	"path/filepath"
	"strings"

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
				ExitWithError("Overlay image must have a .img extension.")
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
			ExitWithError("Path %s already exists.", utils.StylePath(path))
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
			ExitWithError("Invalid size format '%s': %v", sizeStr, err)
		}

		// 3. Create the Overlay
		if fakeroot {
			err = overlay.CreateForRoot(cmd.Context(), path, sizeMB, typeFlag, sparse, fsType, false)
		} else {
			err = overlay.CreateForCurrentUser(cmd.Context(), path, sizeMB, typeFlag, sparse, fsType, false)
		}

		if err != nil {
			if errors.Is(err, context.Canceled) || errors.Is(cmd.Context().Err(), context.Canceled) {
				utils.PrintWarning("Overlay creation cancelled.")
				return
			}
			ExitWithError("%v", err)
		}

		// 4. Initialize with Conda environment if file specified
		if envFile != "" {
			if err := initializeOverlayWithConda(cmd.Context(), path, envFile, fakeroot); err != nil {
				if errors.Is(err, context.Canceled) || errors.Is(cmd.Context().Err(), context.Canceled) {
					utils.PrintWarning("Overlay initialization cancelled.")
					return
				}
				ExitWithError("Failed to initialize overlay with conda environment: %v", err)
			}
		} else {
			// Initialize with minimal conda environment (zlib)
			if err := initializeOverlayWithConda(cmd.Context(), path, "", fakeroot); err != nil {
				if errors.Is(err, context.Canceled) || errors.Is(cmd.Context().Err(), context.Canceled) {
					utils.PrintWarning("Overlay initialization cancelled.")
					return
				}
				ExitWithError("Failed to initialize overlay with conda environment: %v", err)
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
			os.Exit(ExitCodeError)
		}

		// 2. Parse Size
		sizeMB, err := utils.ParseSizeToMB(sizeStr)
		if err != nil {
			ExitWithError("Invalid size format '%s': %v", sizeStr, err)
		}

		// 3. Execute
		absPath, _ := filepath.Abs(path)
		lock, err := overlay.AcquireLock(absPath, true)
		if err != nil {
			ExitWithError("%v", err)
		}
		defer lock.Close()

		err = overlay.Resize(cmd.Context(), path, sizeMB)
		if err != nil {
			ExitWithError("%v", err)
		}
	},
}

// ---------------------------------------------------------
// 3. Info Command
// ---------------------------------------------------------

var infoCmd = &cobra.Command{
	Use:   "info [path]",
	Short: "Display overlay info (disk/inode usage for .img, compression/inode stats for .sqf)",
	Long: `Display detailed information about an overlay file.

For ext3 (.img) overlays: shows filesystem stats, disk/inode usage, block size, and ownership.
For SquashFS (.sqf) overlays: shows compression method, level, inode count, block size, and mount path.

Delegates to the top-level 'condatainer info' command, so both commands produce identical output.`,
	Args:         cobra.ExactArgs(1),
	SilenceUsage: true,

	// Enable Smart Tab Completion for both .img and .sqf files
	ValidArgsFunction: completeOverlayAndSqfFiles,

	RunE: runInfoOverlay,
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

		absPath, _ := filepath.Abs(path)
		lock, err := overlay.AcquireLock(absPath, true)
		if err != nil {
			ExitWithError("%v", err)
		}
		defer lock.Close()

		err = overlay.CheckIntegrity(cmd.Context(), path, force)
		if err != nil {
			ExitWithError("%v", err)
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

Defaults to the current user and paths '/ext3' and '/opt' inside the image.
Use --root to force ownership to 0:0.

Multiple paths can be specified using multiple -p flags.`,
	Example: `  condatainer overlay chown env.img                    # Set /ext3 and /opt to current user
  condatainer overlay chown env.img --root -p /        # Set entire image to root
  condatainer overlay chown env.img -u 1001 -g 1001    # Set to specific ID
  condatainer overlay chown env.img -p /ext3           # Only chown /ext3
  condatainer overlay chown env.img -p /ext3 -p /data  # Chown multiple paths`,
	Args: cobra.ExactArgs(1),

	// Enable Smart Tab Completion for .img files
	ValidArgsFunction: completeImages,

	Run: func(cmd *cobra.Command, args []string) {
		path := args[0]

		// Parse Flags
		isRoot, _ := cmd.Flags().GetBool("root")
		targetUID, _ := cmd.Flags().GetInt("uid")
		targetGID, _ := cmd.Flags().GetInt("gid")
		internalPaths, _ := cmd.Flags().GetStringArray("path")

		// Root flag overrides uid/gid
		if isRoot {
			targetUID = 0
			targetGID = 0
			utils.PrintDebug("Mode: Root (0:0)")
		}

		// Feedback to user
		utils.PrintDebug("Chown Targets: %v inside %s", internalPaths, path)
		utils.PrintDebug("New Owner:     UID=%d GID=%d", targetUID, targetGID)

		absPath, _ := filepath.Abs(path)
		lock, err := overlay.AcquireLock(absPath, true)
		if err != nil {
			ExitWithError("%v", err)
		}
		defer lock.Close()

		// Perform the recursive chown for each path
		for _, internalPath := range internalPaths {
			err = overlay.ChownRecursively(cmd.Context(), path, targetUID, targetGID, internalPath)
			if err != nil {
				if errors.Is(err, context.Canceled) || errors.Is(cmd.Context().Err(), context.Canceled) {
					utils.PrintWarning("Operation cancelled.")
					return
				}
				utils.PrintDebug("Failed to chown path %s: %v", internalPath, err)
			}
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

// completeOverlayAndSqfFiles tells the shell to suggest both .img and .sqf files.
func completeOverlayAndSqfFiles(cmd *cobra.Command, args []string, toComplete string) ([]string, cobra.ShellCompDirective) {
	return []string{"img", "sqf"}, cobra.ShellCompDirectiveFilterFileExt
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
	overlayCreateCmd.Flags().String("fs", "ext3", "Filesystem type: ext3 or ext4 (Now only ext3 is supported)")
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
	chownCmd.Flags().IntP("uid", "u", os.Getuid(), "User ID to set")
	chownCmd.Flags().IntP("gid", "g", os.Getgid(), "Group ID to set")
	chownCmd.Flags().Bool("root", false, "Set UID and GID to 0 (root); overrides -u and -g")
	chownCmd.Flags().StringArrayP("path", "p", []string{"/ext3", "/opt"}, "Path inside the overlay (can specify multiple)")
}

// initializeOverlayWithConda installs a conda environment from a YAML file into an overlay
// If envFile is empty, it creates a minimal environment with zlib
func initializeOverlayWithConda(ctx context.Context, overlayPath, envFile string, fakeroot bool) error {
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
	if err := ensureBaseImage(ctx); err != nil {
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

	if err := exec.Run(ctx, opts); err != nil {
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

	if err := exec.Run(ctx, cleanOpts); err != nil {
		utils.PrintWarning("Failed to clean micromamba cache: %v", err)
		// Not a critical error, continue
	}

	utils.PrintSuccess("Conda env is created inside %s.", utils.StylePath(overlayPath))
	return nil
}
