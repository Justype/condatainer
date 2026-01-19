package cmd

import (
	"fmt"
	"os"
	"strconv"
	"strings"
	"unicode"

	"github.com/Justype/condatainer/internal/overlay"
	"github.com/Justype/condatainer/internal/utils"
	"github.com/spf13/cobra"
)

// ---------------------------------------------------------
// Parent Command
// ---------------------------------------------------------

var OverlayCmd = &cobra.Command{
	Use:   "overlay",
	Short: "Manage persistent overlay images (create, resize, check, info)",
	Long: `Utilities to manage ext3 overlay images for Apptainer.
Allows creating optimized filesystems for Conda, resizing existing images,
and verifying filesystem integrity without mounting.`,
}

// ---------------------------------------------------------
// 1. Create Command
// ---------------------------------------------------------

var overlayCreateCmd = &cobra.Command{
	Use:   "create [image_path]",
	Short: "Create a new sparse overlay image",
	Long: `Creates an ext3 overlay image optimized for specific workloads.
	
If no image path is provided, defaults to 'env.img'.`,
	Example: `  condatainer overlay create                   # Creates default env.img (10GB, conda type)
  condatainer overlay create my_data.img -s 50G -t data
  condatainer overlay create --fakeroot --sparse`,

	Args: cobra.RangeArgs(0, 1),

	Run: func(cmd *cobra.Command, args []string) {
		// 1. Handle Positional Argument (Image Path)
		path := "env.img"
		if len(args) > 0 {
			path = args[0]
		}

		if utils.FileExists(path) {
			utils.PrintError("File already exists: %s", path)
			os.Exit(1)
		}

		// 2. Parse Flags
		sizeStr, _ := cmd.Flags().GetString("size")
		fakeroot, _ := cmd.Flags().GetBool("fakeroot")
		sparse, _ := cmd.Flags().GetBool("sparse")
		typeFlag, _ := cmd.Flags().GetString("type")

		sizeMB, err := parseSizeToMB(sizeStr)
		if err != nil {
			utils.PrintError("Invalid size format '%s': %v", sizeStr, err)
			os.Exit(1)
		}

		// 3. Create the Overlay
		utils.PrintMessage("Creating overlay image '%s' (%d MB)...", utils.StylePath(path), sizeMB)

		if fakeroot {
			err = overlay.CreateForRoot(path, sizeMB, typeFlag, sparse)
		} else {
			err = overlay.CreateForCurrentUser(path, sizeMB, typeFlag, sparse)
		}

		if err != nil {
			utils.PrintError("%v", err)
			os.Exit(1)
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
		sizeMB, err := parseSizeToMB(sizeStr)
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
		inodePct := stats.InodeUsage()

		fmt.Println(utils.StyleTitle("Overlay Statistics"))
		fmt.Printf("  %-15s %s\n", "File:", utils.StylePath(path))
		fmt.Printf("  %-15s %s\n", "State:", stats.FilesystemState)
		fmt.Printf("  %-15s %s\n", "Last Mounted:", stats.LastMounted)
		fmt.Println(utils.StyleTitle("Usage"))
		fmt.Printf("  %-15s %s / %s (%.1f%%)\n", "Disk Space:",
			utils.StyleNumber(usedBytes),
			utils.StyleNumber(stats.TotalBlocks*stats.BlockSize),
			diskPct)
		fmt.Printf("  %-15s %d / %d (%.1f%%)\n", "Inodes:",
			stats.TotalInodes-stats.FreeInodes,
			stats.TotalInodes,
			inodePct)
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
// Helper: Size Parser
// ---------------------------------------------------------

// parseSizeToMB converts strings like "10G", "500M", "1024" into Megabytes.
func parseSizeToMB(input string) (int, error) {
	s := strings.TrimSpace(strings.ToUpper(input))
	if s == "" {
		return 0, fmt.Errorf("empty size")
	}

	// Determine multiplier
	multiplier := 1 // Default assumes MB if no suffix, or raw integer

	if strings.HasSuffix(s, "G") || strings.HasSuffix(s, "GB") {
		multiplier = 1024
		s = strings.TrimRight(s, "GB")
	} else if strings.HasSuffix(s, "M") || strings.HasSuffix(s, "MB") {
		multiplier = 1
		s = strings.TrimRight(s, "MB")
	} else if strings.HasSuffix(s, "K") || strings.HasSuffix(s, "KB") {
		return 0, fmt.Errorf("size too small (KB), minimum 1MB")
	}

	// Clean any remaining non-numeric chars just in case
	s = strings.TrimFunc(s, func(r rune) bool {
		return !unicode.IsDigit(r)
	})

	val, err := strconv.Atoi(s)
	if err != nil {
		return 0, fmt.Errorf("not a number")
	}

	return val * multiplier, nil
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
	// 1. Attach Parent to Root

	// 2. Attach Subcommands
	OverlayCmd.AddCommand(overlayCreateCmd)
	OverlayCmd.AddCommand(resizeCmd)
	OverlayCmd.AddCommand(infoCmd)
	OverlayCmd.AddCommand(checkCmd)
	OverlayCmd.AddCommand(chownCmd)

	// 3. Define Flags

	// --- Create ---
	overlayCreateCmd.Flags().StringP("size", "s", "10G", "Set overlay size (e.g., 500M, 10G)")
	overlayCreateCmd.Flags().StringP("type", "t", "conda", "Overlay type: 'conda' (small files), 'data' (large files)")
	overlayCreateCmd.Flags().Bool("fakeroot", false, "Create a fakeroot-compatible overlay (owned by root)")
	overlayCreateCmd.Flags().Bool("sparse", false, "Create a sparse overlay image")
	// overlayCreateCmd.Flags().StringP("file", "f", "", "Initialize with Conda environment file") (TODO)

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
