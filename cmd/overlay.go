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
	Short: "Manage ext3 overlays (create, resize, check, info)",
	Long: `Utilities to manage ext3 overlays for Apptainer.

- Create a new overlay, or resize an existing one
- Check filesystem integrity, or inspect usage`,
}

// ---------------------------------------------------------
// 1. Create Command
// ---------------------------------------------------------

// overlayCreateHelp is the shared help body for 'overlay create' and its 'o'
// shortcut, so the two descriptions cannot drift apart.
const overlayCreateHelp = `Create an ext3 overlay, sized and tuned by profile (-p).

If no image path is given, defaults to 'env.img'.
Conda packages listed after -- initialize the environment inline.`

// overlayProfiles are the selectable filesystem tuning profiles for a new overlay.
var overlayProfiles = []string{"small", "balanced", "large"}

// registerOverlayCreateFlags registers the full flag set for creating an overlay,
// shared by 'overlay create' and its 'o' shortcut so the two cannot drift apart.
//
// --profile matches internal/overlay's Profile vocabulary; the old --type spelling
// was ambiguous with the overlay file type (.img vs .sqf) and is kept, deprecated,
// on -t so existing invocations keep working.
func registerOverlayCreateFlags(cmd *cobra.Command) {
	cmd.Flags().StringP("size", "s", "10G", "Set overlay size (e.g., 500M, 10g)")
	cmd.Flags().StringP("profile", "p", "balanced", "Overlay profile: small/balanced/large files")
	cmd.Flags().StringP("type", "t", "balanced", "Overlay profile (deprecated)")
	cmd.Flags().MarkDeprecated("type", "use --profile instead") //nolint:errcheck
	cmd.Flags().Bool("fakeroot", false, "Create a fakeroot-compatible overlay (owned by root)")
	cmd.Flags().BoolP("sparse", "S", false, "Create a sparse overlay image (no pre-allocation)")
	cmd.Flags().StringP("file", "f", "", "Initialize with Conda environment file (.yml or .yaml)")
	cmd.Flags().Bool("no-tmp", false, "Create directly at target path (slower on network filesystems)")

	cmd.RegisterFlagCompletionFunc("profile", func(cmd *cobra.Command, args []string, toComplete string) ([]string, cobra.ShellCompDirective) {
		res := make([]string, 0, len(overlayProfiles))
		for _, o := range overlayProfiles {
			if toComplete == "" || strings.HasPrefix(o, toComplete) {
				res = append(res, o)
			}
		}
		return res, cobra.ShellCompDirectiveNoFileComp
	}) //nolint:errcheck

	cmd.RegisterFlagCompletionFunc("file", func(cmd *cobra.Command, args []string, toComplete string) ([]string, cobra.ShellCompDirective) {
		return nil, cobra.ShellCompDirectiveFilterFileExt
	}) //nolint:errcheck
}

// overlayProfileFlag returns the chosen profile, honouring the deprecated --type alias.
func overlayProfileFlag(cmd *cobra.Command) string {
	if cmd.Flags().Changed("type") {
		v, _ := cmd.Flags().GetString("type")
		return v
	}
	v, _ := cmd.Flags().GetString("profile")
	return v
}

// runOverlayCreate is the shared implementation of 'overlay create' and its 'o'
// shortcut. It resolves the target path, parses the flags registered by
// registerOverlayCreateFlags, then builds the overlay either directly at the target
// or via a local tmp copy that is moved into place afterwards.
func runOverlayCreate(cmd *cobra.Command, args []string) {
	// 1. Split args at -- into image path and packages
	var packages []string
	if dashAt := cmd.ArgsLenAtDash(); dashAt >= 0 {
		packages = args[dashAt:]
		args = args[:dashAt]
	}

	// Handle Positional Argument (Image Path)
	path := "env.img"
	if len(args) > 0 {
		path = args[0]
	} else if len(args) > 1 {
		ExitWithError("Too many positional arguments before --.")
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
	profile := overlayProfileFlag(cmd)
	envFile, _ := cmd.Flags().GetString("file")
	noTmp, _ := cmd.Flags().GetBool("no-tmp")

	if envFile != "" && len(packages) > 0 {
		ExitWithError("Cannot use -f/--file and inline packages (--) at the same time.")
	}

	sizeMB, err := utils.ParseSizeToMB(sizeStr)
	if err != nil {
		ExitWithError("Invalid size format '%s': %v", sizeStr, err)
	}

	// Unknown profile names fall back to balanced, but warn so a typo is visible.
	resolvedProfile, err := overlay.ParseProfile(profile)
	if err != nil {
		utils.PrintWarning("Unknown profile %q, using 'balanced'. Valid: %s",
			profile, strings.Join(overlayProfiles, ", "))
		resolvedProfile = overlay.ProfileDefault
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
		Profile:        resolvedProfile,
		Sparse:         sparse,
		FilesystemType: "ext3",
	}

	condaIO := exec.IO{Stdin: os.Stdin, Stdout: os.Stdout, Stderr: os.Stderr}
	if noTmp {
		// 3a. Create directly at target path (no tmp), then conda init there.
		if err := overlay.CreateDirectly(cmd.Context(), opts); err != nil {
			if errors.Is(err, context.Canceled) || errors.Is(cmd.Context().Err(), context.Canceled) {
				utils.PrintWarning("Overlay creation cancelled.")
				return
			}
			ExitWithError("%v", err)
		}
		if err := initCondaInOverlay(cmd.Context(), path, envFile, packages, fakeroot, condaIO); err != nil {
			os.Remove(path)
			if errors.Is(err, context.Canceled) || errors.Is(cmd.Context().Err(), context.Canceled) {
				utils.PrintWarning("Overlay initialization cancelled.")
				return
			}
			ExitWithError("Failed to initialize overlay with conda environment: %v", err)
		}
	} else {
		// 3b. Create sparse at local tmp (fast I/O), conda init there, then move + allocate.
		tmpPath, err := overlay.CreateInTmp(cmd.Context(), opts)
		if err != nil {
			if errors.Is(err, context.Canceled) || errors.Is(cmd.Context().Err(), context.Canceled) {
				utils.PrintWarning("Overlay creation cancelled.")
				return
			}
			ExitWithError("%v", err)
		}

		if err := initCondaInOverlay(cmd.Context(), tmpPath, envFile, packages, fakeroot, condaIO); err != nil {
			os.Remove(tmpPath)
			utils.RemoveDirIfEmpty(filepath.Dir(tmpPath))
			if errors.Is(err, context.Canceled) || errors.Is(cmd.Context().Err(), context.Canceled) {
				utils.PrintWarning("Overlay initialization cancelled.")
				return
			}
			ExitWithError("Failed to initialize overlay with conda environment: %v", err)
		}

		utils.PrintMessage("Moving overlay to %s", utils.StylePath(path))
		copied, err := overlay.MoveOverlayCopied(cmd.Context(), tmpPath, path, sparse)
		if err != nil {
			os.Remove(tmpPath)
			utils.RemoveDirIfEmpty(filepath.Dir(tmpPath))
			ExitWithError("Failed to move overlay to destination: %v", err)
		}
		utils.RemoveDirIfEmpty(filepath.Dir(tmpPath))

		// Skip AllocateOverlay when io.Copy was used: zeros already written physically.
		if !sparse && !copied {
			overlay.AllocateOverlay(cmd.Context(), path, sizeMB)
		}
	}

	utils.PrintSuccess("Created overlay %s", utils.StylePath(path))
}

var overlayCreateCmd = &cobra.Command{
	Use:   "create [flags] [path] [-- packages...]",
	Short: "Create a new ext3 overlay",
	Long:  overlayCreateHelp,
	Example: `  condatainer overlay create # 10G with default inode ratio
  condatainer overlay create my_data.img -s 50g -p large
  condatainer overlay create --fakeroot --sparse
  condatainer overlay create -f environment.yml
  condatainer overlay create myenv.img -- python=3.11`,

	Args: cobra.ArbitraryArgs,

	Run: runOverlayCreate,
}

// ---------------------------------------------------------
// 2. Resize Command
// ---------------------------------------------------------

var resizeCmd = &cobra.Command{
	Use:   "resize [flags] <path>",
	Short: "Expand or shrink an overlay",
	Example: `  condatainer overlay resize env.img -s 20g
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
	Use:          "info [overlay]",
	Short:        "Show details about an overlay",
	Long:         overlayInfoHelp,
	Example:      overlayInfoExample,
	Args:         cobra.ExactArgs(1),
	SilenceUsage: true,

	// Enable Smart Tab Completion for both .img and .sqf files
	ValidArgsFunction: completeInfoArgs,

	RunE: runInfoOverlay,
}

// ---------------------------------------------------------
// 4. Check Command
// ---------------------------------------------------------

var checkCmd = &cobra.Command{
	Use:   "check [flags] <path>",
	Short: "Verify filesystem integrity (e2fsck)",
	Example: `  condatainer overlay check env.img      # Check and repair automatically
  condatainer overlay check env.img -f  # Force a check even if marked clean`,
	Args: cobra.ExactArgs(1),

	// Enable Smart Tab Completion for .img files
	ValidArgsFunction: completeImages,

	Run: func(cmd *cobra.Command, args []string) {
		path := args[0]
		force, _ := cmd.Flags().GetBool("force")

		if err := overlay.CheckIntegrity(cmd.Context(), path, force); err != nil {
			ExitWithError("%v", err)
		}
	},
}

// ---------------------------------------------------------
// 5. Chown Command
// ---------------------------------------------------------

var chownCmd = &cobra.Command{
	Use:   "chown [flags] <path>",
	Short: "Change file ownership inside an overlay",
	Long: `Walks the overlay and updates the UID/GID wihtout mounting.

Defaults to the current user and path '/' inside the image.`,
	Example: `  condatainer overlay chown env.img                   # Set entire image to current user
  condatainer overlay chown env.img --root            # Set entire image to root
  condatainer overlay chown env.img -u 1001 -g 1001   # Set to specific ID
  condatainer overlay chown env.img -p /ext3 -p /data # Chown multiple paths`,
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
// 6. Export Command
// ---------------------------------------------------------

var exportCmd = &cobra.Command{
	Use:   "export [flags] [overlay_path]",
	Short: "Export a conda environment from an overlay",
	Long: `Export a Conda environment from an overlay.

For .img overlays the environment prefix is /ext3/env.
For .sqf module overlays the environment prefix is /cnt/<name>/<version>.`,
	Example: `  condatainer overlay export env.img > environment.yml  # Save to a file
  condatainer overlay export env.img --from-history     # Only explicitly installed
  condatainer overlay export env.img -e                 # Explicit URLs (exact rebuild)`,
	Args: cobra.ExactArgs(1),
	RunE: runExportOverlay,
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
	overlayCmd.AddCommand(exportCmd)

	// 3. Define flags

	// --- Create ---
	registerOverlayCreateFlags(overlayCreateCmd)

	// --- Resize ---
	resizeCmd.Flags().StringP("size", "s", "", "New size (e.g., 20g, 2048M)")
	_ = resizeCmd.MarkFlagRequired("size")

	// --- Check ---
	checkCmd.Flags().BoolP("force", "f", false, "Force check even if filesystem appears clean")

	// --- Chown ---
	chownCmd.Flags().IntP("uid", "u", os.Getuid(), "User ID to set")
	chownCmd.Flags().IntP("gid", "g", os.Getgid(), "Group ID to set")
	chownCmd.Flags().Bool("root", false, "Set UID and GID to 0 (root); overrides -u and -g")
	chownCmd.Flags().StringArrayP("path", "p", []string{"/"}, "Path inside the overlay (repeatable)")

	// --- Export ---
	exportCmd.Flags().BoolP("explicit", "e", false, "Use explicit format")
	exportCmd.Flags().Bool("no-md5", false, "Disable md5")
	exportCmd.Flags().Bool("no-build", false, "Disable the build string in spec")
	exportCmd.Flags().Bool("no-builds", false, "Disable the build string in spec (alias)")
	exportCmd.Flags().Bool("channel-subdir", false, "Enable channel/subdir in spec")
	exportCmd.Flags().Bool("from-history", false, "Build environment spec from history")
	exportCmd.Flags().Bool("json", false, "Report all output as json")
}

// runExportOverlay runs mm-export (micromamba env export) inside a mounted overlay.
func runExportOverlay(cmd *cobra.Command, args []string) error {
	overlayPath := args[0]

	if !utils.FileExists(overlayPath) && !utils.DirExists(overlayPath) {
		return fmt.Errorf("overlay %s not found", overlayPath)
	}

	ResolveFlagAlias(cmd, "no-build", "no-builds")

	// Build micromamba export args
	flagsArgs := []string{}
	explicit, _ := cmd.Flags().GetBool("explicit")
	noMD5, _ := cmd.Flags().GetBool("no-md5")
	noBuild, _ := cmd.Flags().GetBool("no-build")
	channelSubdir, _ := cmd.Flags().GetBool("channel-subdir")
	fromHistory, _ := cmd.Flags().GetBool("from-history")
	jsonOut, _ := cmd.Flags().GetBool("json")

	if explicit {
		flagsArgs = append(flagsArgs, "-e")
	}
	if noMD5 {
		flagsArgs = append(flagsArgs, "--no-md5")
	}
	if noBuild {
		flagsArgs = append(flagsArgs, "--no-builds")
	}
	if channelSubdir {
		flagsArgs = append(flagsArgs, "--channel-subdir")
	}
	if fromHistory {
		flagsArgs = append(flagsArgs, "--from-history")
	}
	if jsonOut {
		flagsArgs = append(flagsArgs, "--json")
	}

	// Determine internal prefix to export from
	envPrefix := ""
	if utils.IsImg(overlayPath) {
		// For .img overlays export from /ext3/env (standard location)
		envPrefix = "/ext3/env"
		if !overlay.PathExists(overlayPath, "/ext3/env/conda-meta") {
			cmd.SilenceUsage = true
			return fmt.Errorf("overlay %s does not contain a conda environment at /ext3/env (missing conda-meta)", overlayPath)
		}
	} else if utils.IsSqf(overlayPath) {
		// For .sqf overlays: derive name/version from filename
		base := strings.TrimSuffix(filepath.Base(overlayPath), filepath.Ext(overlayPath))
		nv := utils.NormalizeNameVersion(base)
		if !overlay.PathExists(overlayPath, "/cnt/"+nv+"/conda-meta") {
			cmd.SilenceUsage = true
			return fmt.Errorf("overlay %s does not contain a conda environment at /cnt/%s (missing conda-meta)", overlayPath, nv)
		}
		envPrefix = "/cnt/" + nv
	} else {
		cmd.SilenceUsage = true
		return fmt.Errorf("unsupported overlay type: %s", overlayPath)
	}

	// Build final command: micromamba -p <prefix> env export [flags]
	cmdArgs := []string{"micromamba", "-p", envPrefix, "env", "export"}
	if len(flagsArgs) > 0 {
		cmdArgs = append(cmdArgs, flagsArgs...)
	}

	opts := exec.Options{
		Overlays:    []string{overlayPath},
		Command:     cmdArgs,
		WritableImg: false,
		EnvSettings: []string{},
		HidePrompt:  true,
	}

	if err := exec.Run(cmd.Context(), opts, exec.IO{Stdout: os.Stdout, Stderr: os.Stderr}); err != nil {
		return fmt.Errorf("failed to export environment: %w", err)
	}
	return nil
}

// initCondaInOverlay initializes a conda environment in an existing overlay.
// Priority: envFile (-f flag) > packages (-- args) > minimal (zlib).
// envFile and packages are mutually exclusive and validated before this call.
func initCondaInOverlay(ctx context.Context, overlayPath, envFile string, packages []string, fakeroot bool, io exec.IO) error {
	if envFile != "" {
		absEnvFile, err := filepath.Abs(envFile)
		if err != nil {
			return fmt.Errorf("failed to get absolute path of environment file: %w", err)
		}
		envFile = absEnvFile
		if !utils.FileExists(envFile) {
			return fmt.Errorf("environment file %s not found", envFile)
		}
		if !utils.IsYaml(envFile) {
			return fmt.Errorf("environment file must be .yml or .yaml")
		}
	}

	if err := ensureBaseImage(ctx); err != nil {
		return err
	}

	absOverlayPath, err := filepath.Abs(overlayPath)
	if err != nil {
		return fmt.Errorf("failed to get absolute path of overlay: %w", err)
	}

	if envFile != "" {
		utils.PrintMessage("Initializing conda environment using %s...", utils.StylePath(envFile))
	} else if len(packages) > 0 {
		utils.PrintMessage("Initializing conda environment with: %s...", strings.Join(packages, " "))
	} else {
		utils.PrintMessage("Initializing %s...", exec.DescribeInitialCondaPackages(nil))
	}

	pkgs := packages
	if envFile != "" {
		pkgs = []string{"-f", envFile}
	}
	if err := exec.InitCondaEnv(ctx, absOverlayPath, pkgs, fakeroot, io); err != nil {
		return fmt.Errorf("failed to run mm-create: %w", err)
	}

	utils.PrintSuccess("Conda env is created inside %s.", utils.StylePath(overlayPath))
	return nil
}
