package exec

import (
	"fmt"
	"os"
	"path/filepath"
	"sort"
	"strconv"
	"strings"

	"github.com/Justype/condatainer/internal/apptainer"
	"github.com/Justype/condatainer/internal/config"
	"github.com/Justype/condatainer/internal/overlay"
	"github.com/Justype/condatainer/internal/utils"
)

// Run executes a command inside a configured Apptainer container.
func Run(options Options) error {
	options = options.ensureDefaults()

	if err := apptainer.SetBin(options.ApptainerBin); err != nil {
		return err
	}

	overlays, err := ResolveOverlayPaths(options.Overlays)
	if err != nil {
		return err
	}

	if err := ensureSingleImage(overlays); err != nil {
		return err
	}
	overlays = putImgToLast(overlays)

	overlayArgs := make([]string, 0, len(overlays))
	var lastImg string
	for _, overlay := range overlays {
		if utils.IsImg(overlay) {
			lastImg = overlay
		}
		overlayArgs = append(overlayArgs, formatOverlayMount(overlay, options.WritableImg))
	}

	configs, notes := collectOverlayEnv(overlays)
	envKeys := make([]string, 0, len(configs))
	for key := range configs {
		envKeys = append(envKeys, key)
	}
	sort.Strings(envKeys)

	envList := make([]string, 0, len(envKeys)+10)
	for _, key := range envKeys {
		envList = append(envList, fmt.Sprintf("%s=%s", key, configs[key]))
	}

	layer := 0
	if raw := os.Getenv("IN_CONDATAINER"); raw != "" {
		if parsed, err := strconv.Atoi(raw); err == nil && parsed >= 1 {
			layer = parsed
		}
	}
	layer++
	ps1Prefix := "CNT"
	if layer > 1 {
		ps1Prefix = fmt.Sprintf("CNT_%d", layer)
	}

	pathEnv := buildPathEnv(overlays)

	// For nested usage: if this is a portable installation, prepend its bin directory to PATH
	// This allows nested condatainer calls to find the same installation
	if config.IsPortable() {
		if portableBin := config.GetPortableBinDir(); portableBin != "" {
			pathEnv = portableBin + ":" + pathEnv
		}
	}

	envList = append(envList,
		fmt.Sprintf("PATH=%s", pathEnv),
		"LC_ALL=C.UTF-8",
		"LANG=C.UTF-8",
		fmt.Sprintf("PS1=%s \\[\\e[0;34m\\]\\w\\[\\e[0m\\]> ", ps1Prefix),
		fmt.Sprintf("IN_CONDATAINER=%d", layer),
	)

	for _, setting := range options.EnvSettings {
		setting = strings.TrimSpace(setting)
		if setting == "" {
			continue
		}
		if !strings.Contains(setting, "=") {
			utils.PrintWarning("Invalid env setting %s. It should be in KEY=VALUE format. Skipping.", utils.StyleName(setting))
			continue
		}
		envList = append(envList, setting)
	}

	imgEnv := []string{}
	fakeroot := options.Fakeroot
	if lastImg != "" {
		if os.Getenv("IN_CONDATAINER") != "" {
			utils.PrintWarning("You are trying to mount an .img overlay inside an existing CondaTainer environment. This may lead to unexpected behavior.")
		}

		imgEnv = append(imgEnv,
			"MAMBA_ROOT_PREFIX=/ext3/env",
			"CONDA_PREFIX=/ext3/env",
			"CONDA_DEFINE_ENV=env",
			"RETICULATE_PYTHON=/ext3/env/bin/python",
		)

		if options.WritableImg {
			if !fakeroot {
				if status := overlay.InspectImageUIDStatus(lastImg); status == overlay.UIDStatusRoot {
					utils.PrintNote("Root overlay %s detected. --fakeroot enabled automatically.", utils.StylePath(filepath.Base(lastImg)))
					fakeroot = true
				} else if status == overlay.UIDStatusDifferentUser {
					utils.PrintWarning("%s's inner UID differs from current user. --fakeroot enabled automatically.", utils.StylePath(filepath.Base(lastImg)))
					fakeroot = true
				}
			}
			imgEnv = append(imgEnv, "CNT_CONDA_PREFIX=/ext3/env")
		}
	}
	envList = append(envList, imgEnv...)

	bindPaths := BindPaths() // BindPaths already collects all base directories internally
	if len(options.BindPaths) > 0 {
		bindPaths = append(bindPaths, options.BindPaths...)
	}

	if utils.DebugMode {
		utils.PrintDebug("[EXEC]Exec overlays: %v", overlays)
		utils.PrintDebug("[EXEC]Bind paths: %v", bindPaths)
		utils.PrintDebug("[EXEC]Overlay mounts: %v", overlayArgs)
		utils.PrintDebug("[EXEC]Env list: %v", envList)
		utils.PrintDebug("[EXEC]Command: %s", strings.Join(options.Command, " "))
		if len(options.EnvSettings) > 0 {
			utils.PrintDebug("[EXEC]Env overrides: %v", options.EnvSettings)
		}
	}

	envNotes := map[string]string{}
	for key, value := range configs {
		note := notes[key]
		if note == "" {
			note = value
		}
		envNotes[key] = note
	}
	if lastImg != "" && options.WritableImg {
		envNotes["CNT_CONDA_PREFIX"] = "/ext3/env"
	}

	if utils.IsInteractiveShell() && !options.HideOutput && !options.HidePrompt {
		if len(envNotes) > 0 {
			utils.PrintMessage("Overlay envs:")
			sortedNotes := make([]string, 0, len(envNotes))
			for key := range envNotes {
				sortedNotes = append(sortedNotes, key)
			}
			sort.Strings(sortedNotes)
			for _, key := range sortedNotes {
				fmt.Printf("  %s: %s\n", utils.StyleName(key), utils.StyleInfo(envNotes[key]))
			}
			fmt.Println("")
		} else if lastImg != "" && options.WritableImg {
			utils.PrintMessage("Overlay env:")
			fmt.Printf("  %s: %s\n", utils.StyleName("CNT_CONDA_PREFIX"), utils.StylePath("/ext3/env"))
			fmt.Println("")
		}
	}

	opts := &apptainer.ExecOptions{
		Bind:       bindPaths,
		Overlay:    overlayArgs,
		Env:        envList,
		Fakeroot:   fakeroot,
		HideOutput: options.HideOutput,
		Additional: []string{},
	}

	// Pass through stdin if requested (for interactive build scripts)
	if options.PassThruStdin {
		if options.Stdin != nil {
			opts.Stdin = options.Stdin
		} else {
			opts.Stdin = os.Stdin
		}
	}

	if err := apptainer.Exec(options.BaseImage, options.Command, opts); err != nil {
		return err
	}
	return nil
}

func formatOverlayMount(path string, writable bool) string {
	// Check if path already has :ro or :rw suffix
	if strings.HasSuffix(path, ":ro") || strings.HasSuffix(path, ":rw") {
		return path
	}

	if utils.IsImg(path) && writable {
		return path
	}
	return path + ":ro"
}

func buildPathEnv(overlays []string) string {
	paths := []string{"/usr/sbin", "/usr/bin"}

	for _, overlay := range overlays {
		// Strip :ro or :rw suffix for path checking
		cleanOverlay := strings.TrimSuffix(strings.TrimSuffix(overlay, ":ro"), ":rw")
		name := strings.TrimSuffix(filepath.Base(cleanOverlay), filepath.Ext(cleanOverlay))
		normalized := utils.NormalizeNameVersion(name)
		if normalized == "" {
			continue
		}
		if strings.Count(normalized, "/") > 1 {
			continue
		}
		var relative string
		if utils.IsImg(cleanOverlay) {
			relative = "/ext3/env/bin"
		} else if utils.IsSqf(cleanOverlay) {
			relative = fmt.Sprintf("/cnt/%s/bin", normalized)
		} else {
			utils.PrintWarning("Unknown overlay file extension for %s. Skipping PATH addition.", utils.StylePath(overlay))
			continue
		}
		paths = append([]string{relative}, paths...)
	}

	return strings.Join(paths, ":")
}
