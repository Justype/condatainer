package container

import (
	"fmt"
	"os"
	"path/filepath"
	"sort"
	"strconv"
	"strings"

	"github.com/Justype/condatainer/internal/artifact/meta"
	"github.com/Justype/condatainer/internal/config"
	"github.com/Justype/condatainer/internal/image"
	"github.com/Justype/condatainer/internal/image/ext3"
	"github.com/Justype/condatainer/internal/utils"
)

var (
	commonEnvVars = []string{
		"LC_ALL=C.UTF-8",
		"LANG=C.UTF-8",
		"CURL_CA_BUNDLE=", // Unset CA bundle to avoid host ENV interference
		"SSL_CERT_FILE=",  // Unset SSL cert file for the same reason
		// Set FPATH to a standard set of zsh completion directories to ensure zsh completions work inside the container
		"FPATH=/usr/local/share/zsh/site-functions:/usr/share/zsh/vendor-completions:/usr/share/zsh/functions/Calendar:/usr/share/zsh/functions/Chpwd:/usr/share/zsh/functions/Completion:/usr/share/zsh/functions/Completion/AIX:/usr/share/zsh/functions/Completion/BSD:/usr/share/zsh/functions/Completion/Base:/usr/share/zsh/functions/Completion/Cygwin:/usr/share/zsh/functions/Completion/Darwin:/usr/share/zsh/functions/Completion/Debian:/usr/share/zsh/functions/Completion/Linux:/usr/share/zsh/functions/Completion/Mandriva:/usr/share/zsh/functions/Completion/Redhat:/usr/share/zsh/functions/Completion/Solaris:/usr/share/zsh/functions/Completion/Unix:/usr/share/zsh/functions/Completion/X:/usr/share/zsh/functions/Completion/Zsh:/usr/share/zsh/functions/Completion/openSUSE:/usr/share/zsh/functions/Exceptions:/usr/share/zsh/functions/MIME:/usr/share/zsh/functions/Math:/usr/share/zsh/functions/Misc:/usr/share/zsh/functions/Newuser:/usr/share/zsh/functions/Prompts:/usr/share/zsh/functions/TCP:/usr/share/zsh/functions/VCS_Info:/usr/share/zsh/functions/VCS_Info/Backends:/usr/share/zsh/functions/Zftp:/usr/share/zsh/functions/Zle",
	}
)

// SetupConfig holds all configuration needed to set up a container
type SetupConfig struct {
	Overlays       []string // Overlay paths (will be resolved)
	WritableImg    bool     // Whether .img overlays should be writable
	EnvSettings    []string // User-specified environment variables (KEY=VALUE format)
	BindPaths      []string // User-specified bind paths
	Fakeroot       bool     // Whether to use fakeroot
	ApptainerFlags []string // Additional apptainer flags to pass through
	GpuRequested   bool     // A script explicitly requires a GPU; forces detection past autoload_gpu:false
}

// SetupResult contains all the processed configuration ready for container execution
type SetupResult struct {
	Overlays       []string          // Resolved and ordered overlay paths
	OverlayArgs    []string          // Overlay paths with :ro/:rw suffixes
	EnvList        []string          // Complete environment variable list
	EnvNotes       map[string]string // Environment variable notes for display
	Diagnostics    []Diagnostic      // Non-fatal messages for callers to present or log
	BindPaths      []string          // Deduplicated bind paths
	Fakeroot       bool              // Final fakeroot setting (may be auto-enabled)
	ApptainerFlags []string          // Apptainer flags including GPU flags
	LastImg        string            // Path to .img overlay if present
}

// Diagnostic is a non-fatal setup message returned to presentation layers.
type Diagnostic struct {
	Level   string
	Message string
}

// Setup processes all container configuration and returns a ready-to-use result
func Setup(cfg SetupConfig) (*SetupResult, error) {
	// Resolve overlay paths
	overlays, err := ResolveOverlayPaths(cfg.Overlays)
	if err != nil {
		return nil, err
	}

	// Ensure at most one .img overlay
	if err := ensureSingleImage(overlays); err != nil {
		return nil, err
	}

	// A writable .img looks beside itself for a paired frozen snapshot before
	// anything else runs, so the snapshot participates in the collision check
	// and the ordering below like any overlay the caller listed explicitly.
	overlays, snapshotDiagnostics := autoloadSnapshot(overlays)

	// Refuse a mount where one payload would disappear under another
	if err := ensureDistinctPrefixes(overlays); err != nil {
		return nil, err
	}

	// Put .img overlay last, with a paired env snapshot immediately beneath it
	overlays = orderOverlays(overlays)

	// Process overlays and check availability
	overlayArgs := make([]string, 0, len(overlays))
	var lastImg string
	for _, ol := range overlays {
		isImg := utils.IsImg(ol)
		if isImg {
			lastImg = ol

			// Only lock .img files as requested
			if utils.FileExists(ol) && !utils.DirExists(ol) {
				// If it's the principal image and writableImg is true, we need an exclusive lock.
				// In orderOverlays, the principal image is always the last one.
				isPrincipalImg := (ol == overlays[len(overlays)-1])
				writeLock := isPrincipalImg && cfg.WritableImg

				if err := image.CheckAvailable(ol, writeLock); err != nil {
					return nil, err
				}
			}
		}

		overlayArgs = append(overlayArgs, FormatOverlayMount(ol, cfg.WritableImg))
	}

	// Build environment variables
	envList, envNotes, diagnostics := buildEnvironment(overlays, lastImg, cfg)
	diagnostics = append(snapshotDiagnostics, diagnostics...)

	// Build bind paths
	bindPaths := BindPaths()
	if len(cfg.BindPaths) > 0 {
		bindPaths = append(bindPaths, cfg.BindPaths...)
	}
	bindPaths = DeduplicateBindPaths(bindPaths)

	// Detect GPU flags
	apptainerFlags := append([]string{}, DetectGPUFlags(cfg.GpuRequested)...)
	apptainerFlags = append(apptainerFlags, cfg.ApptainerFlags...)

	return &SetupResult{
		Overlays:       overlays,
		OverlayArgs:    overlayArgs,
		EnvList:        envList,
		EnvNotes:       envNotes,
		Diagnostics:    diagnostics,
		BindPaths:      bindPaths,
		Fakeroot:       cfg.Fakeroot,
		ApptainerFlags: apptainerFlags,
		LastImg:        lastImg,
	}, nil
}

// buildEnvironment constructs the complete environment variable list
func buildEnvironment(overlays []string, lastImg string, cfg SetupConfig) ([]string, map[string]string, []Diagnostic) {
	// Collect overlay environment variables (from .env files)
	configs, notes, diagnostics := CollectOverlayEnv(overlays)
	envKeys := make([]string, 0, len(configs))
	for key := range configs {
		envKeys = append(envKeys, key)
	}
	sort.Strings(envKeys)

	// Start with overlay-specific vars
	envList := make([]string, 0, len(envKeys)+20)
	for _, key := range envKeys {
		envList = append(envList, fmt.Sprintf("%s=%s", key, configs[key]))
	}

	// Add layer tracking
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

	// Build PATH
	pathEnv := BuildPathEnv(overlays)

	// Add standard environment variables
	envList = append(envList,
		fmt.Sprintf("PATH=%s", pathEnv),
		fmt.Sprintf("PS1=%s \\[\\e[0;34m\\]\\w\\[\\e[0m\\]> ", ps1Prefix),
		fmt.Sprintf("IN_CONDATAINER=%d", layer),
	)
	envList = append(envList, commonEnvVars...)

	// Add .img-specific environment variables
	if lastImg != "" {
		if os.Getenv("IN_CONDATAINER") != "" {
			diagnostics = append(diagnostics, Diagnostic{
				Level:   "warn",
				Message: "You are trying to mount an .img overlay inside an existing CondaTainer environment. This may lead to unexpected behavior.",
			})
		}

		envList = append(envList,
			"CONDA_DEFINE_ENV=env",
			"RETICULATE_PYTHON=/cnt_env/bin/python",
		)

	}

	// Add user-specified environment variables (with validation)
	for _, setting := range cfg.EnvSettings {
		setting = strings.TrimSpace(setting)
		if setting == "" {
			continue
		}
		if !strings.Contains(setting, "=") {
			diagnostics = append(diagnostics, Diagnostic{
				Level:   "warn",
				Message: fmt.Sprintf("Invalid env setting %s. It should be in KEY=VALUE format. Skipping.", setting),
			})
			continue
		}
		envList = append(envList, setting)
	}

	// Runtime-owned markers go last so user-provided environment settings
	// cannot redirect management commands away from the mounted image.
	if lastImg != "" {
		writable := "0"
		if cfg.WritableImg {
			writable = "1"
		}
		if len(config.Global.Build.Channels) > 0 {
			envList = append(envList, "CNT_CONDA_CHANNELS="+strings.Join(config.Global.Build.Channels, "|"))
		}
		envList = append(envList,
			"CNT_CONDA_ROOT=/cnt_env",
			"CONDA_PREFIX=/cnt_env",
			"MAMBA_ROOT_PREFIX=/cnt_env",
			"CNT_CONDA_WRITABLE="+writable,
		)
	}

	// Prepare environment notes for display
	envNotes := make(map[string]string)
	for key, value := range configs {
		note := notes[key]
		if note == "" {
			note = value
		}
		envNotes[key] = note
	}
	if lastImg != "" {
		envNotes["CNT_CONDA_ROOT"] = "/cnt_env"
		if cfg.WritableImg {
			envNotes["CNT_CONDA_WRITABLE"] = "1"
		} else {
			envNotes["CNT_CONDA_WRITABLE"] = "0"
		}
	}

	return envList, envNotes, diagnostics
}

// AutoEnableFakeroot checks if fakeroot should be auto-enabled for writable .img overlays
// Returns the updated fakeroot setting
func AutoEnableFakeroot(lastImg string, writable bool, currentFakeroot bool) (bool, []Diagnostic) {
	if lastImg == "" || !writable || currentFakeroot {
		return currentFakeroot, nil
	}

	// Check UID status and auto-enable fakeroot if needed
	if status := ext3.InspectImageUIDStatus(lastImg); status == ext3.UIDStatusRoot {
		return true, []Diagnostic{{
			Level:   "note",
			Message: fmt.Sprintf("Root overlay %s detected. --fakeroot enabled automatically.", filepath.Base(lastImg)),
		}}
	} else if status == ext3.UIDStatusDifferentUser {
		return true, []Diagnostic{{
			Level:   "warn",
			Message: fmt.Sprintf("%s's inner UID differs from current user. --fakeroot enabled automatically.", filepath.Base(lastImg)),
		}}
	}

	return currentFakeroot, nil
}

// ensureSingleImage checks that at most one .img overlay is specified
func ensureSingleImage(overlays []string) error {
	imgCount := 0
	for _, overlay := range overlays {
		if utils.IsImg(overlay) {
			imgCount++
		}
	}
	if imgCount > 1 {
		return fmt.Errorf("only one .img overlay is allowed, found %d", imgCount)
	}
	return nil
}

// ensureDistinctPrefixes refuses a mount where two images claim one
// /cnt/<name> subtree.
//
// Overlays are disjoint subtrees, not stacked diffs: each payload lives under
// its own prefix and nothing merges them. Two images claiming one prefix
// therefore do not combine — the later mount wins the whole subtree and the
// earlier one silently contributes nothing, while both still appear on PATH and
// in the environment. That is a wrong container that looks like a working one,
// so it is an error rather than a warning.
//
// Two builds of one name are the case this catches: a project's restored copy
// and a flat install of the same name record the same prefix. A base, an OS
// image and anything without readable metadata record no prefix and are exempt
// for free. The same file named twice is redundant, not a collision. A writable
// .img is exempt too: it has no identity of its own, ensureSingleImage already
// guarantees at most one, and it is expected to sit on top of whatever
// env-typed .sqf is present — env.sqf + env.img is not a collision. Two
// env-typed .sqfs together still is: that is still two snapshots with no way to
// tell which one is meant.
func ensureDistinctPrefixes(overlays []string) error {
	return distinctPrefixes(overlays, func(path string) string {
		contribution, _ := resolveImage(path)
		return contribution.Prefix
	})
}

// distinctPrefixes is ensureDistinctPrefixes over a prefix lookup, so the rule
// can be exercised without a real image to read metadata out of.
func distinctPrefixes(overlays []string, prefixOf func(string) string) error {
	claimed := map[string]string{}
	for _, overlay := range overlays {
		path := cleanOverlayPath(overlay)
		if utils.IsImg(path) {
			continue
		}
		prefix := prefixOf(path)
		if prefix == "" {
			continue
		}
		held, taken := claimed[prefix]
		if !taken {
			claimed[prefix] = path
			continue
		}
		if held == path {
			continue
		}
		// Only an environment claims EnvPrefix — a frozen one records it — so
		// the collision there is the one the user already has a word for, and
		// prefixes are not it.
		if prefix == meta.EnvPrefix {
			return fmt.Errorf("%s and %s are both environment snapshots; mount one at a time",
				held, path)
		}
		return fmt.Errorf("%s and %s both install to %s; one would hide the other",
			held, path, prefix)
	}
	return nil
}

// autoloadSnapshot appends a writable .img's paired env-typed .sqf
// (LookupSnapshot) to overlays, unless one is already present in the list or
// the paired slot is Blocked by something that isn't a snapshot — autoload
// never guesses, it just leaves the .img to mount alone in that case.
func autoloadSnapshot(overlays []string) ([]string, []Diagnostic) {
	var imgPath string
	for _, overlay := range overlays {
		if path := cleanOverlayPath(overlay); utils.IsImg(path) {
			imgPath = path
			break
		}
	}
	if imgPath == "" {
		return overlays, nil
	}

	lookup := LookupSnapshot(imgPath)
	if lookup.Path == "" {
		return overlays, nil
	}
	for _, overlay := range overlays {
		if cleanOverlayPath(overlay) == lookup.Path {
			return overlays, nil // already given explicitly
		}
	}

	// lookup.Path is always in the same directory as imgPath (LookupSnapshot's
	// own invariant), so only its filename is shown — the full path would just
	// repeat imgPath's directory back. imgPath stays full: unlike a filename,
	// callers here (CLI, dashboard, job logs) have no shared notion of "cwd"
	// to shorten it against.
	diagnostic := Diagnostic{
		Level: "note",
		Message: fmt.Sprintf("autoloaded snapshot %s beside %s",
			utils.StylePath(filepath.Base(lookup.Path)), utils.StylePath(imgPath)),
	}
	return append(overlays, lookup.Path), []Diagnostic{diagnostic}
}

// orderOverlays puts a writable .img last and, immediately beneath it, the one
// env-typed .sqf present (autoloaded or explicit) — the newest delta on top of
// the snapshot it continues from. Every other overlay's relative order is
// untouched; this is one positioning rule for one specific artifact, not a
// general reordering of os/app/data.
func orderOverlays(overlays []string) []string {
	var img string
	var rest []string
	for _, overlay := range overlays {
		if utils.IsImg(cleanOverlayPath(overlay)) {
			img = overlay
		} else {
			rest = append(rest, overlay)
		}
	}
	if img == "" {
		return rest
	}

	var envSqf string
	var others []string
	for _, overlay := range rest {
		if envSqf == "" && isEnvSnapshotSqf(cleanOverlayPath(overlay)) {
			envSqf = overlay
			continue
		}
		others = append(others, overlay)
	}
	if envSqf != "" {
		others = append(others, envSqf)
	}
	return append(others, img)
}
