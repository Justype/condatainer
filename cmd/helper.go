package cmd

import (
	"compress/gzip"
	"crypto/sha256"
	"encoding/json"
	"fmt"
	"io"
	"net/http"
	"os"
	"os/exec"
	"path/filepath"
	"strings"
	"time"

	"github.com/Justype/condatainer/internal/apptainer"
	"github.com/Justype/condatainer/internal/config"
	"github.com/Justype/condatainer/internal/scheduler"
	"github.com/Justype/condatainer/internal/utils"
	"github.com/spf13/cobra"
)

// ForceRefreshHelpers skips the helper metadata disk cache when set.
var ForceRefreshHelpers bool

// RemoteHelperMetadataCache is the on-disk envelope for cached helper script metadata.
type RemoteHelperMetadataCache struct {
	FetchedAt time.Time                            `json:"fetched_at"`
	SourceURL string                               `json:"source_url"`
	Metadata  map[string]map[string]rawHelperEntry `json:"metadata"`
}

// rawHelperEntry is the wire format for a helper script entry (no SourceURL).
type rawHelperEntry struct {
	Path string `json:"path"`
}

// in-process cache for helper metadata, keyed by base URL
var cachedHelperMetadataByURL = map[string]map[string]map[string]HelperScriptEntry{}

var (
	helperPath   bool
	helperList   bool
	helperUpdate bool
)

// supportedSchedulerTypes lists scheduler types supported for helper scripts.
// Add new supported scheduler types to this slice to enable additional helper categories.
var supportedSchedulerTypes = []scheduler.SchedulerType{
	scheduler.SchedulerSLURM,
}

var helperCmd = &cobra.Command{
	Use:   "helper [flags] [script-name] [script-args...]",
	Short: "Manage and run helper scripts",
	Long: `Manage helper scripts for running services inside CondaTainer on HPC.

Use --list to see available helper scripts.
Use --update to download or refresh helper scripts from remote.

Note: Helper is not available inside a container or a scheduler job.`,
	Example: `  condatainer helper --update                # Update all helper scripts
  condatainer helper --path                  # Show helper scripts directory
  condatainer helper code-server -p 11908    # Run code-server and forward flags`,
	RunE:              runHelper,
	ValidArgsFunction: completeHelperScripts,
}

func init() {
	rootCmd.AddCommand(helperCmd)
	helperCmd.Flags().BoolVar(&helperPath, "path", false, "Show path to helper scripts directory")
	helperCmd.Flags().BoolVarP(&helperList, "list", "l", false, "List available helper scripts")
	helperCmd.Flags().BoolVarP(&helperUpdate, "update", "u", false, "Update helper scripts from remote")

	// Stop flag parsing after the first positional argument so script flags (like -p or -v)
	// are passed through to the helper script rather than being interpreted by cobra.
	helperCmd.Flags().SetInterspersed(false)
}

// completeHelperScripts provides shell completion for helper script names
func completeHelperScripts(cmd *cobra.Command, args []string, toComplete string) ([]string, cobra.ShellCompDirective) {
	// Only complete the first argument (script name), then return to default file completion
	if len(args) > 0 {
		return nil, cobra.ShellCompDirectiveDefault
	}

	// Collect scripts from all search paths (user → system)
	seen := make(map[string]bool)
	var scripts []string

	for _, dir := range config.GetHelperScriptSearchPaths() {
		entries, err := os.ReadDir(dir)
		if err != nil {
			continue
		}

		for _, entry := range entries {
			name := entry.Name()
			// Skip hidden files and directories
			if !entry.IsDir() && !strings.HasPrefix(name, ".") {
				if !seen[name] {
					seen[name] = true
					scripts = append(scripts, name)
				}
			}
		}
	}

	return scripts, cobra.ShellCompDirectiveNoFileComp
}

func runHelper(cmd *cobra.Command, args []string) error {
	// Prevent helper commands when inside a container or scheduler job (except --path / --list)
	if config.IsInsideContainer() && !helperPath && !helperList {
		cmd.SilenceUsage = true
		ExitWithError("helper commands are not available inside a container")
	}
	if scheduler.IsInsideJob() && !helperPath && !helperList {
		cmd.SilenceUsage = true
		ExitWithError("helper commands are not available inside a scheduler job")
	}

	// --- Path Mode ---
	if helperPath {
		if len(args) > 0 {
			// Find specific script in all search paths
			scriptName := args[0]
			scriptPath, err := config.FindHelperScript(scriptName)
			if err != nil {
				cmd.SilenceUsage = true
				ExitWithError("helper script '%s' not found (searched: %v)", scriptName, config.GetHelperScriptSearchPaths())
			}
			fmt.Println(scriptPath)
			return nil
		}

		// Show all helper script search paths
		fmt.Println("Helper script search paths (priority order):")
		for i, dir := range config.GetHelperScriptSearchPaths() {
			exists := ""
			if !config.DirExists(dir) {
				exists = " (not found)"
			}
			fmt.Printf("  %d. %s%s\n", i+1, dir, exists)
		}

		// Show writable directory
		if writableDir, err := config.GetWritableHelperScriptsDir(); err == nil {
			fmt.Printf("\nWritable directory: %s\n", writableDir)
		}
		return nil
	}

	// --- List Mode ---
	if helperList {
		type scriptInfo struct {
			name   string
			whatis string
		}
		seen := make(map[string]bool)
		var scripts []scriptInfo
		maxNameLen := 0

		for _, dir := range config.GetHelperScriptSearchPaths() {
			entries, err := os.ReadDir(dir)
			if err != nil {
				continue
			}
			for _, entry := range entries {
				name := entry.Name()
				if !entry.IsDir() && !strings.HasPrefix(name, ".") && !seen[name] {
					seen[name] = true
					whatis := utils.GetWhatIsFromScript(filepath.Join(dir, name))
					scripts = append(scripts, scriptInfo{name, whatis})
					if len(name) > maxNameLen {
						maxNameLen = len(name)
					}
				}
			}
		}

		for _, s := range scripts {
			if s.whatis != "" {
				fmt.Printf("  %-*s  %s\n", maxNameLen, s.name, s.whatis)
			} else {
				fmt.Printf("  %s\n", s.name)
			}
		}
		return nil
	}

	// --- Update Mode ---
	if helperUpdate {
		// Get writable directory for updates
		helperScriptsDir, err := config.GetWritableHelperScriptsDir()
		if err != nil {
			cmd.SilenceUsage = true
			ExitWithError("failed to find writable helper scripts directory: %v", err)
		}
		if err := updateHelperScripts(args, helperScriptsDir); err != nil {
			cmd.SilenceUsage = true
			ExitWithError("%v", err)
		}
		return nil
	}

	// --- Run Mode ---
	if len(args) == 0 {
		// Parameter error - show usage
		return fmt.Errorf("no helper script name provided - use --update to update scripts")
	}

	scriptName := args[0]
	scriptArgs := args[1:]

	// Find script in all search paths
	scriptPath, err := config.FindHelperScript(scriptName)
	if err != nil {
		cmd.SilenceUsage = true
		ExitWithError("helper script '%s' not found\nSearched: %v\nRun '%s' to fetch available helper scripts",
			scriptName, config.GetHelperScriptSearchPaths(), utils.StyleAction("condatainer helper --update"))
	}

	// Ensure executable
	if err := os.Chmod(scriptPath, utils.PermExec); err != nil {
		utils.PrintDebug("Failed to chmod helper script: %v", err)
	}

	// Check apptainer is available before running helper script
	if err := apptainer.EnsureApptainer(); err != nil {
		cmd.SilenceUsage = true
		ExitWithError("apptainer is required to run helper scripts: %v", err)
	}

	// Run the helper script
	cmdArgs := append([]string{scriptPath}, scriptArgs...)
	utils.PrintDebug("Running helper command: %v", cmdArgs)

	helperCmd := exec.Command(cmdArgs[0], cmdArgs[1:]...)
	helperCmd.Stdin = os.Stdin
	helperCmd.Stdout = os.Stdout
	helperCmd.Stderr = os.Stderr

	if err := helperCmd.Run(); err != nil {
		// Propagate the exit code from the helper script
		if exitErr, ok := err.(*exec.ExitError); ok {
			os.Exit(exitErr.ExitCode())
		}
		return err
	}
	return nil
}

// updateHelperScripts downloads and updates helper scripts from remote metadata
func updateHelperScripts(args []string, helperScriptsDir string) error {
	// Fetch remote metadata
	metadata, err := GetAllRemoteHelperMetadata()
	if err != nil {
		return fmt.Errorf("failed to fetch remote helper metadata: %w", err)
	}

	// Choose category based on scheduler availability
	category := "headless"
	if config.Global.SubmitJob {
		// If user requested submission, ensure a supported scheduler is available
		sched, err := scheduler.DetectSchedulerWithBinary(config.Global.SchedulerBin)
		if err != nil {
			return fmt.Errorf("submit requested but no scheduler found: %v", err)
		}

		if !sched.IsAvailable() || sched.IsInsideJob() {
			return fmt.Errorf("scheduler is not available for submission (available=%v, in_job=%v)", sched.IsAvailable(), sched.IsInsideJob())
		}

		matched := false
		for _, st := range supportedSchedulerTypes {
			if sched.GetType() == st {
				category = strings.ToLower(string(st))
				matched = true
				break
			}
		}

		if !matched {
			return fmt.Errorf("unsupported scheduler type '%s'", sched.GetType())
		}
	}

	category_styled := utils.StyleName(category)
	utils.PrintDebug("Helper category selected: %s", category_styled)

	entries, ok := metadata[category]
	if !ok {
		return fmt.Errorf("no helper scripts found for category '%s'", category_styled)
	}

	// Create helper scripts directory
	if err := os.MkdirAll(helperScriptsDir, utils.PermDir); err != nil {
		return fmt.Errorf("failed to create helper scripts directory: %w", err)
	}

	// Filter to specific script if provided
	if len(args) > 0 {
		scriptName := args[0]
		if entry, ok := entries[scriptName]; ok {
			entries = map[string]HelperScriptEntry{scriptName: entry}
		} else {
			return fmt.Errorf("helper script '%s' not found in remote metadata for category '%s'", scriptName, category_styled)
		}
	} else {
		utils.PrintMessage("Updating all helper scripts for %s...", category_styled)
	}

	// Download each helper script
	for name, entry := range entries {
		if entry.Path == "" {
			continue
		}

		sourceURL := entry.SourceURL
		if sourceURL == "" {
			sourceURL = config.Global.ScriptsLink
		}
		url := fmt.Sprintf("%s/%s", sourceURL, entry.Path)
		destName := filepath.Base(entry.Path)
		dest := filepath.Join(helperScriptsDir, destName)

		utils.PrintMessage("Updating %s/%s", category, name)
		if err := downloadExecutable(url, dest); err != nil {
			utils.PrintWarning("Failed to update %s/%s: %v", category, name, err)
			continue
		}
	}

	utils.PrintSuccess("Helper update finished.")
	return nil
}

// HelperScriptEntry represents a helper script metadata entry.
type HelperScriptEntry struct {
	Path      string `json:"path"`
	SourceURL string `json:"-"` // injected during merge; not stored in remote JSON
}

// helperMetadataCachePathForURL returns the disk cache path for helper metadata from a given base URL.
func helperMetadataCachePathForURL(baseURL string) (string, error) {
	cacheDir, err := config.GetWritableCacheDir()
	if err != nil {
		return "", err
	}
	h := sha256.Sum256([]byte(baseURL))
	name := fmt.Sprintf("helper-scripts-%x.json.gz", h[:6])
	return filepath.Join(cacheDir, name), nil
}

func loadHelperMetadataCacheAny(path, baseURL string, ttl time.Duration, checkTTL bool) (*RemoteHelperMetadataCache, error) {
	var cache RemoteHelperMetadataCache
	if err := utils.ReadGzipJSONFile(path, &cache); err != nil {
		return nil, err
	}
	if cache.SourceURL != baseURL {
		return nil, fmt.Errorf("source URL changed")
	}
	if checkTTL && time.Since(cache.FetchedAt) > ttl {
		return nil, fmt.Errorf("cache expired")
	}
	return &cache, nil
}

// fetchHelperMetadataFromURL fetches (or loads from per-URL cache) helper metadata for one base URL.
// Injects SourceURL into each entry during conversion.
func fetchHelperMetadataFromURL(baseURL string) (map[string]map[string]HelperScriptEntry, error) {
	// In-process cache
	if cached, ok := cachedHelperMetadataByURL[baseURL]; ok {
		return cached, nil
	}

	metaURL := baseURL + "/metadata/helper-scripts.json.gz"
	cachePath, cacheErr := helperMetadataCachePathForURL(baseURL)

	ttl := config.Global.MetadataCacheTTL

	// Try disk cache
	if !ForceRefreshHelpers && cacheErr == nil && ttl > 0 {
		if cache, err := loadHelperMetadataCacheAny(cachePath, baseURL, ttl, true); err == nil {
			utils.PrintDebug("Using cached helper metadata for %s (fetched %s ago)", baseURL, time.Since(cache.FetchedAt).Round(time.Minute))
			result := injectHelperSourceURL(cache.Metadata, baseURL)
			cachedHelperMetadataByURL[baseURL] = result
			return result, nil
		}
	}

	// Network fetch
	utils.PrintDebug("Fetching helper metadata from %s...", baseURL)
	rawMeta, err := fetchRawHelperMetadata(metaURL)
	if err != nil {
		// Stale fallback
		if cacheErr == nil {
			if cache, staleErr := loadHelperMetadataCacheAny(cachePath, baseURL, 0, false); staleErr == nil {
				utils.PrintWarning("Network unavailable for %s; using cached helper metadata (fetched %s ago)",
					baseURL, time.Since(cache.FetchedAt).Round(time.Minute))
				result := injectHelperSourceURL(cache.Metadata, baseURL)
				cachedHelperMetadataByURL[baseURL] = result
				return result, nil
			}
		}
		return nil, err
	}

	// Persist to disk (non-fatal)
	if cacheErr == nil && ttl > 0 {
		envelope := RemoteHelperMetadataCache{
			FetchedAt: time.Now(),
			SourceURL: baseURL,
			Metadata:  rawMeta,
		}
		if err := utils.WriteGzipJSONFileAtomic(cachePath, envelope); err != nil {
			utils.PrintDebug("Failed to write helper metadata cache for %s: %v", baseURL, err)
		}
	}

	result := injectHelperSourceURL(rawMeta, baseURL)
	cachedHelperMetadataByURL[baseURL] = result
	return result, nil
}

// injectHelperSourceURL converts raw entries to HelperScriptEntry, injecting SourceURL.
func injectHelperSourceURL(raw map[string]map[string]rawHelperEntry, baseURL string) map[string]map[string]HelperScriptEntry {
	out := make(map[string]map[string]HelperScriptEntry, len(raw))
	for category, scripts := range raw {
		out[category] = make(map[string]HelperScriptEntry, len(scripts))
		for name, e := range scripts {
			out[category][name] = HelperScriptEntry{Path: e.Path, SourceURL: baseURL}
		}
	}
	return out
}

// fetchRawHelperMetadata performs the HTTP+gzip+JSON fetch for helper metadata.
func fetchRawHelperMetadata(url string) (map[string]map[string]rawHelperEntry, error) {
	resp, err := http.Get(url)
	if err != nil {
		return nil, err
	}
	defer resp.Body.Close()

	if resp.StatusCode != http.StatusOK {
		return nil, fmt.Errorf("failed to fetch metadata: HTTP %d", resp.StatusCode)
	}

	gzReader, err := gzip.NewReader(resp.Body)
	if err != nil {
		return nil, fmt.Errorf("failed to decompress metadata: %w", err)
	}
	defer gzReader.Close()

	data, err := io.ReadAll(gzReader)
	if err != nil {
		return nil, err
	}

	var metadata map[string]map[string]rawHelperEntry
	if err := json.Unmarshal(data, &metadata); err != nil {
		return nil, err
	}
	return metadata, nil
}

// GetAllRemoteHelperMetadata merges helper metadata from all configured remote URLs.
// Earlier URLs in ScriptsLinks take precedence (first match per category/name wins).
func GetAllRemoteHelperMetadata() (map[string]map[string]HelperScriptEntry, error) {
	if ForceRefreshHelpers {
		cachedHelperMetadataByURL = map[string]map[string]map[string]HelperScriptEntry{}
	}

	merged := map[string]map[string]HelperScriptEntry{}
	var firstErr error

	for _, baseURL := range config.Global.ScriptsLinks {
		meta, err := fetchHelperMetadataFromURL(baseURL)
		if err != nil {
			utils.PrintDebug("Failed to fetch helper metadata from %s: %v", baseURL, err)
			if firstErr == nil {
				firstErr = err
			}
			continue
		}
		for category, scripts := range meta {
			if _, ok := merged[category]; !ok {
				merged[category] = map[string]HelperScriptEntry{}
			}
			for name, entry := range scripts {
				if _, exists := merged[category][name]; !exists {
					merged[category][name] = entry // earlier URL wins
				}
			}
		}
	}

	if len(merged) == 0 && firstErr != nil {
		return nil, firstErr
	}
	return merged, nil
}

// downloadExecutable downloads a file and makes it executable
func downloadExecutable(url, destPath string) error {
	resp, err := http.Get(url)
	if err != nil {
		return err
	}
	defer resp.Body.Close()

	if resp.StatusCode != http.StatusOK {
		return fmt.Errorf("HTTP %d", resp.StatusCode)
	}

	// Create parent directory
	if err := os.MkdirAll(filepath.Dir(destPath), utils.PermDir); err != nil {
		return err
	}

	// Write file
	out, err := utils.CreateFileWritable(destPath)
	if err != nil {
		return err
	}
	defer out.Close()

	if _, err := io.Copy(out, resp.Body); err != nil {
		return err
	}

	// Make executable
	return os.Chmod(destPath, utils.PermExec)
}
