package cmd

import (
	"context"
	"os"
	"path/filepath"
	"regexp"
	"sort"
	"strings"

	"github.com/Justype/condatainer/internal/build"
	"github.com/Justype/condatainer/internal/config"
	"github.com/Justype/condatainer/internal/container"
	"github.com/Justype/condatainer/internal/overlay"
	"github.com/Justype/condatainer/internal/utils"
	"github.com/spf13/cobra"
)

// ============================================================================
// Search / filter helpers shared by avail, list, and remove
// ============================================================================

// SearchMode describes how a SearchQuery matches candidate names.
type SearchMode int

const (
	SearchModeAnd     SearchMode = iota // all terms must substring-match (default)
	SearchModeExact                     // each term is an exact full-name match (OR across terms)
	SearchModePattern                   // single term compiled as regex (wildcard or literal regex)
)

// SearchQuery holds the compiled, normalised search state for one command invocation.
type SearchQuery struct {
	Raw     []string // normalised (lowercased) terms
	Mode    SearchMode
	pattern *regexp.Regexp // set only when Mode == SearchModePattern
}

// NewSearchQuery builds a SearchQuery from normalised filter terms.
// exactLookup(term) must return true when term is a case-insensitive exact full
// name in the candidate set (installed overlays or available build scripts).
// Pass nil to skip exact-first detection.
//
// Detection order (single term):
//  1. exactLookup match            → SearchModeExact
//  2. regex metacharacters present → SearchModePattern (regex)
//  3. '*' or '?' present           → SearchModePattern (wildcard anchored)
//  4. plain string                 → SearchModeAnd (substring)
//
// Multiple terms: exact-first detection only (no pattern matching).
func NewSearchQuery(filters []string, exactLookup func(string) bool) *SearchQuery {
	q := &SearchQuery{Raw: filters}
	if len(filters) == 0 {
		q.Mode = SearchModeAnd
		return q
	}

	// Multiple terms: exact-first only
	if len(filters) > 1 {
		if exactLookup != nil && exactLookup(filters[0]) {
			q.Mode = SearchModeExact
		} else {
			q.Mode = SearchModeAnd
		}
		return q
	}

	// Single term
	term := filters[0]

	// 1. Exact lookup takes highest priority
	if exactLookup != nil && exactLookup(term) {
		q.Mode = SearchModeExact
		return q
	}

	// 2. Regex metacharacters before wildcard check (e.g. ^term.*another)
	if strings.ContainsAny(term, `^$([+{|`) {
		re, err := regexp.Compile("(?i)" + term)
		if err != nil {
			utils.PrintWarning("Invalid regex pattern %q, using substring match: %v", term, err)
			q.Mode = SearchModeAnd
			return q
		}
		q.Mode = SearchModePattern
		q.pattern = re
		return q
	}

	// 3. Wildcard (* or ?)
	if strings.ContainsAny(term, "*?") {
		re, err := regexp.Compile("(?i)^" + wildcardToRegex(term) + "$")
		if err != nil {
			utils.PrintWarning("Invalid wildcard pattern %q, using substring match: %v", term, err)
			q.Mode = SearchModeAnd
			return q
		}
		q.Mode = SearchModePattern
		q.pattern = re
		return q
	}

	// 4. Plain substring
	q.Mode = SearchModeAnd
	return q
}

// wildcardToRegex converts a glob-style wildcard pattern to a regex fragment.
// '*' → '.*', '?' → '.', all other characters are escaped.
func wildcardToRegex(term string) string {
	var sb strings.Builder
	for _, ch := range term {
		switch ch {
		case '*':
			sb.WriteString(".*")
		case '?':
			sb.WriteByte('.')
		default:
			sb.WriteString(regexp.QuoteMeta(string(ch)))
		}
	}
	return sb.String()
}

// Matches reports whether name satisfies the query.
func (q *SearchQuery) Matches(name string) bool {
	if len(q.Raw) == 0 {
		return true
	}
	switch q.Mode {
	case SearchModePattern:
		return q.pattern.MatchString(name)
	case SearchModeExact:
		nameLower := strings.ToLower(name)
		for _, term := range q.Raw {
			if nameLower == term {
				return true
			}
		}
		return false
	default: // SearchModeAnd
		nameLower := strings.ToLower(name)
		for _, term := range q.Raw {
			if !strings.Contains(nameLower, term) {
				return false
			}
		}
		return true
	}
}

// MatchesOrAlias is like Matches but, in exact mode only, also accepts alias as an
// alternative name (e.g. "rstudio-server" for "ubuntu24/rstudio-server").
// Pass an empty alias to skip the alias check.
func (q *SearchQuery) MatchesOrAlias(name, alias string) bool {
	if q.Matches(name) {
		return true
	}
	if q.Mode == SearchModeExact && alias != "" {
		return q.Matches(alias)
	}
	return false
}

// HighlightRegexp returns a compiled regex suitable for single-pass term highlighting.
// Returns nil when there are no terms.
// In pattern mode the match regex itself is reused.
// In and/exact modes every raw term is highlighted literally.
func (q *SearchQuery) HighlightRegexp() *regexp.Regexp {
	if len(q.Raw) == 0 {
		return nil
	}
	if q.Mode == SearchModePattern {
		return q.pattern
	}
	sorted := make([]string, len(q.Raw))
	copy(sorted, q.Raw)
	sort.Slice(sorted, func(i, j int) bool { return len(sorted[i]) > len(sorted[j]) })
	parts := make([]string, len(sorted))
	for i, t := range sorted {
		parts[i] = regexp.QuoteMeta(t)
	}
	re, _ := regexp.Compile("(?i)(" + strings.Join(parts, "|") + ")")
	return re
}

// normalizeFilters lowercases and normalises a slice of raw user-supplied filter strings.
func normalizeFilters(filters []string) []string {
	normalized := make([]string, 0, len(filters))
	for _, filter := range filters {
		trimmed := strings.TrimSpace(filter)
		if trimmed == "" {
			continue
		}
		normalized = append(normalized, strings.ToLower(utils.NormalizeNameVersion(trimmed)))
	}
	return normalized
}

// CommonFlags holds the common flags used by exec, instance start, and instance exec
type CommonFlags struct {
	Overlays    []string
	WritableImg bool
	EnvSettings []string
	BaseImage   string
	BindPaths   []string
	Fakeroot    bool
}

// RegisterCommonFlags registers common flags on a cobra command
func RegisterCommonFlags(cmd *cobra.Command, flags *CommonFlags) {
	cmd.Flags().StringSliceVarP(&flags.Overlays, "overlay", "o", nil, "overlay file to mount (can be used multiple times)")
	cmd.Flags().BoolVarP(&flags.WritableImg, "writable", "w", false, "mount .img overlays as writable (default: read-only)")
	cmd.Flags().BoolVar(&flags.WritableImg, "writable-img", false, "Alias for --writable")
	cmd.Flags().StringSliceVar(&flags.EnvSettings, "env", nil, "set environment variable 'KEY=VALUE' (can be used multiple times)")
	cmd.Flags().StringVarP(&flags.BaseImage, "base-image", "b", "", "base image to use instead of default")
	cmd.Flags().StringSliceVar(&flags.BindPaths, "bind", nil, "bind path 'HOST:CONTAINER' (can be used multiple times)")
	cmd.Flags().BoolVarP(&flags.Fakeroot, "fakeroot", "f", false, "run container with fakeroot privileges")

	// Register completions
	cmd.RegisterFlagCompletionFunc("overlay", overlayFlagCompletion(true, true))
	cmd.RegisterFlagCompletionFunc("base-image", baseImageFlagCompletion())

	// Stop flag parsing after the first positional argument
	cmd.Flags().SetInterspersed(false)

	// Allow unknown flags to pass through (for apptainer)
	cmd.FParseErrWhitelist.UnknownFlags = true
}

// KnownFlags returns a map of all known condatainer flags
func KnownFlags() map[string]bool {
	return map[string]bool{
		"--overlay": true, "-o": true,
		"--writable": true, "--writable-img": true, "-w": true,
		"--env":        true,
		"--bind":       true,
		"--base-image": true, "-b": true,
		"--fakeroot": true, "-f": true,
		"--debug": true,
		"--local": true,
		"--quiet": true, "-q": true,
		"--yes": true, "-y": true,
	}
}

// Exit codes used by various commands
const (
	// Generic error code
	ExitCodeError = 1
	// Wrong or missing arguments (misuse of command)
	ExitCodeUsage = 2
	// Returned when jobs were submitted to a scheduler (overlays will be created asynchronously)
	ExitCodeJobsSubmitted = 3
)

// ExitIfJobsSubmitted exits the process with ExitCodeJobsSubmitted if a scheduler
// submission was requested and the graph shows job IDs were created.
func ExitIfJobsSubmitted(graph *build.BuildGraph) {
	if config.Global.SubmitJob && graph != nil && len(graph.GetJobIDs()) > 0 {
		utils.PrintNote("%d scheduler job(s) submitted. exiting with code %d",
			len(graph.GetJobIDs()), ExitCodeJobsSubmitted)
		os.Exit(ExitCodeJobsSubmitted)
	}
}

// ExitWithError prints an error and exits with ExitCodeError
func ExitWithError(format string, a ...interface{}) {
	utils.PrintError(format, a...)
	os.Exit(ExitCodeError)
}

// ExitWithUsageError prints an error and exits with ExitCodeUsage (2)
func ExitWithUsageError(format string, a ...interface{}) {
	utils.PrintError(format, a...)
	os.Exit(ExitCodeUsage)
}

// ParseCommandArgs parses arguments from os.Args after a given subcommand
// Returns: commands (positional args), apptainerFlags (unknown flags for apptainer)
func ParseCommandArgs(subcommand string) ([]string, []string) {
	var commands, apptainerFlags []string

	// Find where subcommand appears in os.Args
	cmdIdx := -1
	for i, arg := range os.Args {
		if arg == subcommand {
			cmdIdx = i
			break
		}
	}

	if cmdIdx == -1 {
		// Fallback: no args found
		return commands, apptainerFlags
	}

	knownFlags := KnownFlags()

	// Parse from after subcommand in os.Args
	commandStarted := false
	for i := cmdIdx + 1; i < len(os.Args); i++ {
		arg := os.Args[i]

		// Once we hit the command, everything after is part of the command
		if commandStarted {
			commands = append(commands, arg)
			continue
		}

		// Skip known flags (already handled by cobra)
		if knownFlags[arg] || isKnownFlagWithEquals(knownFlags, arg) {
			// Value flags (space-separated)
			if knownFlags[arg] && needsValue(arg) && i+1 < len(os.Args) {
				i++ // Skip value
			}
			continue
		}

		// Unknown flag → must use --flag=value format for apptainer pass-through
		if strings.HasPrefix(arg, "-") {
			apptainerFlags = append(apptainerFlags, arg)
			continue
		}

		// First non-flag argument starts the command
		commandStarted = true
		commands = append(commands, arg)
	}

	return commands, apptainerFlags
}

// ParseInstanceArgs parses arguments for instance start/run commands
// Returns: instanceName, apptainerFlags
func ParseInstanceArgs(subcommand string) (string, []string) {
	var instanceName string
	var apptainerFlags []string

	// Find where "instance <subcommand>" appears in os.Args
	cmdIdx := -1
	for i := 0; i < len(os.Args)-1; i++ {
		if os.Args[i] == "instance" && os.Args[i+1] == subcommand {
			cmdIdx = i + 1 // Point to the subcommand
			break
		}
	}

	if cmdIdx == -1 {
		// Fallback: no args found
		return instanceName, apptainerFlags
	}

	knownFlags := KnownFlags()

	// Parse from after subcommand in os.Args
	for i := cmdIdx + 1; i < len(os.Args); i++ {
		arg := os.Args[i]

		// Skip known flags (already handled by cobra)
		if knownFlags[arg] || isKnownFlagWithEquals(knownFlags, arg) {
			// Value flags (space-separated)
			if knownFlags[arg] && needsValue(arg) && i+1 < len(os.Args) {
				i++ // Skip value
			}
			continue
		}

		// Unknown flag → pass through to apptainer
		if strings.HasPrefix(arg, "-") {
			apptainerFlags = append(apptainerFlags, arg)
			continue
		}

		// First (and only) non-flag argument is the instance name
		if instanceName == "" {
			instanceName = arg
		}
		// Ignore any additional positional arguments - overlays must use -o flag
	}

	return instanceName, apptainerFlags
}

func isKnownFlagWithEquals(knownFlags map[string]bool, arg string) bool {
	if !strings.Contains(arg, "=") {
		return false
	}
	// Check if prefix matches a known flag
	parts := strings.SplitN(arg, "=", 2)
	return knownFlags[parts[0]]
}

func needsValue(flag string) bool {
	valueFlags := map[string]bool{
		"-o": true, "--overlay": true,
		"--env": true, "--bind": true,
		"-b": true, "--base-image": true,
	}
	return valueFlags[flag]
}

// ensureBaseImage ensures the base image exists
func ensureBaseImage(ctx context.Context) error {
	return build.EnsureBaseImage(ctx, false)
}

// ResolveBaseImage resolves a base image path to an absolute path
// Returns empty string if baseImage is empty
func ResolveBaseImage(baseImage string) string {
	if baseImage == "" {
		return ""
	}

	if utils.FileExists(baseImage) {
		return baseImage
	}

	resolvedBase, err := container.ResolveOverlayPaths([]string{baseImage})
	if err == nil && len(resolvedBase) > 0 {
		return resolvedBase[0]
	}

	return baseImage
}

// PrepareCommandAndHidePrompt prepares the command array and determines if prompt should be hidden
// If commands is empty, defaults to ["bash"] with hidePrompt=false
// If commands has 1 element, returns hidePrompt=false
// Otherwise returns hidePrompt=true
func PrepareCommandAndHidePrompt(commands []string) ([]string, bool) {
	hidePrompt := true
	if len(commands) == 0 {
		commands = []string{"bash"}
		hidePrompt = false
	} else if len(commands) == 1 {
		hidePrompt = false
	}
	return commands, hidePrompt
}

// ============================================================================
// Shell Completion Functions
// ============================================================================

// overlayFlagCompletion returns completion function for -o/--overlay flag
func overlayFlagCompletion(includeData bool, includeImg bool) func(*cobra.Command, []string, string) ([]string, cobra.ShellCompDirective) {
	return func(cmd *cobra.Command, args []string, toComplete string) ([]string, cobra.ShellCompDirective) {
		return overlaySuggestions(includeData, includeImg, toComplete)
	}
}

// baseImageFlagCompletion returns completion function for -b/--base-image flag
func baseImageFlagCompletion() func(*cobra.Command, []string, string) ([]string, cobra.ShellCompDirective) {
	return func(cmd *cobra.Command, args []string, toComplete string) ([]string, cobra.ShellCompDirective) {
		return systemOverlaySuggestions(toComplete)
	}
}

// systemOverlaySuggestions returns only OS overlays and local OS overlay files
func systemOverlaySuggestions(toComplete string) ([]string, cobra.ShellCompDirective) {
	installed, err := container.InstalledOverlays()
	if err != nil {
		return nil, cobra.ShellCompDirectiveError
	}

	choices := map[string]struct{}{}

	// Only include installed overlays that are OS overlays (contain .singularity.d)
	for name, path := range installed {
		if isOSOverlay(path) {
			if toComplete == "" || strings.HasPrefix(name, toComplete) {
				choices[name] = struct{}{}
			}
		}
	}

	addDistroAliasChoices(installed, choices, toComplete)

	// Add local .sif files and local .sqf files that are OS overlays - for -b flag
	for _, candidate := range localImageSuggestions(toComplete) {
		if utils.IsSif(candidate) {
			choices[candidate] = struct{}{}
			continue
		}
		absPath, err := filepath.Abs(candidate)
		if err == nil && isOSOverlay(absPath) {
			choices[candidate] = struct{}{}
		}
	}

	suggestions := make([]string, 0, len(choices))
	for choice := range choices {
		suggestions = append(suggestions, choice)
	}
	sort.Strings(suggestions)

	return suggestions, cobra.ShellCompDirectiveNoFileComp
}

// overlaySuggestions returns overlay suggestions including installed overlays and local files
func overlaySuggestions(includeData bool, includeImg bool, toComplete string) ([]string, cobra.ShellCompDirective) {
	installed, err := container.InstalledOverlays()
	if err != nil {
		return nil, cobra.ShellCompDirectiveError
	}

	choices := map[string]struct{}{}
	for name, path := range installed {
		if includeData || isAppOverlay(path) {
			if toComplete == "" || strings.HasPrefix(name, toComplete) {
				choices[name] = struct{}{}
			}
		}
	}

	addDistroAliasChoices(installed, choices, toComplete)

	for _, candidate := range localOverlaySuggestions(toComplete, includeImg) {
		choices[candidate] = struct{}{}
	}

	suggestions := make([]string, 0, len(choices))
	for choice := range choices {
		suggestions = append(suggestions, choice)
	}
	sort.Strings(suggestions)

	// Add space after completions (all suggestions are files now, no folders)
	return suggestions, cobra.ShellCompDirectiveNoFileComp
}

// findLocalFilesWithFilter is a shared helper for finding local files recursively
// Includes directories for navigation and recursively finds files up to maxDepth
func findLocalFilesWithFilter(toComplete string, maxDepth int, fileFilter func(name string) bool) []string {
	pathDir, _ := filepath.Split(toComplete)
	dirForRead := pathDir
	if dirForRead == "" {
		dirForRead = "."
	}

	suggestions := []string{}

	// Recursively find files up to maxDepth
	var findFiles func(dir string, prefix string, currentDepth int)
	findFiles = func(dir string, prefix string, currentDepth int) {
		entries, err := os.ReadDir(dir)
		if err != nil {
			return
		}

		for _, entry := range entries {
			name := entry.Name()

			// Skip hidden files/directories (starting with .)
			if strings.HasPrefix(name, ".") {
				continue
			}

			candidate := name
			if prefix != "" {
				candidate = prefix + name
			}

			if entry.IsDir() {
				// Recurse into subdirectories if within depth limit (don't add directories to suggestions)
				if currentDepth < maxDepth {
					findFiles(filepath.Join(dir, name), candidate+"/", currentDepth+1)
				}
			} else {
				// Apply file filter
				if fileFilter(name) {
					if toComplete == "" || strings.HasPrefix(candidate, toComplete) {
						suggestions = append(suggestions, candidate)
					}
				}
			}
		}
	}

	findFiles(dirForRead, pathDir, 0)
	return suggestions
}

// localImageSuggestions returns local .sqf and .sif files (for -b flag)
// Includes directories for navigation and recursively finds files up to 1 level deep
func localImageSuggestions(toComplete string) []string {
	return findLocalFilesWithFilter(toComplete, 1, func(name string) bool {
		// Include .sqf and .sif files only (not .img)
		return utils.IsSqf(name) || utils.IsSif(name)
	})
}

// localOverlaySuggestions returns local overlay files (for -o flag)
// Includes directories for navigation and recursively finds files up to 1 level deep
func localOverlaySuggestions(toComplete string, includeImg bool) []string {
	return findLocalFilesWithFilter(toComplete, 1, func(name string) bool {
		// Filter based on file type
		if !utils.IsOverlay(name) {
			return false
		}

		// Filter .sif and .img based on includeImg parameter
		if utils.IsSif(name) && includeImg {
			// Skip .sif files when includeImg is true (for -o flag)
			return false
		}
		if utils.IsImg(name) && !includeImg {
			// Skip .img files when includeImg is false (for -b flag)
			return false
		}

		return true
	})
}

// addDistroAliasChoices adds shorthand aliases for OS overlays matching the default distro.
// For each installed OS overlay named "<distro>/<name>", also suggests "<name>".
func addDistroAliasChoices(installed map[string]string, choices map[string]struct{}, toComplete string) {
	distro := config.Global.DefaultDistro
	if distro == "" {
		return
	}
	prefix := distro + "/"
	for name, path := range installed {
		if !strings.HasPrefix(name, prefix) || !isOSOverlay(path) {
			continue
		}
		alias := strings.TrimPrefix(name, prefix)
		if toComplete == "" || strings.HasPrefix(alias, toComplete) {
			choices[alias] = struct{}{}
		}
	}
}

// isOSOverlay reports whether a SquashFS overlay is an OS overlay.
// Delegates to overlay.IsOSType which checks for .singularity.d in the archive.
func isOSOverlay(overlayPath string) bool {
	return overlay.IsOSType(overlayPath)
}

// isAppOverlay checks if an overlay path is considered an "app" overlay
func isAppOverlay(path string) bool {
	// OS overlays (Apptainer-built SquashFS containing .singularity.d) are
	// considered "app" overlays regardless of their location
	if isOSOverlay(path) {
		return true
	}

	// Check if path is in any of the image search directories
	// Regular overlays (created via mksquashfs) are only considered "app" if
	// they're in the images directory
	for _, imagesDir := range config.GetImageSearchPaths() {
		if pathWithinDir(path, imagesDir) {
			return true
		}
	}
	return false
}

// pathWithinDir checks if a path is within a directory
func pathWithinDir(path, dir string) bool {
	if dir == "" {
		return false
	}
	rel, err := filepath.Rel(dir, path)
	if err != nil {
		return false
	}
	if rel == "." {
		return true
	}
	return !(rel == ".." || strings.HasPrefix(rel, ".."+string(os.PathSeparator)))
}
