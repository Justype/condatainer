package helper

import (
	"bufio"
	"fmt"
	"os"
	"strconv"
	"strings"
	"time"

	"github.com/Justype/condatainer/internal/utils"
)

// HelperParam represents a single #PARAM: declaration from a helper script.
// Format: #PARAM: KEY=default --long-flag,-s "Description"
// Use KEY=? to mark the param as optional (passes through empty; script handles it).
type HelperParam struct {
	Key       string `json:"key"`
	Default   string `json:"default,omitempty"`
	Optional  bool   `json:"optional,omitempty"` // true when declared as KEY=?
	LongFlag  string `json:"long_flag,omitempty"`
	ShortFlag string `json:"short_flag,omitempty"`
	Desc      string `json:"desc,omitempty"`
}

// HelperScriptMeta holds overlay and singleton metadata extracted from a helper script.
type HelperScriptMeta struct {
	// ImgPackages is a template of conda packages for guided overlay creation and
	// pre-submission package verification. {KEY} tokens are substituted from
	// resolved #PARAM: values before use.
	ImgPackages string
	// RequiredOverlays is a space-separated template of named condatainer overlays
	// (SquashFS images) to load before the helper runs. {KEY} tokens are substituted
	// from resolved #PARAM: values (e.g. "rstudio-server build-essential r{POSIT_R}").
	// Parsed from #REQUIRED_OVERLAYS: (alias: legacy #OVERLAY_PACKAGES:).
	RequiredOverlays string
	// PostInstallCmd is run inside the container after micromamba installs ImgPackages.
	PostInstallCmd string
	// Singleton is true when at most one running instance of this helper is allowed.
	Singleton bool
	// Binds is a list of extra bind-mount specs from #BIND: headers, in "src:dest" form.
	// {KEY} tokens are substituted from resolved #PARAM: values before use.
	// $VAR references (e.g. $CNT_HELPER_STATE_DIR, $HOME) are passed verbatim and
	// expanded by the shell when the wrapper runs on the compute node.
	Binds []string
	// ParamValues holds the valid-value lists parsed from #VALUE: headers.
	// Keys are param names (e.g. "POSIT_R", "CONDA_PYTHON"); values are the parsed list.
	ParamValues map[string][]string
	// NCPUs is the default CPU count from #NCPUS:. Zero means not set.
	NCPUs int
	// MemMB is the default memory in MB from #MEM:. Zero means not set.
	MemMB int64
	// Walltime is the default walltime from #TIME:. Zero means not set.
	Walltime time.Duration
	// GPU is the default GPU spec from #GPU:. Empty means not set.
	GPU string
}

// ParseValueSpec parses a #VALUE: value specification into a string slice.
// Formats:
//
//	"opt1 | opt2 | opt3"     → pipe-separated option list (spaces allowed in labels)
//	"1.0,2.0,3.0"            → comma-separated version list
//	"EXPR:range:1-41"        → expression (stored as single-element slice for future use)
func ParseValueSpec(spec string) []string {
	spec = strings.TrimSpace(spec)
	if spec == "" {
		return nil
	}
	if strings.HasPrefix(spec, "EXPR:") {
		return []string{spec}
	}
	if strings.Contains(spec, "|") {
		parts := strings.Split(spec, "|")
		result := make([]string, 0, len(parts))
		for _, p := range parts {
			if t := strings.TrimSpace(p); t != "" {
				result = append(result, t)
			}
		}
		return result
	}
	parts := strings.Split(spec, ",")
	result := make([]string, 0, len(parts))
	for _, p := range parts {
		if t := strings.TrimSpace(p); t != "" {
			result = append(result, t)
		}
	}
	return utils.SortVersionsDescending(result)
}

// ParseHelperParams extracts #PARAM: declarations from a helper script.
// Format per line: #PARAM: KEY=default --long-flag,-s "Description"
// - KEY=default: variable name and default value (empty after = means required)
// - --long-flag,-s: optional flag spec (comma-separated; either part may be omitted)
// - "Description": quoted description string
func ParseHelperParams(scriptPath string) ([]HelperParam, error) {
	if !utils.FileExists(scriptPath) {
		return nil, fmt.Errorf("helper script not found at %s", scriptPath)
	}
	file, err := os.Open(scriptPath)
	if err != nil {
		return nil, fmt.Errorf("failed to open script: %w", err)
	}
	defer file.Close()

	var params []HelperParam
	scanner := bufio.NewScanner(file)
	for scanner.Scan() {
		line := scanner.Text()
		if !strings.HasPrefix(line, "#PARAM:") {
			continue
		}
		rest := strings.TrimSpace(line[len("#PARAM:"):])
		if rest == "" {
			continue
		}
		p, err := parseHelperParamLine(rest)
		if err != nil {
			return nil, fmt.Errorf("invalid #PARAM: line %q: %w", line, err)
		}
		params = append(params, p)
	}
	if err := scanner.Err(); err != nil {
		return nil, fmt.Errorf("error reading script: %w", err)
	}
	return params, nil
}

// parseHelperParamLine parses the content after "#PARAM: " into a HelperParam.
// Expected: KEY=default [--long-flag[,-s]] "Description"
func parseHelperParamLine(s string) (HelperParam, error) {
	// Extract trailing quoted description first.
	desc, rest, err := extractQuotedSuffix(s)
	if err != nil {
		return HelperParam{}, fmt.Errorf("missing quoted description: %w", err)
	}

	// Split remaining into whitespace tokens.
	tokens := strings.Fields(rest)
	if len(tokens) == 0 {
		return HelperParam{}, fmt.Errorf("missing KEY=default")
	}

	// First token must be KEY=value (no spaces allowed).
	keyVal := tokens[0]
	eqIdx := strings.IndexByte(keyVal, '=')
	if eqIdx < 0 {
		return HelperParam{}, fmt.Errorf("KEY=default token missing '=': %q", keyVal)
	}
	key := keyVal[:eqIdx]
	def := keyVal[eqIdx+1:]
	optional := false
	if def == "?" {
		def = ""
		optional = true
	}
	if key == "" {
		return HelperParam{}, fmt.Errorf("empty KEY in: %q", keyVal)
	}

	// Remaining tokens are optional flag specs: "--long,-s" or "--long" or "-s"
	var longFlag, shortFlag string
	for _, tok := range tokens[1:] {
		parts := strings.SplitN(tok, ",", 2)
		for _, p := range parts {
			p = strings.TrimSpace(p)
			if strings.HasPrefix(p, "--") {
				longFlag = p
			} else if strings.HasPrefix(p, "-") && len(p) == 2 {
				shortFlag = p
			}
		}
	}

	return HelperParam{
		Key:       key,
		Default:   def,
		Optional:  optional,
		LongFlag:  longFlag,
		ShortFlag: shortFlag,
		Desc:      desc,
	}, nil
}

// extractQuotedSuffix finds and removes the last double-quoted string from s.
// Returns (quoted content, remainder, error).
func extractQuotedSuffix(s string) (string, string, error) {
	s = strings.TrimSpace(s)
	end := len(s) - 1
	if end < 1 || s[end] != '"' {
		return "", "", fmt.Errorf("no closing quote found")
	}
	start := strings.LastIndex(s[:end], "\"")
	if start < 0 {
		return "", "", fmt.Errorf("no opening quote found")
	}
	quoted := s[start+1 : end]
	rest := strings.TrimSpace(s[:start])
	return quoted, rest, nil
}

// ParseHelperScriptMeta extracts metadata from a helper script.
// Recognised headers: #NCPUS:, #MEM:, #TIME:, #GPU:, #IMG_PACKAGES:, #REQUIRED_OVERLAYS:,
// #OVERLAY_PACKAGES: (legacy), #POST_INSTALL_CMD:, #SINGLETON:, #BIND:, #VALUE:.
func ParseHelperScriptMeta(scriptPath string) (HelperScriptMeta, error) {
	if !utils.FileExists(scriptPath) {
		return HelperScriptMeta{}, fmt.Errorf("helper script not found at %s", scriptPath)
	}
	file, err := os.Open(scriptPath)
	if err != nil {
		return HelperScriptMeta{}, fmt.Errorf("failed to open script: %w", err)
	}
	defer file.Close()

	meta := HelperScriptMeta{ParamValues: make(map[string][]string)}
	scanner := bufio.NewScanner(file)
	for scanner.Scan() {
		line := scanner.Text()
		switch {
		case strings.HasPrefix(line, "#NCPUS:"):
			if n, err := strconv.Atoi(strings.TrimSpace(line[len("#NCPUS:"):])); err == nil && n > 0 {
				meta.NCPUs = n
			}
		case strings.HasPrefix(line, "#MEM:"):
			// Default unit is GB (bare number like "16" means 16 GB).
			rawMem := strings.TrimSpace(line[len("#MEM:"):])
			if mb, err := utils.ParseMemoryMBWithDefault(rawMem, "GB"); err == nil && mb > 0 {
				meta.MemMB = mb
			}
		case strings.HasPrefix(line, "#TIME:"):
			raw := strings.TrimSpace(line[len("#TIME:"):])
			// Bare integer defaults to hours (e.g. "12" → "12h").
			if _, err := strconv.Atoi(raw); err == nil {
				raw = raw + "h"
			}
			if d, err := utils.ParseWalltime(raw); err == nil && d > 0 {
				meta.Walltime = d
			}
		case strings.HasPrefix(line, "#GPU:"):
			meta.GPU = strings.TrimSpace(line[len("#GPU:"):])
		case strings.HasPrefix(line, "#IMG_PACKAGES:"):
			meta.ImgPackages = strings.TrimSpace(line[len("#IMG_PACKAGES:"):])
		case strings.HasPrefix(line, "#REQUIRED_OVERLAYS:"),
			strings.HasPrefix(line, "#OVERLAY_PACKAGES:"):
			// #OVERLAY_PACKAGES: is a legacy alias
			idx := strings.Index(line, ":") + 1
			meta.RequiredOverlays = strings.TrimSpace(line[idx:])
		case strings.HasPrefix(line, "#BIND:"):
			if b := strings.TrimSpace(line[len("#BIND:"):]); b != "" {
				meta.Binds = append(meta.Binds, b)
			}
		case strings.HasPrefix(line, "#POST_INSTALL_CMD:"):
			meta.PostInstallCmd = strings.TrimSpace(line[len("#POST_INSTALL_CMD:"):])
		case strings.HasPrefix(line, "#SINGLETON:"):
			v := strings.TrimSpace(line[len("#SINGLETON:"):])
			meta.Singleton = strings.EqualFold(v, "true") || v == "1"
		case strings.HasPrefix(line, "#VALUE:"):
			rest := strings.TrimSpace(line[len("#VALUE:"):])
			key, spec, ok := strings.Cut(rest, "=")
			if ok && strings.TrimSpace(key) != "" {
				meta.ParamValues[strings.TrimSpace(key)] = ParseValueSpec(spec)
			}
		}
	}
	if err := scanner.Err(); err != nil {
		return HelperScriptMeta{}, fmt.Errorf("error reading script: %w", err)
	}
	return meta, nil
}
