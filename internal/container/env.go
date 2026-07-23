package container

import (
	"bufio"
	"fmt"
	"os"
	"path/filepath"
	"strings"

	"github.com/Justype/condatainer/internal/overlay"
	"github.com/Justype/condatainer/internal/utils"
)

func splitKeyNote(content string) (string, string) {
	content = strings.TrimSpace(content)
	if content == "" {
		return "", ""
	}
	if idx := strings.Index(content, "="); idx >= 0 {
		key := strings.TrimSpace(content[:idx])
		note := strings.TrimSpace(content[idx+1:])
		return key, note
	}
	return content, ""
}

// overlayAppRoot derives an overlay's in-container mount root (/cnt/<name>/<version>)
// from its filename, e.g. images/orad--2.7.0.sqf -> /cnt/orad/2.7.0. This is the
// value that replaces $app_root in the overlay's #ENV: declarations.
func overlayAppRoot(overlayPath string) string {
	base := filepath.Base(overlayPath)
	base = strings.TrimSuffix(strings.TrimSuffix(base, ".sqf"), ".img")
	return "/cnt/" + strings.ReplaceAll(base, "--", "/")
}

// CollectOverlayEnv resolves environment variables for the given overlay paths and
// returns environment configs, notes, and non-fatal diagnostics for callers to
// present or log.
//
// For each overlay the embedded build script (/cnt/<name>/<version>/.cnt-build-script)
// is the source of truth, keeping the .sqf self-contained; a sidecar <overlay>.env,
// if present, shadows it for local overrides. $app_root is substituted with the
// overlay's mount root at load time. When the same variable is set by more than one
// overlay, the later overlay wins and a diagnostic is recorded.
func CollectOverlayEnv(paths []string) (map[string]string, map[string]string, []Diagnostic) {
	configs := map[string]string{}
	notes := map[string]string{}
	var diagnostics []Diagnostic

	for _, overlay := range paths {
		if overlay == "" {
			continue
		}
		// Strip :ro/:rw suffix — ResolveOverlayPaths preserves it but the embedded
		// script and .env file live at the bare path (e.g. dep.sqf, not dep.sqf:ro).
		cleanOverlay := strings.TrimSuffix(strings.TrimSuffix(overlay, ":ro"), ":rw")
		appRoot := overlayAppRoot(cleanOverlay)

		// Embedded build-script env first, then let the sidecar shadow it. Both are
		// merged per-overlay so a sidecar overriding the embedded value is not
		// reported as a cross-overlay conflict below.
		_, ovConfigs, ovNotes := readEmbeddedEnv(cleanOverlay, appRoot)
		_, diags := readSidecarEnv(cleanOverlay, appRoot, ovConfigs, ovNotes)
		diagnostics = append(diagnostics, diags...)

		for key, value := range ovConfigs {
			if _, exists := configs[key]; exists {
				diagnostics = append(diagnostics, Diagnostic{
					Level:   "info",
					Message: fmt.Sprintf("Environment variable %s is defined in multiple overlays. Using value from %s.", key, overlay),
				})
			}
			configs[key] = value
		}
		for key, note := range ovNotes {
			notes[key] = note
		}
	}

	return configs, notes, diagnostics
}

// readEmbeddedEnv parses #WHATIS:/#ENV:/#ENVNOTE: from the .cnt-build-script
// embedded in a .sqf overlay and returns the description, resolved env values, and
// notes, with $app_root replaced by appRoot. Non-.sqf overlays and overlays without
// an embedded script yield empty results. The build script is read directly from
// the archive via unsquashfs -cat, so no mount is needed.
func readEmbeddedEnv(overlayPath, appRoot string) (whatis string, configs, notes map[string]string) {
	configs = map[string]string{}
	notes = map[string]string{}

	if !strings.HasSuffix(overlayPath, ".sqf") {
		return whatis, configs, notes
	}

	// The payload (and its embedded script) lives at cnt/<name>/<version>/ inside
	// the archive, mirroring the /cnt/<name>/<version> mount root.
	scriptPath := strings.TrimPrefix(appRoot, "/") + "/" + utils.BuildScriptName
	data := overlay.Cat(overlayPath, scriptPath)
	if len(data) == 0 {
		return whatis, configs, notes
	}

	lines := strings.Split(string(data), "\n")
	for i, line := range lines {
		line = strings.TrimSpace(line)
		if after, ok := strings.CutPrefix(line, "#WHATIS:"); ok {
			whatis = strings.TrimSpace(after)
			continue
		}
		if !strings.HasPrefix(line, "#ENV:") {
			continue
		}
		content := utils.StripInlineComment(line[len("#ENV:"):])
		if !strings.Contains(content, "=") {
			continue
		}
		parts := strings.SplitN(content, "=", 2)
		key := strings.TrimSpace(parts[0])
		value := strings.TrimSpace(parts[1])
		if key == "" {
			continue
		}
		configs[key] = strings.ReplaceAll(value, "$app_root", appRoot)

		// An #ENVNOTE: on the following line describes the variable just parsed.
		if i+1 < len(lines) {
			next := strings.TrimSpace(lines[i+1])
			if strings.HasPrefix(next, "#ENVNOTE:") {
				notes[key] = strings.TrimSpace(next[len("#ENVNOTE:"):])
			}
		}
	}

	return whatis, configs, notes
}

// readSidecarEnv overlays a sidecar <overlay>.env file on top of configs/notes,
// shadowing the embedded build-script env for local overrides. $app_root is
// substituted with appRoot. Returns the sidecar's #WHATIS: (empty if none) and
// non-fatal diagnostics. A missing sidecar is not an error.
func readSidecarEnv(cleanOverlay, appRoot string, configs, notes map[string]string) (whatis string, diagnostics []Diagnostic) {
	envPath := cleanOverlay + ".env"
	file, err := os.Open(envPath)
	if err != nil {
		if !os.IsNotExist(err) {
			diagnostics = append(diagnostics, Diagnostic{
				Level:   "warn",
				Message: fmt.Sprintf("Unable to read overlay env %s: %v", envPath, err),
			})
		}
		return whatis, diagnostics
	}
	defer file.Close()

	scanner := bufio.NewScanner(file)
	for scanner.Scan() {
		line := strings.TrimSpace(scanner.Text())
		if line == "" {
			continue
		}
		if after, ok := strings.CutPrefix(line, "#WHATIS:"); ok {
			whatis = strings.TrimSpace(after)
			continue
		}
		if strings.HasPrefix(line, "#ENVNOTE:") {
			key, note := splitKeyNote(line[len("#ENVNOTE:"):])
			if key != "" {
				notes[key] = note
			}
			continue
		}
		if strings.HasPrefix(line, "#") {
			continue
		}
		pair := strings.SplitN(line, "=", 2)
		if len(pair) != 2 {
			continue
		}
		key := strings.TrimSpace(pair[0])
		value := strings.TrimSpace(pair[1])
		if key == "" {
			continue
		}
		configs[key] = strings.ReplaceAll(value, "$app_root", appRoot)
	}
	if err := scanner.Err(); err != nil {
		diagnostics = append(diagnostics, Diagnostic{
			Level:   "warn",
			Message: fmt.Sprintf("Failed to scan %s: %v", envPath, err),
		})
	}

	return whatis, diagnostics
}

// ResolveOverlayEnv resolves a single overlay's description (#WHATIS:), env vars,
// and notes for display. The embedded build script is the source of truth (.sqf);
// a sidecar <overlay>.env shadows it. $app_root is substituted with the overlay's
// mount root. A :ro/:rw suffix on the path is ignored.
func ResolveOverlayEnv(overlayPath string) (whatis string, configs, notes map[string]string) {
	cleanOverlay := strings.TrimSuffix(strings.TrimSuffix(overlayPath, ":ro"), ":rw")
	appRoot := overlayAppRoot(cleanOverlay)

	whatis, configs, notes = readEmbeddedEnv(cleanOverlay, appRoot)
	sidecarWhatis, _ := readSidecarEnv(cleanOverlay, appRoot, configs, notes)
	if sidecarWhatis != "" {
		whatis = sidecarWhatis
	}
	return whatis, configs, notes
}
