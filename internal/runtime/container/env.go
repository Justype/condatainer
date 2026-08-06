package container

import (
	"bufio"
	"bytes"
	"fmt"
	"os"
	"path/filepath"
	"strings"

	"github.com/Justype/condatainer/catalog"
	"github.com/Justype/condatainer/internal/image/squashfs"
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

// overlayPrefix derives an overlay's in-container mount root (/cnt/<name>/<version>)
// from its filename, e.g. images/orad--2.7.0.sqf -> /cnt/orad/2.7.0. This is the
// value {prefix} takes in the overlay's #ENV: declarations — the same location
// the recipe wrote to as $CNT_PREFIX, seen at load time.
func overlayPrefix(overlayPath string) string {
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
		prefix := overlayPrefix(cleanOverlay)

		// Embedded build-script env first, then let the sidecar shadow it. Both are
		// merged per-overlay so a sidecar overriding the embedded value is not
		// reported as a cross-overlay conflict below.
		_, ovConfigs, ovNotes := readEmbeddedEnv(cleanOverlay, prefix)
		_, diags := readSidecarEnv(cleanOverlay, prefix, ovConfigs, ovNotes)
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

// readEmbeddedEnv reads #DESCRIPTION:/#ENV: from the .cnt-build-script embedded in a
// .sqf overlay and returns the description, resolved env values, and notes, with
// {prefix} filled in. Non-.sqf overlays and overlays without an embedded script
// yield empty results. The script is read straight out of the archive via
// unsquashfs -cat, so no mount is needed.
func readEmbeddedEnv(overlayPath, prefix string) (description string, configs, notes map[string]string) {
	configs = map[string]string{}
	notes = map[string]string{}

	if !strings.HasSuffix(overlayPath, ".sqf") {
		return description, configs, notes
	}

	// The payload (and its embedded script) lives at cnt/<name>/<version>/ inside
	// the archive, mirroring the /cnt/<name>/<version> mount root.
	scriptPath := strings.TrimPrefix(prefix, "/") + "/" + utils.BuildScriptName
	data := squashfs.Cat(overlayPath, scriptPath)
	if len(data) == 0 {
		return description, configs, notes
	}

	// Parsed by the same reader the catalog uses, so a recipe cannot mean one
	// thing when it is resolved and another when its artifact is loaded.
	recipe, err := catalog.ParseRecipe(strings.TrimPrefix(prefix, "/cnt/"), bytes.NewReader(data))
	if err != nil {
		return description, configs, notes
	}
	vars := map[string]string{"prefix": prefix}
	description = recipe.Description
	for _, env := range recipe.Env {
		configs[env.Key] = env.Value(vars)
		if env.Note != "" {
			notes[env.Key] = env.Note
		}
	}

	return description, configs, notes
}

// readSidecarEnv overlays a sidecar <overlay>.env file on top of configs/notes,
// shadowing the embedded build-script env for local overrides. {prefix} is
// substituted with the mount root. Returns the sidecar's #DESCRIPTION: (empty if none) and
// non-fatal diagnostics. A missing sidecar is not an error.
func readSidecarEnv(cleanOverlay, prefix string, configs, notes map[string]string) (description string, diagnostics []Diagnostic) {
	envPath := cleanOverlay + ".env"
	file, err := os.Open(envPath)
	if err != nil {
		if !os.IsNotExist(err) {
			diagnostics = append(diagnostics, Diagnostic{
				Level:   "warn",
				Message: fmt.Sprintf("Unable to read overlay env %s: %v", envPath, err),
			})
		}
		return description, diagnostics
	}
	defer file.Close()

	scanner := bufio.NewScanner(file)
	for scanner.Scan() {
		line := strings.TrimSpace(scanner.Text())
		if line == "" {
			continue
		}
		if after, ok := strings.CutPrefix(line, "#DESCRIPTION:"); ok {
			description = strings.TrimSpace(after)
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
		configs[key] = strings.ReplaceAll(value, "{prefix}", prefix)
	}
	if err := scanner.Err(); err != nil {
		diagnostics = append(diagnostics, Diagnostic{
			Level:   "warn",
			Message: fmt.Sprintf("Failed to scan %s: %v", envPath, err),
		})
	}

	return description, diagnostics
}

// ResolveOverlayEnv resolves a single overlay's description (#DESCRIPTION:), env vars,
// and notes for display. The embedded build script is the source of truth (.sqf);
// a sidecar <overlay>.env shadows it. {prefix} is substituted with the overlay's
// mount root. A :ro/:rw suffix on the path is ignored.
func ResolveOverlayEnv(overlayPath string) (description string, configs, notes map[string]string) {
	cleanOverlay := strings.TrimSuffix(strings.TrimSuffix(overlayPath, ":ro"), ":rw")
	prefix := overlayPrefix(cleanOverlay)

	description, configs, notes = readEmbeddedEnv(cleanOverlay, prefix)
	sidecarDescription, _ := readSidecarEnv(cleanOverlay, prefix, configs, notes)
	if sidecarDescription != "" {
		description = sidecarDescription
	}
	return description, configs, notes
}
