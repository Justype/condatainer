package container

import (
	"bufio"
	"fmt"
	"os"
	"strings"
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

// CollectOverlayEnv reads .env files from overlay paths and returns environment configs, notes,
// and non-fatal diagnostics for callers to present or log.
func CollectOverlayEnv(paths []string) (map[string]string, map[string]string, []Diagnostic) {
	configs := map[string]string{}
	notes := map[string]string{}
	var diagnostics []Diagnostic

	for _, overlay := range paths {
		if overlay == "" {
			continue
		}
		// Strip :ro/:rw suffix — ResolveOverlayPaths preserves it but the .env
		// file lives next to the bare path (e.g. dep.sqf.env, not dep.sqf:ro.env).
		cleanOverlay := strings.TrimSuffix(strings.TrimSuffix(overlay, ":ro"), ":rw")
		envPath := cleanOverlay + ".env"
		file, err := os.Open(envPath)
		if err != nil {
			if !os.IsNotExist(err) {
				diagnostics = append(diagnostics, Diagnostic{
					Level:   "warn",
					Message: fmt.Sprintf("Unable to read overlay env %s: %v", envPath, err),
				})
			}
			continue
		}

		scanner := bufio.NewScanner(file)
		for scanner.Scan() {
			line := strings.TrimSpace(scanner.Text())
			if line == "" {
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
			if _, exists := configs[key]; exists {
				diagnostics = append(diagnostics, Diagnostic{
					Level:   "info",
					Message: fmt.Sprintf("Environment variable %s is defined in multiple overlays. Using value from %s.", key, overlay),
				})
			}
			configs[key] = value
		}
		if err := scanner.Err(); err != nil {
			diagnostics = append(diagnostics, Diagnostic{
				Level:   "warn",
				Message: fmt.Sprintf("Failed to scan %s: %v", envPath, err),
			})
		}
		file.Close()
	}

	return configs, notes, diagnostics
}
