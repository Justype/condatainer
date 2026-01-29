package exec

import (
	"bufio"
	"os"
	"strings"

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

func collectOverlayEnv(paths []string) (map[string]string, map[string]string) {
	configs := map[string]string{}
	notes := map[string]string{}

	for _, overlay := range paths {
		if overlay == "" {
			continue
		}
		envPath := overlay + ".env"
		if !utils.FileExists(envPath) {
			continue
		}

		file, err := os.Open(envPath)
		if err != nil {
			utils.PrintWarning("Unable to read overlay env %s: %v", utils.StylePath(envPath), err)
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
				utils.PrintMessage("Environment variable %s is defined in multiple overlays. Using value from %s.",
					utils.StyleName(key), utils.StylePath(overlay))
			}
			configs[key] = value
		}
		if err := scanner.Err(); err != nil {
			utils.PrintWarning("Failed to scan %s: %v", utils.StylePath(envPath), err)
		}
		file.Close()
	}

	return configs, notes
}
