package overlay

import (
	"encoding/json"
	"os/exec"
	"strings"

	"github.com/Justype/condatainer/internal/utils"
)

// CondaInfo holds the channels and explicitly-requested packages parsed
// from a conda-meta/history file inside an overlay.
type CondaInfo struct {
	Channels []string
	Specs    []string
}

// ReadFile reads a file at innerPath from inside an overlay image.
// For .sqf overlays it uses unsquashfs -cat; for .img overlays it uses debugfs.
// Returns nil if the file cannot be read.
func ReadFile(overlayPath, innerPath string) []byte {
	if utils.IsSqf(overlayPath) {
		return Cat(overlayPath, innerPath)
	}
	if utils.IsImg(overlayPath) {
		return imgCat(overlayPath, innerPath)
	}
	return nil
}

// imgCat reads a file from inside an ext3 overlay image using debugfs.
// The overlay stores writable content under the "upper/" directory,
// so innerPath (e.g. "/ext3/env/conda-meta/history") is resolved as
// "upper/ext3/env/conda-meta/history" in the debugfs filesystem.
func imgCat(imgPath, innerPath string) []byte {
	dbg, err := exec.LookPath("debugfs")
	if err != nil {
		return nil
	}
	inner := strings.TrimPrefix(innerPath, "/")
	catArg := "cat upper/" + inner
	cmd := exec.Command(dbg, "-R", catArg, imgPath)
	out, err := cmd.Output()
	if err != nil {
		return nil
	}
	return out
}

// CondarcChannels reads the channels list (in priority order) from the .condarc
// at envPrefix inside the overlay. This is what conda reads as context.channels
// and orders `env export` by; micromamba ignores it and sorts alphabetically.
// Returns nil if the file is absent or lists no channels.
func CondarcChannels(overlayPath, envPrefix string) []string {
	data := ReadFile(overlayPath, strings.TrimSuffix(envPrefix, "/")+"/.condarc")
	if data == nil {
		return nil
	}
	return parseCondarcChannels(string(data))
}

// parseCondarcChannels extracts the block-list under a top-level `channels:` key
// from .condarc YAML, preserving order and stopping at the next top-level key.
func parseCondarcChannels(content string) []string {
	var channels []string
	inBlock := false
	for _, line := range strings.Split(content, "\n") {
		trimmed := strings.TrimSpace(line)
		if !inBlock {
			if trimmed == "channels:" {
				inBlock = true
			}
			continue
		}
		if strings.HasPrefix(trimmed, "- ") {
			if ch := strings.Trim(strings.TrimSpace(trimmed[2:]), `"'`); ch != "" {
				channels = append(channels, ch)
			}
			continue
		}
		if trimmed == "" {
			continue
		}
		break // a non-list, non-empty line ends the channels block
	}
	return channels
}

// ReadCondaInfo reads conda-meta/history from envPrefix inside the overlay
// and returns the channels and explicitly-installed package specs.
// Returns nil if the history file is absent or contains no specs.
func ReadCondaInfo(overlayPath, envPrefix string) *CondaInfo {
	histPath := strings.TrimSuffix(envPrefix, "/") + "/conda-meta/history"
	data := ReadFile(overlayPath, histPath)
	if data == nil {
		return nil
	}
	return parseCondaHistory(string(data))
}

// parseCondaHistory parses the text of a conda-meta/history file.
// Channels are extracted from installed package URLs (+https://conda.anaconda.org/<channel>/...).
// Explicitly-requested specs come from "# update specs:" JSON arrays, deduplicating
// by package name and keeping the last-seen spec for each name.
// Packages listed in "# remove specs:" are removed from the result.
func parseCondaHistory(content string) *CondaInfo {
	channelsSeen := map[string]bool{}
	var channels []string

	// specNames preserves insertion order; specMap holds the latest spec per name.
	var specNames []string
	specNamesSeen := map[string]bool{}
	specMap := map[string]string{}

	for _, line := range strings.Split(content, "\n") {
		line = strings.TrimSpace(line)
		const condaAnacondaOrg = "+https://conda.anaconda.org/"
		switch {
		case strings.HasPrefix(line, condaAnacondaOrg):
			// URL format: +https://conda.anaconda.org/<channel>/<subdir>::<pkg>
			rest := strings.TrimPrefix(line, condaAnacondaOrg)
			if ch, _, ok := strings.Cut(rest, "/"); ok {
				if !channelsSeen[ch] {
					channelsSeen[ch] = true
					channels = append(channels, ch)
				}
			}

		case strings.HasPrefix(line, "# update specs:"):
			after := strings.TrimPrefix(line, "# update specs:")
			var specs []string
			if err := json.Unmarshal([]byte(strings.TrimSpace(after)), &specs); err != nil {
				continue
			}
			for _, spec := range specs {
				// Base package name: everything before the first version operator.
				name := strings.FieldsFunc(spec, func(r rune) bool {
					return r == '=' || r == '>' || r == '<' || r == '!'
				})[0]
				if !specNamesSeen[name] {
					specNamesSeen[name] = true
					specNames = append(specNames, name)
				}
				specMap[name] = spec
			}

		case strings.HasPrefix(line, "# remove specs:"):
			after := strings.TrimPrefix(line, "# remove specs:")
			var specs []string
			if err := json.Unmarshal([]byte(strings.TrimSpace(after)), &specs); err != nil {
				continue
			}
			for _, spec := range specs {
				name := strings.FieldsFunc(spec, func(r rune) bool {
					return r == '=' || r == '>' || r == '<' || r == '!'
				})[0]
				delete(specMap, name)
			}
		}
	}

	if len(specNames) == 0 {
		return nil
	}

	var specs []string
	for _, name := range specNames {
		if s, ok := specMap[name]; ok {
			specs = append(specs, s)
		}
	}
	if len(specs) == 0 {
		return nil
	}
	return &CondaInfo{Channels: channels, Specs: specs}
}
