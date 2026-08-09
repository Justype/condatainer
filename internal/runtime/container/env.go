package container

import (
	"bufio"
	"errors"
	"fmt"
	"os"
	"strings"

	"github.com/Justype/condatainer/catalog"
	"github.com/Justype/condatainer/internal/image/meta"
	"github.com/Justype/condatainer/internal/image/tool"
	"github.com/Justype/condatainer/internal/utils"
)

// EnvPrefix is where a writable .img's payload lives inside the container.
// An .img is a working image with no manifest, so it has no recorded prefix.
const EnvPrefix = "/cnt_env"

// Contribution is what one image adds to the container at load time.
//
// An image with no usable metadata yields the zero value and contributes
// nothing — no variables, no PATH entry, no description. It still mounts.
type Contribution struct {
	Name        string
	Type        catalog.Type
	Description string
	Prefix      string
	Configs     map[string]string
	Notes       map[string]string
}

// resolveImage reads what one image contributes, plus a diagnostic when it
// contributes nothing.
//
// The manifest is the only source for a .sqf or .sif; an adjacent .env sidecar
// is ignored for those, because an installed image is immutable and its metadata
// travels inside it. A writable .img is the exception and reads its sidecar.
func resolveImage(cleanPath string) (Contribution, *Diagnostic) {
	if utils.IsImg(cleanPath) {
		return imgContribution(cleanPath)
	}

	manifest, err := meta.Read(cleanPath)
	if err != nil {
		return Contribution{}, degradeDiagnostic(cleanPath, err)
	}

	c := Contribution{
		Name:        manifest.Name,
		Type:        manifest.Type,
		Description: manifest.Description,
		Prefix:      manifest.Runtime.Prefix,
		Configs:     map[string]string{},
		Notes:       map[string]string{},
	}
	for _, env := range manifest.Runtime.Env {
		c.Configs[env.Key] = env.Resolved(manifest.Runtime.Prefix)
		if env.Note != "" {
			c.Notes[env.Key] = strings.ReplaceAll(env.Note, "{prefix}", manifest.Runtime.Prefix)
		}
	}
	return c, nil
}

// degradeDiagnostic turns a failed manifest read into the message the user sees.
//
// The three cases are kept apart so a broken host never reads as a broken image:
// an absent manifest is expected for anything built before manifests existed and
// for a plain Apptainer .sif, while a present-but-unreadable one, or a missing
// tool, is something the user can act on.
func degradeDiagnostic(path string, err error) *Diagnostic {
	name := utils.StylePath(path)
	switch {
	case errors.Is(err, meta.ErrNoManifest):
		return &Diagnostic{
			Level:   "info",
			Message: fmt.Sprintf("%s has no CondaTainer metadata; mounted, but contributes no environment. Rebuild to add it.", name),
		}
	case errors.Is(err, meta.ErrInvalid), errors.Is(err, meta.ErrUnsupportedSchema):
		return &Diagnostic{
			Level:   "warn",
			Message: fmt.Sprintf("%s has unreadable CondaTainer metadata: %v. Mounted, but contributes no environment.", name, err),
		}
	case errors.Is(err, tool.ErrToolMissing):
		return &Diagnostic{
			Level:   "warn",
			Message: fmt.Sprintf("cannot read metadata from %s: %v. Install squashfs-tools to restore environment setup.", name, err),
		}
	default:
		return &Diagnostic{
			Level:   "warn",
			Message: fmt.Sprintf("cannot read metadata from %s: %v. Mounted, but contributes no environment.", name, err),
		}
	}
}

// imgContribution reads a writable .img's sidecar.
//
// An .img is mutable working state, so its environment lives beside it rather
// than inside it and can be edited without a rebuild.
func imgContribution(imgPath string) (Contribution, *Diagnostic) {
	c := Contribution{
		Type:    catalog.TypeApp,
		Prefix:  EnvPrefix,
		Configs: map[string]string{},
		Notes:   map[string]string{},
	}

	file, err := os.Open(imgPath + ".env")
	if err != nil {
		if os.IsNotExist(err) {
			return c, nil // no sidecar is normal; the .img still gets its bin/ on PATH
		}
		return c, &Diagnostic{
			Level:   "warn",
			Message: fmt.Sprintf("Unable to read overlay env %s.env: %v", imgPath, err),
		}
	}
	defer file.Close()

	scanner := bufio.NewScanner(file)
	for scanner.Scan() {
		key, value, note, ok := parseEnvLine(scanner.Text())
		if !ok {
			continue
		}
		c.Configs[key] = strings.ReplaceAll(value, "{prefix}", EnvPrefix)
		if note != "" {
			c.Notes[key] = strings.ReplaceAll(note, "{prefix}", EnvPrefix)
		}
	}
	if err := scanner.Err(); err != nil {
		return c, &Diagnostic{
			Level:   "warn",
			Message: fmt.Sprintf("Failed to scan %s.env: %v", imgPath, err),
		}
	}
	return c, nil
}

// parseEnvLine parses one `KEY=value ## note` sidecar line. Blank lines and
// lines starting with # are skipped, as is anything without a valid key.
func parseEnvLine(line string) (key, value, note string, ok bool) {
	line = strings.TrimSpace(line)
	if line == "" || strings.HasPrefix(line, "#") {
		return "", "", "", false
	}
	key, rest, found := strings.Cut(line, "=")
	key = strings.TrimSpace(key)
	if !found || key == "" {
		return "", "", "", false
	}
	// The value may itself contain '=', so only the first one separates.
	if v, n, hasNote := strings.Cut(rest, "##"); hasNote {
		return key, strings.TrimSpace(v), strings.TrimSpace(n), true
	}
	return key, strings.TrimSpace(rest), "", true
}

// CollectOverlayEnv resolves environment variables for the given overlay paths,
// returning configs, notes, and non-fatal diagnostics to present or log. It is
// the one place degradation is reported. See the README's Environment Variables.
func CollectOverlayEnv(paths []string) (map[string]string, map[string]string, []Diagnostic) {
	configs := map[string]string{}
	notes := map[string]string{}
	var diagnostics []Diagnostic

	for _, overlay := range paths {
		if overlay == "" {
			continue
		}
		contribution, diag := resolveImage(cleanOverlayPath(overlay))
		if diag != nil {
			diagnostics = append(diagnostics, *diag)
		}

		for key, value := range contribution.Configs {
			if _, exists := configs[key]; exists {
				diagnostics = append(diagnostics, Diagnostic{
					Level:   "info",
					Message: fmt.Sprintf("Environment variable %s is defined in multiple overlays. Using value from %s.", key, overlay),
				})
			}
			configs[key] = value
		}
		for key, note := range contribution.Notes {
			notes[key] = note
		}
	}

	return configs, notes, diagnostics
}

// cleanOverlayPath strips the :ro/:rw suffix ResolveOverlayPaths preserves; the
// image and its sidecar live at the bare path.
func cleanOverlayPath(overlay string) string {
	return strings.TrimSuffix(strings.TrimSuffix(overlay, ":ro"), ":rw")
}

// ResolveOverlayEnv resolves a single image's description, env vars, and notes
// for display. A :ro/:rw suffix on the path is ignored.
func ResolveOverlayEnv(overlayPath string) (description string, configs, notes map[string]string) {
	contribution, _ := resolveImage(cleanOverlayPath(overlayPath))
	if contribution.Configs == nil {
		contribution.Configs = map[string]string{}
	}
	if contribution.Notes == nil {
		contribution.Notes = map[string]string{}
	}
	return contribution.Description, contribution.Configs, contribution.Notes
}

// ResolveOverlayInfo returns an image's recorded identity for display, and
// whether it had usable metadata.
func ResolveOverlayInfo(overlayPath string) (Contribution, bool) {
	contribution, diag := resolveImage(cleanOverlayPath(overlayPath))
	return contribution, diag == nil
}
