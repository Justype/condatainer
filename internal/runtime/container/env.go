package container

import (
	"bufio"
	"errors"
	"fmt"
	"os"
	"strings"

	"github.com/Justype/condatainer/catalog"
	"github.com/Justype/condatainer/internal/artifact/compare"
	"github.com/Justype/condatainer/internal/artifact/meta"
	"github.com/Justype/condatainer/internal/image/tool"
	"github.com/Justype/condatainer/internal/utils"
)

// EnvPrefix is where a writable .img's payload lives inside the container.
// An .img is a working image with no embedded metadata, so it has no recorded
// prefix; a frozen environment records meta.EnvPrefix and so collides with the
// .img it came from through the ordinary prefix check.
const EnvPrefix = meta.EnvPrefix

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
// runtime.json is the only source for a .sqf or .sif, and the only metadata this
// path reads: an adjacent .env sidecar is ignored, because an installed image is
// immutable and its metadata travels inside it, and the manifest is never opened
// here so that provenance can grow without costing every mount. A writable .img
// is the exception and reads its sidecar.
func resolveImage(cleanPath string) (Contribution, *Diagnostic) {
	if utils.IsImg(cleanPath) {
		return imgContribution(cleanPath)
	}

	rt, err := meta.ReadRuntime(cleanPath)
	if err != nil {
		return Contribution{}, degradeDiagnostic(cleanPath, err)
	}
	// The one comparison rule on the execution path, and only because
	// runtime.json is already in hand: a wrong-architecture image mounts cleanly
	// and fails somewhere further downstream, where the cause is unrecognizable.
	if err := compare.MountAllowed(rt); err != nil {
		return Contribution{}, &Diagnostic{
			Level: "warn",
			Message: fmt.Sprintf("%s was %v; mounted, but contributes no environment.",
				utils.StylePath(cleanPath), err),
		}
	}

	c := Contribution{
		Name:        rt.Name,
		Type:        rt.Type,
		Description: rt.Description,
		Prefix:      rt.Prefix,
		Configs:     map[string]string{},
		Notes:       map[string]string{},
	}
	for _, env := range rt.Env {
		c.Configs[env.Key] = env.Resolved(rt.Prefix)
		if env.Note != "" {
			c.Notes[env.Key] = strings.ReplaceAll(env.Note, "{prefix}", rt.Prefix)
		}
	}
	return c, nil
}

// degradeDiagnostic turns a failed runtime read into the message the user sees.
//
// The three cases are kept apart so a broken host never reads as a broken image:
// absent runtime metadata is expected for anything built before the format and
// for a plain Apptainer .sif, while a present-but-unreadable document, or a
// missing tool, is something the user can act on.
func degradeDiagnostic(path string, err error) *Diagnostic {
	name := utils.StylePath(path)
	switch {
	case errors.Is(err, meta.ErrNoRuntime):
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

	vars, err := ReadEnvSidecar(imgPath)
	if err != nil {
		return c, &Diagnostic{
			Level:   "warn",
			Message: fmt.Sprintf("Unable to read overlay env %s.env: %v", imgPath, err),
		}
	}
	for _, v := range vars {
		c.Configs[v.Key] = strings.ReplaceAll(v.Value, "{prefix}", EnvPrefix)
		if v.Note != "" {
			c.Notes[v.Key] = strings.ReplaceAll(v.Note, "{prefix}", EnvPrefix)
		}
	}
	return c, nil
}

// ReadEnvSidecar reads the .env beside a writable overlay, with {prefix} left
// intact — substituting it is the caller's job, because a frozen artifact keeps
// the token and resolves it when the image is loaded.
//
// A missing sidecar is not an error: an .img without one still gets its bin/ on
// PATH. This is exported because freeze has to carry these variables into the
// artifact: an installed image is immutable and its metadata travels inside it,
// so an environment that lost them at the freeze would lose them for good.
func ReadEnvSidecar(imgPath string) ([]meta.EnvVar, error) {
	file, err := os.Open(imgPath + ".env")
	if err != nil {
		if os.IsNotExist(err) {
			return nil, nil
		}
		return nil, err
	}
	defer file.Close()

	var vars []meta.EnvVar
	scanner := bufio.NewScanner(file)
	for scanner.Scan() {
		key, value, note, ok := parseEnvLine(scanner.Text())
		if !ok {
			continue
		}
		vars = append(vars, meta.EnvVar{Key: key, Value: value, Note: note})
	}
	if err := scanner.Err(); err != nil {
		return nil, err
	}
	return vars, nil
}

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
