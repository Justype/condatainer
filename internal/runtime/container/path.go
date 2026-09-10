package container

import (
	"strings"

	"github.com/Justype/condatainer/internal/utils"
)

// FormatOverlayMount formats an overlay path with the appropriate :ro or :rw suffix
func FormatOverlayMount(path string, writable bool) string {
	// Check if path already has :ro or :rw suffix
	if strings.HasSuffix(path, ":ro") || strings.HasSuffix(path, ":rw") {
		return path
	}

	if utils.IsImg(path) {
		// For .img files, add :ro or :rw suffix based on writable flag
		if writable {
			return path + ":rw"
		}
		return path + ":ro"
	} else if utils.IsSqf(path) {
		// For .sqf files, always add :ro suffix (they're always read-only)
		return path + ":ro"
	}
	// For .sif and other files, no suffix needed
	return path
}

// BuildPathEnv constructs the PATH environment variable from the overlays:
// anything ContributesBin() (an app, or a mounted conda environment) puts
// <prefix>/bin on PATH, every other type nothing, deduplicated by prefix —
// a writable .img and its paired env.sqf snapshot both claim EnvPrefix, but
// merge into the one physical directory at mount time, so it needs only one
// PATH entry. MPI_DIR's bin/ is prepended when set. See the README's
// Environment Variables.
func BuildPathEnv(overlays []string) string {
	// paths := []string{"/usr/sbin", "/usr/bin"}
	paths := []string{"$PATH"} // $PATH here is the PATH from the base image

	seen := map[string]bool{}
	for _, ov := range overlays {
		contribution, _ := resolveImage(cleanOverlayPath(ov))
		if !contribution.ContributesBin() || seen[contribution.Prefix] {
			continue
		}
		seen[contribution.Prefix] = true
		paths = append([]string{contribution.Prefix + "/bin"}, paths...)
	}

	// Last, so nothing in an image is shadowed by the bound executable.
	paths = append(paths, BoundExecDir)

	return strings.Join(paths, ":")
}
