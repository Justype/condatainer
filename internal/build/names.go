package build

import (
	"strings"

	"github.com/Justype/condatainer/internal/config"
	"github.com/Justype/condatainer/internal/utils"
)

// ExpandDefaultDistroName normalizes a build-script overlay name and expands a
// bare name to <default_distro>/<name> only when that build script exists.
func ExpandDefaultDistroName(nameVersion string) (string, bool) {
	normalized := utils.NormalizeNameVersion(nameVersion)
	if !strings.Contains(normalized, "/") && !strings.Contains(normalized, "::") && config.Global.DefaultDistro != "" {
		candidate := config.Global.DefaultDistro + "/" + normalized
		if _, found := FindBuildScript(candidate); found {
			return candidate, true
		}
	}
	return normalized, false
}
