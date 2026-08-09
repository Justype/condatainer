package build

import (
	"fmt"

	"github.com/Justype/condatainer/internal/runtime/container"
)

// dependencyOverlays resolves each build dependency to the overlay providing it,
// as read-only overlay specs ready to mount while the recipe runs.
func dependencyOverlays(dependencies []string) ([]string, error) {
	if len(dependencies) == 0 {
		return nil, nil
	}

	depPaths, err := container.ResolveOverlayPaths(dependencies)
	if err != nil {
		return nil, fmt.Errorf("failed to resolve dependency paths: %w", err)
	}

	overlays := make([]string, len(depPaths))
	for i, depPath := range depPaths {
		overlays[i] = depPath + ":ro"
	}
	return overlays, nil
}
