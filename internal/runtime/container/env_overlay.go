package container

import (
	"path/filepath"

	"github.com/Justype/condatainer/internal/utils"
)

// ResolveEnvOverlay resolves the project's environment overlay in cwd: the
// writable .img if one exists, else the read-only env-typed .sqf beside
// where it would go. Never creates anything. Callers check utils.IsImg on
// the result to tell which form was found.
//
// Lives here rather than in internal/helper so internal/project can call it
// too (checking whether a frozen env.sqf is pinned) without importing
// internal/helper, which itself imports internal/project.
func ResolveEnvOverlay(envImg, cwd string) string {
	if p := utils.FindEnvOverlay(envImg, cwd); p != "" {
		return p
	}
	if envImg != "" && envImg != "env.img" {
		return ""
	}
	return FindEnvSnapshot(cwd)
}

// FindEnvSnapshot looks in cwd for an env-typed .sqf: a personal
// env-$USER.sqf line checked before the shared env.sqf line. Returns "" if
// neither exists.
func FindEnvSnapshot(cwd string) string {
	wd := utils.ResolveWD(cwd)
	return LookupSnapshot(filepath.Join(wd, "env.img")).Path
}
