package helper

import (
	"context"
	"os"

	"github.com/Justype/condatainer/internal/overlay"
	"github.com/Justype/condatainer/internal/scheduler"
	"github.com/Justype/condatainer/internal/utils"
)

// EnvStatus describes the state of an env overlay path for UI display.
// Returned by CheckEnv. Both CLI and server consume this directly.
type EnvStatus struct {
	Path     string `json:"path"`
	Exists   bool   `json:"exists"`
	InUse    bool   `json:"in_use"`
	Writable bool   `json:"writable"`
	SizeMB   int64  `json:"size_mb,omitempty"`
}

// ResolveEnv returns the env overlay path that the folder-based convention
// would auto-pick for the given cwd, or "" if none exists.
// Thin wrapper around ResolveEnvOverlayInDir for symmetry with CheckEnv.
func ResolveEnv(cwd string) string {
	return ResolveEnvOverlayInDir("", cwd)
}

// CheckEnv inspects a single overlay path: existence, lock status, type, size.
// Returns (zero EnvStatus, nil) when path is empty so callers can probe an
// unresolved path without a special case. The lock probe uses the same kind
// (exclusive vs. shared) that an actual mount would request: exclusive for
// writable .img, shared for read-only .sqf.
func CheckEnv(ctx context.Context, path string) (EnvStatus, error) {
	if path == "" {
		return EnvStatus{}, nil
	}
	st := EnvStatus{Path: path}
	info, err := os.Stat(path)
	if err != nil {
		if os.IsNotExist(err) {
			return st, nil
		}
		return st, err
	}
	st.Exists = true
	st.SizeMB = info.Size() / (1024 * 1024)
	st.Writable = utils.IsImg(path)

	if err := overlay.CheckAvailable(path, st.Writable); err != nil {
		st.InUse = true
	}
	return st, nil
}

// MissingParams returns the #PARAM: entries from scriptPath that need user
// input — i.e. not present in `supplied` and without an auto-usable default.
// Version-token defaults (e.g. "{CONDA_PYTHON}") always count as missing
// because the user must pick a concrete version.
//
// Pure: no prompts, no I/O beyond reading the script file.
func MissingParams(scriptPath string, supplied map[string]string) ([]HelperParam, error) {
	params, err := ParseHelperParams(scriptPath)
	if err != nil {
		return nil, err
	}
	var missing []HelperParam
	for _, p := range params {
		if v, ok := supplied[p.Key]; ok && v != "" {
			continue
		}
		// Optional params (KEY=?) pass through empty — script handles them.
		if p.Optional {
			continue
		}
		// Non-version defaults auto-fill — not missing.
		if p.Default != "" && !versionTokenRe.MatchString(p.Default) {
			continue
		}
		missing = append(missing, p)
	}
	return missing, nil
}

// Running returns the active HelperRun records for name (empty = any).
// Symmetric alias for RunningHelpers.
func Running(name string) ([]*HelperRun, error) {
	return RunningHelpers(name)
}

// SingletonBlocked reports whether the script's #SINGLETON: meta forbids a
// new launch given the running set.
func SingletonBlocked(meta HelperScriptMeta, running []*HelperRun) bool {
	return meta.Singleton && len(running) > 0
}

// CheckRequiredOverlays expands {tokens} in meta.RequiredOverlays using
// params, ensures each overlay exists on disk (building any missing ones via
// `condatainer create`), and returns the resolved absolute paths.
//
// Returns (nil, nil) when the template is empty. Public wrapper around the
// previously unexported checkAndInstallNamedOverlays.
func CheckRequiredOverlays(ctx context.Context, requiredTemplate string, params map[string]string) ([]string, error) {
	if requiredTemplate == "" {
		return nil, nil
	}
	names := resolveOverlayTemplate(requiredTemplate, params)
	return checkAndInstallNamedOverlays(ctx, names)
}

// CheckHelperPackages verifies that every package in meta.ImgPackages is
// installed in the conda environment inside envImg. {KEY} tokens in ImgPackages
// are substituted from params before checking.
// Public wrapper around checkPackages so the server can run pre-submission
// validation without going through PlanRun.
func CheckHelperPackages(meta HelperScriptMeta, envImg string, params map[string]string) error {
	return checkPackages(meta, envImg, params)
}

// ResolveResources merges script headers with config defaults and explicit overrides.
// Pass a non-nil overrides to apply user-supplied resource values on top. Pure (no logger calls).
func ResolveResources(ctx context.Context, scriptPath string, overrides *scheduler.ResourceSpec) *scheduler.ResourceSpec {
	return resolveSpec(scriptPath, overrides)
}
