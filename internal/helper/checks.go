package helper

import (
	"context"
	"os"

	"github.com/Justype/condatainer/internal/image"
	"github.com/Justype/condatainer/internal/logging"
	"github.com/Justype/condatainer/internal/project"
	"github.com/Justype/condatainer/internal/runtime/container"
	"github.com/Justype/condatainer/internal/scheduler"
	"github.com/Justype/condatainer/internal/utils"
)

// DefaultOverlaySize is the size offered for a new writable conda overlay
// when nothing else is specified — the CLI's guided-creation prompt default,
// the dashboard's create-form default, and the size used for the no-prompt
// overlay created on top of a found snapshot (the third overlay state) all
// read this one value rather than each hardcoding "20G" separately.
const DefaultOverlaySize = "20G"

// EnvStatus describes the state of an env overlay path for UI display.
// Returned by CheckEnv. Both CLI and server consume this directly.
//
// SizeMB and Snapshot describe the pair, not the .img alone: a thin .img on
// top of a multi-gigabyte snapshot would otherwise report a misleadingly
// small size, and there would be no way to tell which snapshot it continues
// from at all.
type EnvStatus struct {
	Path     string `json:"path"`
	Exists   bool   `json:"exists"`
	InUse    bool   `json:"in_use"`
	Writable bool   `json:"writable"`
	SizeMB   int64  `json:"size_mb,omitempty"`
	Snapshot string `json:"snapshot,omitempty"`
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
//
// SizeMB and Snapshot come from container.PairedSize, which is looked up
// whether or not path exists on disk — this is what lets a caller tell the
// third overlay state apart from "nothing here at all": no .img, but a
// fully-populated snapshot right beside where one would go.
func CheckEnv(ctx context.Context, path string) (EnvStatus, error) {
	if path == "" {
		return EnvStatus{}, nil
	}
	st := EnvStatus{Path: path}
	isImg := utils.IsImg(path)

	if _, err := os.Stat(path); err == nil {
		st.Exists = true
		st.Writable = isImg
		if err := image.CheckAvailable(path, st.Writable); err != nil {
			st.InUse = true
		}
	} else if !os.IsNotExist(err) {
		return st, err
	}

	sizeBytes, snapshot := container.PairedSize(path)
	st.SizeMB = sizeBytes / (1024 * 1024)
	st.Snapshot = snapshot
	return st, nil
}

// MissingParams returns the #PARAM: entries from scriptPath that need user
// input — i.e. not present in `supplied`, not optional (KEY=?), and without
// a literal default.
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
		// Optional params (KEY=?) auto-fill from #VALUE: or pass empty.
		if p.Optional {
			continue
		}
		// Literal defaults auto-fill — not missing.
		if p.Default != "" {
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
// params and returns the resolved absolute paths.
//
// cwd standing in a project resolves every name through that project's lock
// instead of by installed name, building nothing — a helper's required
// overlays reach the lock only as manual pins, so an unresolved one names
// `project restore`. Outside a project, when cwd is not one, or when
// noProject is set, each name is ensured on disk the ordinary way, building
// any missing one via `condatainer create`.
//
// Returns (nil, nil) when the template is empty. Public wrapper around the
// previously unexported checkAndInstallNamedOverlays.
func CheckRequiredOverlays(ctx context.Context, cwd, requiredTemplate string, params map[string]string, noProject bool) ([]string, error) {
	if requiredTemplate == "" {
		return nil, nil
	}
	logger := logging.FromContext(ctx)
	names := resolveOverlayTemplate(requiredTemplate, params)
	if !noProject {
		standing, err := project.StandingAt(cwd)
		if err != nil {
			return nil, err
		}
		if standing != nil {
			logger.Info("Checking required overlays", "project", standing.Root)
			return standing.ResolveNames(names)
		}
	}
	logger.Info("Checking required overlays")
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
