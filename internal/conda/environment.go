package conda

import (
	"fmt"
	"os"
	"path/filepath"
	"strings"

	"github.com/Justype/condatainer/internal/toolpath"
)

// Environment is the mounted conda prefix managed by an in-container command.
type Environment struct {
	Root            string
	CondarcPath     string
	PinnedPath      string
	Micromamba      string
	DefaultChannels []string
	Writable        bool
}

// ResolveEnvironment validates the in-container runtime contract. Read-only
// commands need a mounted prefix; mutating commands also require the separate
// writable marker. CNT_CONDA_ROOT is the authoritative project path.
func ResolveEnvironment(requireWritable bool) (*Environment, error) {
	return resolveEnvironment(requireWritable, false)
}

// ResolveInstallEnvironment allows install to create the project root in a
// writable blank image. No other command may initialize the environment as
// a side effect.
func ResolveInstallEnvironment() (*Environment, error) {
	return resolveEnvironment(true, true)
}

func resolveEnvironment(requireWritable, allowInitialize bool) (*Environment, error) {
	if os.Getenv("IN_CONDATAINER") == "" {
		return nil, fmt.Errorf("env commands must run inside CondaTainer")
	}

	root := filepath.Clean(os.Getenv("CNT_CONDA_ROOT"))
	if root == "." || root == "" {
		return nil, fmt.Errorf("no environment overlay is mounted (CNT_CONDA_ROOT is unset)")
	}
	writable := os.Getenv("CNT_CONDA_WRITABLE") == "1"
	if requireWritable && !writable {
		return nil, fmt.Errorf("the environment is read-only")
	}

	info, err := os.Stat(root)
	if os.IsNotExist(err) && allowInitialize {
		if err := os.MkdirAll(root, 0o755); err != nil {
			return nil, fmt.Errorf("initialize environment directory at %s: %w", root, err)
		}
		info, err = os.Stat(root)
	}
	if os.IsNotExist(err) {
		return nil, fmt.Errorf("the environment is not initialized; initialize it in a writable container with: mm install <package>")
	}
	if err != nil {
		return nil, fmt.Errorf("cannot access environment at %s: %w", root, err)
	}
	if !info.IsDir() {
		return nil, fmt.Errorf("conda environment path is not a directory: %s", root)
	}

	mmPath, err := toolpath.Resolve("micromamba")
	if err != nil {
		return nil, err
	}

	env := &Environment{
		Root:            root,
		CondarcPath:     filepath.Join(root, ".condarc"),
		PinnedPath:      filepath.Join(root, "conda-meta", "pinned"),
		Micromamba:      mmPath,
		DefaultChannels: splitRuntimeChannels(os.Getenv("CNT_CONDA_CHANNELS")),
		Writable:        writable,
	}
	if !allowInitialize {
		state, err := env.State()
		if err != nil {
			return nil, fmt.Errorf("inspect environment at %s: %w", root, err)
		}
		switch state {
		case PrefixFresh:
			return nil, fmt.Errorf("the environment is not initialized; initialize it in a writable container with: mm install <package>")
		case PrefixPartial:
			return nil, fmt.Errorf("the environment is only partially initialized; recreate the overlay or remove %s and run mm install again", root)
		}
	}
	return env, nil
}

func splitRuntimeChannels(value string) []string {
	var channels []string
	for _, channel := range strings.Split(value, "|") {
		if channel = strings.TrimSpace(channel); channel != "" {
			channels = append(channels, channel)
		}
	}
	return deduplicate(channels)
}

// PrefixState reports whether a prefix is fresh, initialized, or partial.
type PrefixState int

const (
	PrefixFresh PrefixState = iota
	PrefixInitialized
	PrefixPartial
)

func (e *Environment) State() (PrefixState, error) {
	meta := filepath.Join(e.Root, "conda-meta")
	info, err := os.Stat(meta)
	if os.IsNotExist(err) {
		return PrefixFresh, nil
	}
	if err != nil {
		return PrefixPartial, err
	}
	if !info.IsDir() {
		return PrefixPartial, nil
	}
	if _, err := os.Stat(filepath.Join(meta, "history")); err == nil {
		return PrefixInitialized, nil
	} else if !os.IsNotExist(err) {
		return PrefixPartial, err
	}
	return PrefixPartial, nil
}
