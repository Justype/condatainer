package conda

import (
	"bytes"
	"context"
	"encoding/json"
	"fmt"
	"io"
	"os"
	"os/exec"
	"slices"
	"strings"

	"go.yaml.in/yaml/v3"
)

type IO struct {
	Stdin  io.Reader
	Stdout io.Writer
	Stderr io.Writer
}

func (e *Environment) commandEnv() []string {
	managed := map[string]bool{
		"CONDA_PREFIX": true, "CONDARC": true, "MAMBA_NO_RC": true,
		"MAMBA_ROOT_PREFIX": true, "MAMBA_TARGET_PREFIX": true,
	}
	env := make([]string, 0, len(os.Environ())+4)
	for _, entry := range os.Environ() {
		name, _, _ := strings.Cut(entry, "=")
		if !managed[name] {
			env = append(env, entry)
		}
	}
	env = append(env,
		"CONDA_PREFIX="+e.Root,
		"MAMBA_ROOT_PREFIX="+e.Root,
	)
	return env
}

func (e *Environment) command(ctx context.Context, command ...string) *exec.Cmd {
	configArgs := []string{"--rc-file", e.CondarcPath}
	if _, err := os.Stat(e.CondarcPath); os.IsNotExist(err) {
		configArgs = []string{"--no-rc"}
	}
	commandArgs := append(configArgs, command...)
	cmd := exec.CommandContext(ctx, e.Micromamba, commandArgs...)
	cmd.Env = e.commandEnv()
	return cmd
}

func (e *Environment) Run(ctx context.Context, streams IO, command ...string) error {
	cmd := e.command(ctx, command...)
	cmd.Stdin = streams.Stdin
	cmd.Stdout = streams.Stdout
	cmd.Stderr = streams.Stderr
	if err := cmd.Run(); err != nil {
		if ctx.Err() != nil {
			return ctx.Err()
		}
		return err
	}
	return nil
}

func (e *Environment) runOutput(ctx context.Context, command ...string) ([]byte, error) {
	cmd := e.command(ctx, command...)
	var stderr bytes.Buffer
	cmd.Stderr = &stderr
	out, err := cmd.Output()
	if err != nil {
		if ctx.Err() != nil {
			return nil, ctx.Err()
		}
		if stderr.Len() > 0 {
			return nil, fmt.Errorf("%s: %w", strings.TrimSpace(stderr.String()), err)
		}
		return nil, err
	}
	return out, nil
}

// Install lets Micromamba initialize a missing root or modify an existing one.
func (e *Environment) Install(ctx context.Context, args []string, streams IO) error {
	if _, err := os.Stat(e.CondarcPath); err == nil {
		command := append([]string{"install"}, args...)
		return e.Run(ctx, streams, command...)
	} else if !os.IsNotExist(err) {
		return err
	}

	channels := projectChannels(args, e.DefaultChannels)
	command := []string{"install"}
	for _, channel := range channels {
		command = append(command, "-c", channel)
	}
	command = append(command, args...)
	if err := e.Run(ctx, streams, command...); err != nil {
		return err
	}
	return WriteChannels(e.CondarcPath, channels)
}

// projectChannels derives the persistent project configuration without
// changing the arguments passed to Micromamba. A YAML environment file owns
// its channel list; otherwise explicit CLI channels are saved. The defaults
// match the original mm-create policy.
func projectChannels(args, defaultChannels []string) []string {
	var cliChannels []string
	var environmentFile string
	for i := 0; i < len(args); i++ {
		switch {
		case args[i] == "-c" || args[i] == "--channel":
			if i+1 < len(args) {
				i++
				cliChannels = append(cliChannels, args[i])
			}
		case strings.HasPrefix(args[i], "--channel="):
			cliChannels = append(cliChannels, strings.TrimPrefix(args[i], "--channel="))
		case strings.HasPrefix(args[i], "-c="):
			cliChannels = append(cliChannels, strings.TrimPrefix(args[i], "-c="))
		case args[i] == "-f" || args[i] == "--file":
			if i+1 < len(args) {
				i++
				environmentFile = args[i]
			}
		case strings.HasPrefix(args[i], "--file="):
			environmentFile = strings.TrimPrefix(args[i], "--file=")
		case strings.HasPrefix(args[i], "-f="):
			environmentFile = strings.TrimPrefix(args[i], "-f=")
		}
	}

	if environmentFile != "" {
		if data, err := os.ReadFile(environmentFile); err == nil {
			var spec struct {
				Channels []string `yaml:"channels"`
			}
			if yaml.Unmarshal(data, &spec) == nil && len(spec.Channels) > 0 {
				return deduplicate(spec.Channels)
			}
		}
		return deduplicate(defaultChannels)
	}

	channels := deduplicate(cliChannels)
	if len(channels) == 0 {
		return deduplicate(defaultChannels)
	}
	for _, channel := range deduplicate(defaultChannels) {
		if !slices.Contains(channels, channel) {
			channels = append(channels, channel)
		}
	}
	return channels
}

func deduplicate(values []string) []string {
	seen := make(map[string]bool, len(values))
	out := make([]string, 0, len(values))
	for _, value := range values {
		if value != "" && !seen[value] {
			seen[value] = true
			out = append(out, value)
		}
	}
	return out
}

func (e *Environment) InstalledVersion(ctx context.Context, packageName string) (string, error) {
	out, err := e.runOutput(ctx, "list", "--json")
	if err != nil {
		return "", err
	}
	var packages []struct {
		Name    string `json:"name"`
		Version string `json:"version"`
	}
	if err := json.Unmarshal(out, &packages); err != nil {
		return "", fmt.Errorf("parse micromamba package list: %w", err)
	}
	for _, pkg := range packages {
		if pkg.Name == packageName {
			return pkg.Version, nil
		}
	}
	return "", fmt.Errorf("package %q is not installed", packageName)
}
