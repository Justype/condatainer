package cmd

import (
	"fmt"
	"os"

	condapkg "github.com/Justype/condatainer/internal/conda"
	"github.com/spf13/cobra"
)

var envCmd = &cobra.Command{
	Use:   "env",
	Short: "Manage the Conda environment of the mounted overlay",
	Long: `Install, update and inspect the Conda packages of the environment overlay
mounted in the current container.

These commands run inside CondaTainer only. Enter a writable container first:
  condatainer exec -w -o env.img bash`,
	SilenceUsage: true,
	RunE: func(cmd *cobra.Command, args []string) error {
		return cmd.Help()
	},
}

func init() {
	rootCmd.AddCommand(envCmd)

	envCmd.AddCommand(newMicromambaCommand("install", "Install packages into the environment", runEnvInstall))
	envCmd.AddCommand(newMicromambaCommand("update", "Update installed packages", runEnvMutation("update")))
	envCmd.AddCommand(newMicromambaCommand("remove", "Remove installed packages", runEnvMutation("remove")))
	envCmd.AddCommand(newMicromambaCommand("list", "List installed packages", runEnvRead("list")))
	envCmd.AddCommand(newMicromambaCommand("export", "Print the environment specification", runEnvRead("env", "export")))
	envCmd.AddCommand(newMicromambaCommand("search", "Search the channels for a package", runEnvMutation("search")))
	envCmd.AddCommand(newMicromambaCommand("clean", "Clean the package caches", runEnvMutation("clean")))

	channelsCmd := &cobra.Command{Use: "channels", Short: "Manage the environment's Conda channels"}
	channelsCmd.AddCommand(
		envStructuredCommand("list", "List the channels, highest priority first", cobra.NoArgs, runEnvChannelsList),
		envStructuredCommand("append <channel>", "Add a channel with the lowest priority", cobra.ExactArgs(1), runEnvChannelsAppend),
		envStructuredCommand("prepend <channel>", "Add a channel with the highest priority", cobra.ExactArgs(1), runEnvChannelsPrepend),
		envStructuredCommand("remove <channel>", "Remove a channel", cobra.ExactArgs(1), runEnvChannelsRemove),
	)
	envCmd.AddCommand(channelsCmd)

	pinCmd := &cobra.Command{Use: "pin", Short: "Keep packages at a chosen version"}
	pinCmd.AddCommand(
		envStructuredCommand("list", "List the pinned packages", cobra.NoArgs, runEnvPinList),
		envStructuredCommand("add <spec>", "Pin a package; a bare name pins the installed version", cobra.ExactArgs(1), runEnvPinAdd),
		envStructuredCommand("remove <package>", "Unpin a package", cobra.ExactArgs(1), runEnvPinRemove),
	)
	envCmd.AddCommand(pinCmd)
}

func envStructuredCommand(use, short string, args cobra.PositionalArgs, run envRunFunc) *cobra.Command {
	return &cobra.Command{Use: use, Short: short, Args: args, RunE: run, SilenceUsage: true}
}

type envRunFunc func(*cobra.Command, []string) error

func newMicromambaCommand(name, short string, run envRunFunc) *cobra.Command {
	return &cobra.Command{
		Use:                name + " [micromamba arguments...]",
		Short:              short,
		DisableFlagParsing: true,
		SilenceUsage:       true,
		RunE:               run,
	}
}

func envIO(cmd *cobra.Command) condapkg.IO {
	return condapkg.IO{Stdin: os.Stdin, Stdout: cmd.OutOrStdout(), Stderr: cmd.ErrOrStderr()}
}

func runEnvInstall(cmd *cobra.Command, args []string) error {
	env, err := condapkg.ResolveInstallEnvironment()
	if err != nil {
		return err
	}
	return env.Install(cmd.Context(), args, envIO(cmd))
}

func runEnvMutation(operation ...string) envRunFunc {
	return func(cmd *cobra.Command, args []string) error {
		env, err := condapkg.ResolveEnvironment(true)
		if err != nil {
			return err
		}
		commandArgs := append([]string{}, operation...)
		commandArgs = append(commandArgs, args...)
		return env.Run(cmd.Context(), envIO(cmd), commandArgs...)
	}
}

func runEnvRead(operation ...string) envRunFunc {
	return func(cmd *cobra.Command, args []string) error {
		env, err := condapkg.ResolveEnvironment(false)
		if err != nil {
			return err
		}
		commandArgs := append([]string{}, operation...)
		commandArgs = append(commandArgs, args...)
		return env.Run(cmd.Context(), envIO(cmd), commandArgs...)
	}
}

func runEnvChannelsList(cmd *cobra.Command, _ []string) error {
	env, err := condapkg.ResolveEnvironment(false)
	if err != nil {
		return err
	}
	channels, err := condapkg.ReadChannels(env.CondarcPath)
	if err != nil {
		return err
	}
	if len(channels) == 0 {
		return fmt.Errorf("no channels are configured in %s", env.CondarcPath)
	}
	for _, channel := range channels {
		fmt.Fprintln(cmd.OutOrStdout(), channel)
	}
	return nil
}

func runEnvChannelsAppend(cmd *cobra.Command, args []string) error {
	return changeChannels(cmd, args[0], false, false)
}

func runEnvChannelsPrepend(cmd *cobra.Command, args []string) error {
	return changeChannels(cmd, args[0], true, false)
}

func runEnvChannelsRemove(cmd *cobra.Command, args []string) error {
	return changeChannels(cmd, args[0], false, true)
}

func changeChannels(cmd *cobra.Command, channel string, prepend, remove bool) error {
	env, err := condapkg.ResolveEnvironment(true)
	if err != nil {
		return err
	}
	change := condapkg.ChannelAppend
	if remove {
		change = condapkg.ChannelRemove
	} else if prepend {
		change = condapkg.ChannelPrepend
	}
	found, err := env.ChangeChannel(channel, change)
	if err != nil {
		return err
	}
	switch {
	case remove:
		fmt.Fprintf(cmd.OutOrStdout(), "Removed channel: %s\n", channel)
	case found:
		fmt.Fprintf(cmd.OutOrStdout(), "Moved channel: %s\n", channel)
	default:
		fmt.Fprintf(cmd.OutOrStdout(), "Added channel: %s\n", channel)
	}
	return nil
}

func runEnvPinList(cmd *cobra.Command, _ []string) error {
	env, err := condapkg.ResolveEnvironment(false)
	if err != nil {
		return err
	}
	pins, err := condapkg.ReadPins(env.PinnedPath)
	if err != nil {
		return err
	}
	if len(pins) == 0 {
		fmt.Fprintln(cmd.OutOrStdout(), "No packages are pinned.")
		return nil
	}
	for _, pin := range pins {
		fmt.Fprintln(cmd.OutOrStdout(), pin)
	}
	return nil
}

func runEnvPinAdd(cmd *cobra.Command, args []string) error {
	env, err := condapkg.ResolveEnvironment(true)
	if err != nil {
		return err
	}
	spec, err := env.Pin(cmd.Context(), args[0])
	if err != nil {
		return err
	}
	fmt.Fprintf(cmd.OutOrStdout(), "Pinned: %s\n", spec)
	return nil
}

func runEnvPinRemove(cmd *cobra.Command, args []string) error {
	env, err := condapkg.ResolveEnvironment(true)
	if err != nil {
		return err
	}
	if err := env.Unpin(args[0]); err != nil {
		return err
	}
	fmt.Fprintf(cmd.OutOrStdout(), "Unpinned: %s\n", args[0])
	return nil
}
