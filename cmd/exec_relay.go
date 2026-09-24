package cmd

import (
	"os"

	"github.com/Justype/condatainer/internal/runtime/proxy"
	"github.com/spf13/cobra"
)

func init() {
	rootCmd.AddCommand(execRelayCmd)
}

var execRelayCmd = &cobra.Command{
	Use:    "_exec_relay",
	Hidden: true,
	Short:  "Relay loopback connections for the dashboard (run inside a helper's job)",
	Args:   cobra.NoArgs,
	RunE: func(cmd *cobra.Command, args []string) error {
		cmd.SilenceUsage = true
		return proxy.RunExecRelay(os.Stdin, os.Stdout)
	},
}
