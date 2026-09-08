package cmd

import (
	"os"

	"github.com/Justype/condatainer/internal/image/freeze"
	"github.com/spf13/cobra"
)

func init() {
	rootCmd.AddCommand(mountSentinelCmd)
}

var mountSentinelCmd = &cobra.Command{
	Use:    "_mount_sentinel <fuseBin> <mnt> [fuseArgs...]",
	Hidden: true,
	Short:  "Run the FUSE-mount sentinel condatainer re-execs itself into (internal)",
	Args:   cobra.MinimumNArgs(2),
	RunE: func(cmd *cobra.Command, args []string) error {
		cmd.SilenceUsage = true
		fuseBin, mnt := args[0], args[1]
		fuseArgs := args[2:]
		return freeze.RunSentinel(fuseBin, mnt, os.Getenv(freeze.EnvSentinelWork), fuseArgs)
	},
}
