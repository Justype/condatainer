package cmd

import (
	"fmt"
	"strconv"

	"github.com/Justype/condatainer/internal/artifact/key"
	"github.com/spf13/cobra"
)

func init() {
	rootCmd.AddCommand(payloadKeyCmd)
}

var payloadKeyCmd = &cobra.Command{
	Use:    "_payload_key <dir> <base> <workers>",
	Hidden: true,
	Short:  "Print the payload key of a directory (internal)",
	Args:   cobra.ExactArgs(3),
	RunE: func(cmd *cobra.Command, args []string) error {
		cmd.SilenceUsage = true
		workers, err := strconv.Atoi(args[2])
		if err != nil {
			return fmt.Errorf("workers: %w", err)
		}
		records, err := key.TreeOf(cmd.Context(), args[0], args[1], workers)
		if err != nil {
			return err
		}
		ref, err := key.PayloadKey(records)
		if err != nil {
			return err
		}
		fmt.Fprintln(cmd.OutOrStdout(), ref.SHA256)
		return nil
	},
}
