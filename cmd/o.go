package cmd

import (
	"github.com/spf13/cobra"
)

// oCmd is a shortcut for "overlay create"
var oCmd = &cobra.Command{
	Use:   "o [flags] [path] [-- packages...]",
	Short: "Shortcut for 'overlay create'",
	Long: overlayCreateHelp + `

This is a shortcut for 'overlay create'.`,
	Example: `  condatainer o # 10G with default inode ratio
  condatainer o my_data.img -s 50g -p large
  condatainer o --fakeroot --sparse
  condatainer o -f environment.yml
  condatainer o myenv.img -- python=3.11`,

	Args: cobra.ArbitraryArgs,

	Run: runOverlayCreate,
}

func init() {
	// 1. Attach to root command
	rootCmd.AddCommand(oCmd)

	// 2. Define flags (identical to overlay create)
	registerOverlayCreateFlags(oCmd)
}
