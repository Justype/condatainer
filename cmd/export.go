package cmd

import "github.com/spf13/cobra"

// overlayExportHelp is shared by 'export' and 'overlay export', which run the same command.
const overlayExportHelp = `Export the recipe that produced an overlay:

  - a def or build-script overlay prints its embedded recipe (Conda flags ignored)
  - a conda overlay is exported with micromamba env export

Output goes to stdout, or to <prefix>.<ext> with --prefix
  Apptainer (.def), build script (.sh), Conda (.yml or .txt)`

// overlayExportExample is shared by 'export' and 'overlay export'.
const overlayExportExample = `  condatainer export samtools/1.22                 # Conda overlay -> environment.yml on stdout
  condatainer export samtools/1.22 --from-history  # Conda: only explicitly-installed specs
  condatainer export cellranger/9.0.1              # Script overlay -> its build script
  condatainer export ubuntu24/base_image           # Def overlay -> its definition
  condatainer export ./env.img -p ./env            # Write ./env.yml (extension by type)`

var exportOverlayCmd = &cobra.Command{
	Use:               "export [flags] <overlay>",
	Short:             "Export the recipe that produced an overlay",
	Long:              overlayExportHelp,
	Example:           overlayExportExample,
	Args:              cobra.ExactArgs(1),
	SilenceUsage:      true,
	ValidArgsFunction: completeOverlayArg,
	RunE:              runExportOverlay,
}

// registerExportFlags registers the micromamba env-export flags shared by
// 'export' and 'overlay export'.
func registerExportFlags(cmd *cobra.Command) {
	// Keep registration order so --prefix (applies to every export type) leads,
	// ahead of the Conda-only flags, instead of cobra's alphabetical sort.
	cmd.Flags().SortFlags = false
	cmd.Flags().StringP("prefix", "p", "", "Write to <prefix>.<ext>; extension set by type")
	cmd.Flags().BoolP("explicit", "e", false, "Conda: use explicit format")
	cmd.Flags().Bool("no-md5", false, "Conda: disable md5")
	cmd.Flags().Bool("no-build", false, "Conda: disable the build string in spec")
	cmd.Flags().Bool("no-builds", false, "Conda: disable the build string in spec (alias)")
	cmd.Flags().Bool("channel-subdir", false, "Conda: enable channel/subdir in spec")
	cmd.Flags().Bool("from-history", false, "Conda: build environment spec from history")
}

func init() {
	rootCmd.AddCommand(exportOverlayCmd)
	registerExportFlags(exportOverlayCmd)
}
