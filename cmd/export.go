package cmd

import "github.com/spf13/cobra"

// overlayExportHelp is the help body for 'overlay export'.
const overlayExportHelp = `Export the Conda environment in a writable .img overlay.

Runs micromamba env export against /cnt_env, so it captures the environment as
it is now — including anything installed into the overlay since it was created.

Installed .sqf and .sif images are not exportable: they are immutable and carry
their own metadata, so rebuild them from their recipe instead.

Output goes to stdout, or to <prefix>.yml (or .txt with --explicit) with --prefix.`

// overlayExportExample is the example block for 'overlay export'.
const overlayExportExample = `  condatainer overlay export ./env.img                    # environment.yml on stdout
  condatainer overlay export ./env.img --from-history     # only explicitly-installed specs
  condatainer overlay export ./env.img -e                 # explicit spec (.txt)
  condatainer overlay export ./env.img -p ./env           # write ./env.yml`

// registerExportFlags registers the micromamba env-export flags for 'overlay export'.
func registerExportFlags(cmd *cobra.Command) {
	// Keep registration order so --prefix leads, ahead of the Conda-only flags,
	// instead of cobra's alphabetical sort.
	cmd.Flags().SortFlags = false
	cmd.Flags().StringP("prefix", "p", "", "Write to <prefix>.<ext>; extension set by format")
	cmd.Flags().BoolP("explicit", "e", false, "Conda: use explicit format")
	cmd.Flags().Bool("no-md5", false, "Conda: disable md5")
	cmd.Flags().Bool("no-build", false, "Conda: disable the build string in spec")
	cmd.Flags().Bool("no-builds", false, "Conda: disable the build string in spec (alias)")
	cmd.Flags().Bool("channel-subdir", false, "Conda: enable channel/subdir in spec")
	cmd.Flags().Bool("from-history", false, "Conda: build environment spec from history")
}
