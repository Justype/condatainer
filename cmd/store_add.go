package cmd

import (
	"fmt"

	"github.com/Justype/condatainer/internal/artifact/compare"
	"github.com/Justype/condatainer/internal/artifact/meta"
	"github.com/Justype/condatainer/internal/config"
	"github.com/Justype/condatainer/internal/store"
	"github.com/Justype/condatainer/internal/utils"
	"github.com/spf13/cobra"
)

func newStoreAddCmd() *cobra.Command {
	var layer string
	cmd := &cobra.Command{
		Use:   "add <file.sqf>",
		Short: "Install an existing overlay file into the store",
		Long: `Copies one .sqf into the store, filed under the identity the file itself
records — the source filename is never a naming claim.

--layer chooses which store, not what the file is called: inside it the name
comes from the keys. To put an overlay at a path of your choosing, use
'create -p' or 'registry pull -p'.

Copied, never moved or linked: the source file stays exactly where it is.`,
		Example: `  condatainer store add ./star-2.7.11b.sqf
  condatainer store add /scratch/builds/star.sqf --layer user`,
		Args:         cobra.ExactArgs(1),
		SilenceUsage: true,
		RunE: func(cmd *cobra.Command, args []string) error {
			source := args[0]
			if !utils.IsSqf(source) {
				return fmt.Errorf("%s is not a .sqf: only a read-only overlay has an identity to file it under", source)
			}
			dir, err := storeAddDir(layer)
			if err != nil {
				return err
			}
			artifact, err := compare.Read(source)
			if err != nil {
				return fmt.Errorf("cannot read the identity of %s: %w", source, err)
			}
			// The store files an artifact under its keys, and a frozen environment
			// has none — it is a project's own environment, addressed by path.
			// Filed here it would land under an empty identity, and every project's
			// would land under the same name.
			if artifact.Format == meta.BuildTypeSnapshot {
				return fmt.Errorf("%s is a frozen environment: it belongs to its project and is addressed by path, so it is never filed in the store", source)
			}

			// A named layer scopes the search for an existing copy to that layer.
			// Left to search every root, an identity installed anywhere would be
			// adopted and nothing written — which is right when no destination was
			// asked for, and wrong when one was: copying a build into the directory
			// that owns a name is exactly how `store use` says to lift a shadow.
			opts := store.BeginOptions{ImagesDir: dir, Equiv: artifact.EquivRef(), StoreOnly: true}
			if cmd.Flags().Changed("layer") {
				opts.SearchDirs = []string{dir}
			}

			candidate, err := store.InstallFile(artifact.Name, artifact.IdentityRef(), source, opts)
			if err != nil {
				return err
			}
			if candidate.Path != source && candidate.Layout == store.LayoutStored &&
				candidate.Root == dir {
				utils.PrintSuccess("Added %s", utils.StylePath(candidate.Path))
			} else {
				// Begin adopted a copy that was already installed, which may be in
				// another directory entirely. Say where, so "nothing was written"
				// is never read as "installed here".
				utils.PrintSuccess("Already installed at %s", utils.StylePath(candidate.Path))
			}
			utils.PrintMessage("  %s  %s", utils.StyleName(candidate.Name), store.FormatKeyRef(candidate.Identity))
			return nil
		},
	}
	cmd.Flags().StringVarP(&layer, "layer", "l", "", "Add to this data layer: u/user, r/app-root, e/extra-root")
	cmd.RegisterFlagCompletionFunc("layer", //nolint:errcheck
		func(*cobra.Command, []string, string) ([]string, cobra.ShellCompDirective) {
			return []string{"user", "app-root", "extra-root"}, cobra.ShellCompDirectiveNoFileComp
		})
	return cmd
}

// storeAddDir resolves which images directory receives the artifact.
func storeAddDir(layer string) (string, error) {
	if layer == "" {
		dir, err := config.GetWritableImagesDir()
		if err != nil {
			return "", fmt.Errorf("no writable images directory found: %w", err)
		}
		return dir, nil
	}
	selected, err := config.ParseDataLayer(layer)
	if err != nil {
		return "", err
	}
	return config.GetWritableImagesDirIn(selected)
}
