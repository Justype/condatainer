package cmd

import (
	"errors"
	"fmt"

	"github.com/Justype/condatainer/internal/config"
	"github.com/Justype/condatainer/internal/store"
	"github.com/Justype/condatainer/internal/utils"
	"github.com/spf13/cobra"
)

func newStoreUseCmd() *cobra.Command {
	var identity string
	var jsonOutput bool
	cmd := &cobra.Command{
		Use:   "use <name>@<identity>",
		Short: "Make one build the default for its name",
		Long: `Makes the named identity the build that a bare name resolves to: it takes the
plain name, and whichever build held it moves into the same directory's store.

Only renames, and only inside the directory the chosen build already lives in.
Nothing leaves that directory, so every identity stays available to everyone who
reads it — a build is only ever addressed differently, never removed.

A build in a directory that a nearer one already shadows is refused, because
promoting it there would change nothing.`,
		Example: `  condatainer store use star/2.7.11b@9f2c1ab
  condatainer store use star/2.7.11b --identity 9f2c1ab`,
		Args:         cobra.ExactArgs(1),
		SilenceUsage: true,
		RunE: func(cmd *cobra.Command, args []string) error {
			name, query, err := splitStoreAddress(args[0], identity)
			if err != nil {
				return err
			}
			if err := confirmSharedPromotion(cmd, name, query); err != nil {
				return err
			}
			result, err := store.Promote(name, query, nil)
			if err != nil {
				return err
			}
			if jsonOutput {
				return printJSON(result)
			}
			reportPromotion(result)
			return nil
		},
	}
	cmd.Flags().StringVar(&identity, "identity", "", "Complete identity or unambiguous digest prefix")
	cmd.Flags().BoolVar(&jsonOutput, "json", false, "Print JSON")
	return cmd
}

// confirmSharedPromotion asks before changing what a name means for other people.
//
// A promotion breaks no project lock — those pin identities, and every identity
// stays resolvable — but everyone addressing the name gets a different build
// from then on. In a personal directory there is nobody else to surprise.
func confirmSharedPromotion(cmd *cobra.Command, name string, query store.IdentityQuery) error {
	if utils.ShouldAnswerYes() {
		return nil
	}
	candidate, _, err := store.ResolveIdentity(name, query, nil)
	if err != nil || candidate.Layout == store.LayoutFlat || config.IsPersonalImagesDir(candidate.Root) {
		// A resolution failure is Promote's to report, with its own message.
		return nil
	}
	utils.PrintWarning("%s is shared: everyone reading it gets this build for %s from now on.",
		utils.StylePath(candidate.Root), utils.StyleName(candidate.Name))
	fmt.Printf("Promote anyway? [y/N]: ")
	choice, err := utils.ReadLineContext(cmd.Context())
	if err != nil || (choice != "y" && choice != "yes") {
		return errors.New("cancelled")
	}
	return nil
}

func reportPromotion(result store.Promotion) {
	if result.AlreadyFlat {
		utils.PrintSuccess("%s already resolves to this build (%s).",
			utils.StyleName(result.Name), store.FormatKeyRef(result.Promoted.Identity))
		return
	}
	utils.PrintSuccess("%s now resolves to %s", utils.StyleName(result.Name),
		utils.StylePath(result.Promoted.Path))
	utils.PrintMessage("  identity %s", store.FormatKeyRef(result.Promoted.Identity))
	if result.DemotedTo != "" {
		utils.PrintMessage("  was %s, now at %s",
			store.FormatKeyRef(result.Demoted.Identity), utils.StylePath(result.DemotedTo))
	}
	// Named because it is still there: a farther copy is shadowed, not removed,
	// and someone reading only that directory still gets it.
	if result.Shadows != "" {
		utils.PrintNote("%s is now shadowed by this one.", utils.StylePath(result.Shadows))
	}
}
