package cmd

import (
	"fmt"
	"os"
	"strings"
	"time"

	"github.com/Justype/condatainer/internal/config"
	"github.com/Justype/condatainer/internal/store"
	"github.com/Justype/condatainer/internal/utils"
	"github.com/spf13/cobra"
)

func newStoreRmCmd() *cobra.Command {
	var identity string
	cmd := &cobra.Command{
		Use:   "rm <name>@<identity>",
		Short: "Remove one store entry",
		Long: `Deletes one build from the store, chosen by name and identity.

Store entries only. An overlay under a plain name belongs to
'condatainer remove'.

The entry has to be free: removal fails while a container is reading it, and an
overlay whose write bit was cleared to pin it is never touched.`,
		Example: `  condatainer store rm star/2.7.11b@9f2c1ab
  condatainer store rm star/2.7.11b --identity 9f2c1ab`,
		Args:         cobra.ExactArgs(1),
		SilenceUsage: true,
		RunE: func(_ *cobra.Command, args []string) error {
			name, query, err := splitStoreAddress(args[0], identity)
			if err != nil {
				return err
			}
			removed, err := store.Remove(name, query, nil)
			if err != nil {
				return err
			}
			utils.PrintSuccess("Removed %s (%s)", utils.StylePath(removed.Path),
				store.FormatKeyRef(removed.Identity))
			return nil
		},
	}
	cmd.Flags().StringVar(&identity, "identity", "", "Complete identity or unambiguous digest prefix")
	return cmd
}

func newStoreGCCmd() *cobra.Command {
	var (
		apply    bool
		graceDay int
		dir      string
		layer    string
		jsonOut  bool
	)
	cmd := &cobra.Command{
		Use:   "gc",
		Short: "Report reclaimable store entries, or delete them",
		Long: `Reports which store entries could be deleted, and how much space that would
free. Nothing is deleted without --apply, which repeats the whole judgement
before removing anything.

An entry is collectable when it is older than the grace period, nothing is
reading it, and it still matches its recorded identity. Anything uncertain is
kept.

Age is the newest of the file's access, modification and change times, and each
row names which one decided it. A report where every row says "ctime" means the
filesystem does not record access times, not that the artifacts are cold.

-D/--dir and -l/--layer narrow the run, exactly as they do for 'list'. --apply
requires one of them, so a group directory shared with people who are not at
the keyboard is never collected by omission.`,
		Example: `  condatainer store gc
  condatainer store gc --layer user --grace 60
  condatainer store gc --dir scratch --apply`,
		Args:         cobra.NoArgs,
		SilenceUsage: true,
		RunE: func(_ *cobra.Command, _ []string) error {
			dirs, err := storeGCDirs(dir, layer)
			if err != nil {
				return err
			}
			report, err := store.GC(store.GCOptions{
				Dirs:  dirs,
				Grace: time.Duration(graceDay) * 24 * time.Hour,
				Apply: apply,
			})
			if err != nil {
				return err
			}
			return printGCReport(report, jsonOut)
		},
	}
	cmd.Flags().BoolVar(&apply, "apply", false, "Delete the collectable entries; requires --dir or --layer")
	cmd.Flags().IntVar(&graceDay, "grace", 0, "Age in days below which nothing is collectable (default: store_gc_grace)")
	cmd.Flags().StringVarP(&dir, "dir", "D", "", "Limit to image directories matching this substring")
	cmd.Flags().StringVarP(&layer, "layer", "l", "", "Limit to a data layer: u/user, r/app-root, e/extra-root")
	cmd.Flags().BoolVar(&jsonOut, "json", false, "Print JSON")
	cmd.RegisterFlagCompletionFunc("layer", //nolint:errcheck
		func(*cobra.Command, []string, string) ([]string, cobra.ShellCompDirective) {
			return []string{"user", "app-root", "extra-root"}, cobra.ShellCompDirectiveNoFileComp
		})
	return cmd
}

// storeGCDirs narrows the writable image roots to what --dir and --layer name.
//
// Empty means every writable root, which only a report may have: store.GC
// refuses an applied run with no scope. A filter that matches nothing is an
// error rather than an empty report, since "nothing to collect" and "you named
// a root that does not exist" are answers someone would act on differently.
func storeGCDirs(dir, layer string) ([]string, error) {
	if dir == "" && layer == "" {
		return nil, nil
	}
	dirs := config.GetImageSearchPaths()
	if layer != "" {
		selected, err := config.ParseDataLayer(layer)
		if err != nil {
			return nil, err
		}
		dirs = config.FilterDirsByLayer(dirs, selected)
		if len(dirs) == 0 {
			return nil, fmt.Errorf("no image directory belongs to the %s layer", selected)
		}
	}
	if dir != "" {
		var matched []string
		for _, candidate := range dirs {
			if strings.Contains(candidate, dir) {
				matched = append(matched, candidate)
			}
		}
		if len(matched) == 0 {
			return nil, fmt.Errorf("no image directory matches %q", dir)
		}
		dirs = matched
	}
	return dirs, nil
}

// printGCReport shows what was found, largest first.
//
// Retained entries are printed too, and not only as a count: "0 collectable"
// with no explanation is indistinguishable from a broken scan, and the reasons
// are what tell a locked artifact from a pinned one from a young one.
func printGCReport(report *store.GCReport, jsonOut bool) error {
	if jsonOut {
		return printJSON(report)
	}
	verb := "collectable"
	if report.Applied {
		verb = "removed"
	}
	for _, entry := range report.Collectable {
		label := entry.Name
		if entry.Staging {
			label = "abandoned staging file"
		}
		fmt.Fprintf(os.Stdout, "  %-10s %-28s %8s  %s (%s)\n", verb, label,
			utils.FormatSize(entry.Size), utils.StylePath(entry.Path),
			fmt.Sprintf("%s %s", entry.Basis, utils.FormatDuration(entry.Age)))
	}
	for _, entry := range report.Retained {
		utils.PrintDebug("kept %s: %s", entry.Path, entry.Reason)
	}
	if len(report.Collectable) == 0 {
		utils.PrintSuccess("Nothing to collect (%d entries kept, grace %s).",
			len(report.Retained), utils.FormatDuration(report.Grace))
		return nil
	}
	if report.Applied {
		utils.PrintSuccess("Removed %d entrie(s), %s freed.",
			len(report.Collectable), utils.FormatSize(report.Reclaimable))
		return nil
	}
	utils.PrintSuccess("%d entrie(s) collectable, %s reclaimable.",
		len(report.Collectable), utils.FormatSize(report.Reclaimable))
	utils.PrintHint("re-run with --apply and --dir or --layer to delete them")
	return nil
}
