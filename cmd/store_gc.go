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
		Use:   "rm <name>",
		Short: "Remove one store entry",
		Long: `Delete a single identity from the store.

Only a store entry. A flat artifact answers to a bare name and belongs to
` + "`condatainer remove`" + `; taking a name out of service through a command
that reads as identity housekeeping would be a surprise.

The entry must be free: an exclusive lock fails while a container is reading it,
and a cleared write bit means someone pinned that identity deliberately.`,
		Args:         cobra.ExactArgs(1),
		SilenceUsage: true,
		RunE: func(_ *cobra.Command, args []string) error {
			if identity == "" {
				return fmt.Errorf("--identity is required")
			}
			query, err := store.ParseIdentityQuery(identity)
			if err != nil {
				return err
			}
			removed, err := store.Remove(args[0], query, nil)
			if err != nil {
				return err
			}
			utils.PrintSuccess("Removed %s (%s)", utils.StylePath(removed.Path),
				store.FormatKeyRef(removed.Identity))
			return nil
		},
	}
	cmd.Flags().StringVar(&identity, "identity", "", "complete identity or unambiguous digest prefix")
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
		Long: `Report which store entries could be deleted, and how much that would free.

Reporting is the default because no timestamp separates "installed and never
opened" from "opened last week by a job about to run again", and a store entry
costs hours to reproduce. --apply repeats the whole judgement and deletes what
still passes.

Age is the newest of atime, mtime and ctime, and the report names which one
decided each entry. atime is the signal that tracks use; where a mount disables
it the other two answer, which is why every row saying "ctime" means the
filesystem is noatime rather than the artifacts being cold.

-D/--dir and -l/--layer scope the run, exactly as they do for ` + "`list`" + `.
--apply requires one of them: deleting across every writable tier at once,
including a group root shared by people who are not at the keyboard, is not
something to do by omission.`,
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
	cmd.Flags().BoolVar(&apply, "apply", false, "delete the entries that still pass; requires --dir or --layer")
	cmd.Flags().IntVar(&graceDay, "grace", 0, "age in days below which nothing is collectable (default: store_gc_grace)")
	cmd.Flags().StringVarP(&dir, "dir", "D", "", "limit to image directories matching this substring")
	cmd.Flags().StringVarP(&layer, "layer", "l", "", "limit to a data layer: u/user, r/app-root, e/extra-root")
	cmd.Flags().BoolVar(&jsonOut, "json", false, "print JSON")
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
