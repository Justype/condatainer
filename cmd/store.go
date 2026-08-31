package cmd

import (
	"encoding/json"
	"errors"
	"fmt"
	"os"
	"path/filepath"
	"strings"

	"github.com/Justype/condatainer/internal/artifact/meta"
	"github.com/Justype/condatainer/internal/store"
	"github.com/Justype/condatainer/internal/utils"
	"github.com/spf13/cobra"
)

var storeCmd = &cobra.Command{
	Use:   "store",
	Short: "Inspect overlays kept under an exact build identity",
	Long: `The store holds the overlays that could not take a plain name, because another
build of that name was already installed. Each entry is addressed by its name
plus one exact identity, so several builds of a name live side by side and a
project can pin the one it needs.

Overlays under a plain name are listed by 'condatainer list' and removed by
'condatainer remove'.`,
}

func init() {
	rootCmd.AddCommand(storeCmd)
	storeCmd.AddCommand(newStoreAddCmd(), newStoreListCmd(), newStorePathCmd(), newStoreValidateCmd(),
		newStoreUseCmd(), newStoreRmCmd(), newStoreGCCmd())
}

func newStoreListCmd() *cobra.Command {
	var equivalence string
	var jsonOutput bool
	var all bool
	var detail bool
	cmd := &cobra.Command{
		Use:   "list [name]",
		Short: "List store entries, or the builds that can stand in for one",
		Long: `Lists the store entries by address, in columns. Give a name to list only
its builds.

Entries are grouped by images directory, whose header names the directory and
its layer.

The address is what every other store command accepts, so a line here can be
pasted into 'store use'. An identity shown in grey belongs to the build already
installed under the plain name, where the identity is informational rather than
how you address it.

--detail adds the layout, the size, and each file's path relative to its
directory. --json carries the complete keys, absolute paths and exact sizes.

--all includes the overlay installed under the plain name, which is the one
'list' and 'exec -o' resolve. With --equiv, lists the installed builds that can
stand in for an equivalence key instead.`,
		Example: `  condatainer store list
  condatainer store list star/2.7.11b --all --detail
  condatainer store list star/2.7.11b --equiv 9f2c1ab`,
		Args:         cobra.MaximumNArgs(1),
		SilenceUsage: true,
		RunE: func(_ *cobra.Command, args []string) error {
			name := ""
			if len(args) == 1 {
				name = args[0]
			}
			var report store.Report
			if equivalence != "" {
				if name == "" {
					return errors.New("store list --equiv requires a name")
				}
				query, err := store.ParseIdentityQuery(equivalence)
				if err != nil {
					return err
				}
				report, err = store.Equivalent(name, query, nil)
				if err != nil {
					return err
				}
			} else {
				// Flat is opt-in: `store list` lists the store. But an overlay under
				// the plain name is the one that currently answers to it, so any
				// question about which build a name means needs --all to be
				// answerable at all.
				report = store.Scan(store.ScanOptions{Name: name, Stored: true, Flat: all})
			}
			store.SortReport(&report)
			return printStoreReport(report, jsonOutput, detail)
		},
	}
	cmd.Flags().BoolVar(&all, "all", false, "Also list the overlay installed under the plain name")
	cmd.Flags().BoolVar(&detail, "detail", false, "Show layout, size and path for each entry")
	cmd.Flags().StringVar(&equivalence, "equiv", "", "List the builds matching an equivalence key or prefix")
	cmd.Flags().BoolVar(&jsonOutput, "json", false, "Print JSON")
	return cmd
}

func newStorePathCmd() *cobra.Command {
	var identity string
	cmd := &cobra.Command{
		Use:   "path <name>@<identity>",
		Short: "Print the file path of one exact build",
		Long: `Prints the path of the build with the given identity, for use in a script.

Fails when the name and identity do not select exactly one artifact.`,
		Example: `  condatainer store path star/2.7.11b@9f2c1ab
  condatainer store path star/2.7.11b --identity 9f2c1ab`,
		Args:         cobra.ExactArgs(1),
		SilenceUsage: true,
		RunE: func(_ *cobra.Command, args []string) error {
			name, query, err := splitStoreAddress(args[0], identity)
			if err != nil {
				return err
			}
			candidate, _, err := store.ResolveIdentity(name, query, nil)
			if err != nil {
				return err
			}
			fmt.Fprintln(os.Stdout, candidate.Path)
			return nil
		},
	}
	cmd.Flags().StringVar(&identity, "identity", "", "Complete identity or unambiguous digest prefix")
	return cmd
}

func newStoreValidateCmd() *cobra.Command {
	return &cobra.Command{
		Use:   "validate [name]",
		Short: "Check that every entry still matches the identity it is filed under",
		Long: `Recomputes each entry's keys from the artifact itself, ignoring any cache, and
compares them with the name it is filed under. Give a name to check only its
builds.

Prints one line per entry and fails if any entry does not match.`,
		Args:         cobra.MaximumNArgs(1),
		SilenceUsage: true,
		RunE: func(_ *cobra.Command, args []string) error {
			name := ""
			if len(args) == 1 {
				name = args[0]
			}
			report := store.Scan(store.ScanOptions{Name: name, Stored: true, Uncached: true})
			store.SortReport(&report)
			for _, candidate := range report.Candidates {
				fmt.Fprintf(os.Stdout, "ok\t%s\t%s\n", store.FormatKeyRef(candidate.Identity), candidate.Path)
			}
			for _, issue := range report.Issues {
				fmt.Fprintf(os.Stderr, "invalid\t%s\t%s\n", issue.Path, issue.Error)
			}
			if len(report.Issues) != 0 {
				return fmt.Errorf("store validation failed for %d entries", len(report.Issues))
			}
			return nil
		},
	}
}

// printStoreReport lists entries grouped by images directory, as `list` does.
//
// By default the address is the whole row, laid out in columns: it is what the
// other store commands accept, and everything else is answerable from it. --detail
// adds what an inventory wants — layout, size, and the filename relative to the
// directory heading it. The full scheme-backed key, absolute paths and exact byte
// counts are in --json, where something reads them rather than someone.
func printStoreReport(report store.Report, jsonOutput, detail bool) error {
	if jsonOutput {
		encoder := json.NewEncoder(os.Stdout)
		encoder.SetIndent("", "  ")
		return encoder.Encode(report)
	}

	byRoot := map[string][]store.Candidate{}
	var roots []string
	for _, candidate := range report.Candidates {
		if _, seen := byRoot[candidate.Root]; !seen {
			roots = append(roots, candidate.Root)
		}
		byRoot[candidate.Root] = append(byRoot[candidate.Root], candidate)
	}

	for i, root := range roots {
		if i > 0 {
			fmt.Println()
		}
		fmt.Println(dirHeader(root))
		if !detail {
			plain := make([]string, len(byRoot[root]))
			styled := make([]string, len(byRoot[root]))
			for j, candidate := range byRoot[root] {
				plain[j] = storeAddress(candidate)
				styled[j] = styleStoreAddress(candidate)
			}
			printColumns(plain, styled, 0, terminalWidth())
			continue
		}
		for _, candidate := range byRoot[root] {
			fmt.Fprintf(os.Stdout, "%-34s %-6s %9s  %s\n", storeAddress(candidate),
				candidate.Layout, utils.FormatSize(candidate.Size),
				utils.StylePath(relativeTo(root, candidate.Path)))
		}
	}
	for _, issue := range report.Issues {
		fmt.Fprintf(os.Stderr, "%-34s %s: %s\n", "invalid", issue.Path, issue.Error)
	}
	return nil
}

// styleStoreAddress greys the identity of the build installed under the plain
// name: that one answers to its name, so its identity is informational.
func styleStoreAddress(candidate store.Candidate) string {
	suffix := "@" + shortDigest(candidate.Identity)
	if candidate.Layout == store.LayoutFlat {
		return candidate.Name + utils.StyleDebug(suffix)
	}
	return candidate.Name + suffix
}

// relativeTo renders a path against the directory heading it, falling back to
// the absolute path when it lies elsewhere.
func relativeTo(root, path string) string {
	if rel, err := filepath.Rel(root, path); err == nil && !strings.HasPrefix(rel, "..") {
		return rel
	}
	return path
}

// splitStoreAddress reads "<name>" or "<name>@<identity>", with --identity as
// the alternative spelling of the second half.
//
// Splitting on the first @ is unambiguous: an artifact name never contains one,
// while an identity may contain two (scheme@sha256:...), so everything after the
// first @ belongs to the identity. The combined form exists so a line of
// `store list` output can be pasted straight into `store use`.
func splitStoreAddress(arg, flag string) (string, store.IdentityQuery, error) {
	name, inline, hasInline := strings.Cut(arg, "@")
	switch {
	case hasInline && flag != "" && inline != flag:
		return "", store.IdentityQuery{}, fmt.Errorf(
			"two identities given: %q in the address and %q in --identity", inline, flag)
	case !hasInline:
		inline = flag
	}
	if strings.TrimSpace(inline) == "" {
		return "", store.IdentityQuery{}, errors.New(
			"an identity is required: give it as <name>@<identity>, or with --identity")
	}
	query, err := store.ParseIdentityQuery(inline)
	if err != nil {
		return "", store.IdentityQuery{}, err
	}
	return name, query, nil
}

// storeAddress renders the address that selects one artifact, which is the form
// every store command accepts back.
func storeAddress(candidate store.Candidate) string {
	return candidate.Name + "@" + shortDigest(candidate.Identity)
}

// shortDigest renders the identity prefix that names the file on disk.
func shortDigest(ref meta.KeyRef) string {
	if len(ref.SHA256) < store.DefaultPrefixChars {
		return ref.SHA256
	}
	return ref.SHA256[:store.DefaultPrefixChars]
}
