package cmd

import (
	"encoding/json"
	"errors"
	"fmt"
	"os"

	"github.com/Justype/condatainer/internal/store"
	"github.com/spf13/cobra"
)

var storeCmd = &cobra.Command{
	Use:   "store",
	Short: "Inspect immutable identity-addressed overlays",
}

func init() {
	rootCmd.AddCommand(storeCmd)
	storeCmd.AddCommand(newStoreListCmd(), newStorePathCmd(), newStoreValidateCmd(),
		newStoreRmCmd(), newStoreGCCmd())
}

func newStoreListCmd() *cobra.Command {
	var equivalence string
	var jsonOutput bool
	cmd := &cobra.Command{
		Use:          "list [name]",
		Short:        "List verified store entries or equivalent candidates",
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
				report = store.Scan(store.ScanOptions{Name: name, Stored: true})
			}
			store.SortReport(&report)
			return printStoreReport(report, jsonOutput)
		},
	}
	cmd.Flags().StringVar(&equivalence, "equiv", "", "list candidates matching an equivalence key or prefix")
	cmd.Flags().BoolVar(&jsonOutput, "json", false, "print JSON")
	return cmd
}

func newStorePathCmd() *cobra.Command {
	var identity string
	cmd := &cobra.Command{
		Use:          "path <name>",
		Short:        "Resolve an exact verified artifact path",
		Args:         cobra.ExactArgs(1),
		SilenceUsage: true,
		RunE: func(_ *cobra.Command, args []string) error {
			if identity == "" {
				return errors.New("--identity is required")
			}
			query, err := store.ParseIdentityQuery(identity)
			if err != nil {
				return err
			}
			candidate, _, err := store.ResolveIdentity(args[0], query, nil)
			if err != nil {
				return err
			}
			fmt.Fprintln(os.Stdout, candidate.Path)
			return nil
		},
	}
	cmd.Flags().StringVar(&identity, "identity", "", "complete identity or unambiguous digest prefix")
	return cmd
}

func newStoreValidateCmd() *cobra.Command {
	return &cobra.Command{
		Use:          "validate [name]",
		Short:        "Regenerate store keys and validate their filenames",
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

func printStoreReport(report store.Report, jsonOutput bool) error {
	if jsonOutput {
		encoder := json.NewEncoder(os.Stdout)
		encoder.SetIndent("", "  ")
		return encoder.Encode(report)
	}
	for _, candidate := range report.Candidates {
		fmt.Fprintf(os.Stdout, "%s\t%s\t%d\t%s\n", candidate.Name, store.FormatKeyRef(candidate.Identity), candidate.Size, candidate.Path)
	}
	for _, issue := range report.Issues {
		fmt.Fprintf(os.Stderr, "invalid\t%s\t%s\n", issue.Path, issue.Error)
	}
	return nil
}
