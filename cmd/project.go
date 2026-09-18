package cmd

import (
	"encoding/json"
	"errors"
	"fmt"
	"os"
	"sort"
	"strings"

	"github.com/Justype/condatainer/internal/config"
	"github.com/Justype/condatainer/internal/logging"
	"github.com/Justype/condatainer/internal/project"
	"github.com/Justype/condatainer/internal/project/lock"
	"github.com/Justype/condatainer/internal/project/orchestrate"
	"github.com/Justype/condatainer/internal/project/publish"
	"github.com/Justype/condatainer/internal/project/restore"
	"github.com/Justype/condatainer/internal/registry"
	"github.com/Justype/condatainer/internal/utils"
	"github.com/spf13/cobra"
)

// projectDir is the --project flag, shared by every subcommand.
var projectDir string

var projectCmd = &cobra.Command{
	Use:   "project",
	Short: "Pin and restore a project's artifacts",
	Long: `A project is a folder holding a cnt-lock/ directory. Locking pins every #DEP:
in its scripts to one exact artifact and records how to rebuild it, so the same
scripts run the same way on another machine.

Every artifact carries two keys:
- identity    names one exact build;
- equivalence is shared by any build that can stand in for it.

A restore accepts an equivalent build unless you ask for the identity.`,
}

func init() {
	rootCmd.AddCommand(projectCmd)
	projectCmd.PersistentFlags().StringVar(&projectDir, "project", "", "Project root (default: the current directory)")
	projectCmd.AddCommand(newProjectLockCmd(), newProjectStatusCmd(), newProjectPinCmd(), newProjectUnpinCmd(),
		newProjectSelectDistroCmd(), newProjectListCmd(), newProjectValidateCmd(),
		newProjectRestoreCmd(), newProjectRegistryCmd(), newProjectPushCmd())
}

// projectRoot resolves the root a command acts on — the current directory, or
// the one --project names — and reports it once.
func projectRoot(announce bool) (string, error) {
	cwd, err := os.Getwd()
	if err != nil {
		return "", err
	}
	root, err := lock.RootFor(projectDir, cwd)
	if errors.Is(err, lock.ErrNoProject) {
		return "", fmt.Errorf("%w\na project is a directory containing %s; create one with %s, or name it with %s",
			err, "cnt-lock/", "condatainer project lock", "--project DIR")
	}
	if err != nil {
		return "", err
	}
	if announce {
		utils.PrintMessage("Project: %s", utils.StylePath(root))
	}
	return root, nil
}

// projectRootOrInit resolves the root like projectRoot, but treats a directory
// with no cnt-lock/ as one to create rather than an error.
//
// Only 'project lock' uses it: scanning a directory for declarations is what
// makes it a project, so there is nothing for it to refuse. Every other
// subcommand acts on pins that must already exist, and would be creating
// an empty project as a side effect of a mistyped path.
func projectRootOrInit(announce bool) (string, error) {
	root, err := projectRoot(false)
	if err != nil && !errors.Is(err, lock.ErrNoProject) {
		return "", err
	}
	created := err != nil
	if created {
		if root, err = os.Getwd(); err != nil {
			return "", err
		}
	}
	if announce {
		if created {
			utils.PrintMessage("New project: %s", utils.StylePath(root))
		} else {
			utils.PrintMessage("Project: %s", utils.StylePath(root))
		}
	}
	return root, nil
}

func newProjectLockCmd() *cobra.Command {
	var jsonOutput bool
	cmd := &cobra.Command{
		Use:   "lock",
		Short: "Create or update the lock from #DEP: declarations",
		Long: `Reads the .sh and .bash scripts in the project and pins the overlay each
#DEP: names. Creates cnt-lock/ when there is none.

- Pins every declaration, and drops entries nothing declares any more.
- Skips cnt-lock/, overlays/ and dot directories, and never follows symlinks.

A declaration whose overlay is not installed fails the lock. Other installed
builds of the same name are listed rather than pinned; pick one of those with
'condatainer project pin'.`,
		Args:         cobra.NoArgs,
		SilenceUsage: true,
		RunE: func(cmd *cobra.Command, _ []string) error {
			root, err := projectRootOrInit(!jsonOutput)
			if err != nil {
				return err
			}
			result, err := orchestrate.Lock(cmd.Context(), root, orchestrate.LockOptions{})
			if err != nil {
				return err
			}
			return reportLockState(result, jsonOutput)
		},
	}
	cmd.Flags().BoolVar(&jsonOutput, "json", false, "Print JSON")
	return cmd
}

func newProjectStatusCmd() *cobra.Command {
	var jsonOutput bool
	cmd := &cobra.Command{
		Use:   "status",
		Short: "Show whether this directory is a project, and what it pins",
		Long: `Reports whether the current directory (or --project DIR) is a project,
how many artifacts it pins, and which #DEP: declarations have no pin yet.

Read-only: it never scans-and-writes the way 'condatainer project lock' does.`,
		Args:         cobra.NoArgs,
		SilenceUsage: true,
		RunE: func(_ *cobra.Command, _ []string) error {
			// Resolved the same non-walking way every other subcommand in
			// this family resolves --project/ambient cwd (projectRoot →
			// lock.RootFor → the strict lock.RootAt) — never project.Status's
			// own StandingAt walk-up, which exists for an ambient caller with
			// no --project concept at all (a helper's #REQUIRED_OVERLAYS:
			// resolution, the dashboard's raw cwd field), not for this
			// command. See internal/project/README.md, "Acting in a
			// project": reaching an ancestor project from a --project
			// directory that looks empty is a hazard an ambient read never
			// is, and every sibling subcommand already refuses it.
			root, err := projectRoot(false)
			if err != nil && !errors.Is(err, lock.ErrNoProject) {
				return err
			}
			if err != nil {
				return printNoProjectStatus(jsonOutput)
			}
			// projectRoot's --project branch takes the directory as named
			// without verifying cnt-lock/ is actually there (RootFor's own
			// contract); StatusAt is what actually checks, so a --project
			// pointing at a bare subdirectory still reports "no project"
			// instead of walking up to an unrelated ancestor.
			status, err := project.StatusAt(root)
			if errors.Is(err, lock.ErrNoProject) {
				return printNoProjectStatus(jsonOutput)
			}
			if err != nil {
				return err
			}
			if jsonOutput {
				return printJSON(status)
			}
			utils.PrintMessage("Project: %s", utils.StylePath(status.Root))
			utils.PrintMessage("%d pin(s)", status.PinCount)
			if len(status.Unpinned) > 0 {
				utils.PrintWarning("%d declaration(s) not pinned:", len(status.Unpinned))
				for _, key := range status.Unpinned {
					utils.PrintMessage("  %s", key)
				}
				utils.PrintHint("Run %s to pin them.", utils.StyleAction("condatainer project lock"))
			}
			printUnpinnedHelperOverlays(status.UnpinnedHelperOverlays)
			printManualPinUsage(status.ManualPinUsage)
			return nil
		},
	}
	cmd.Flags().BoolVar(&jsonOutput, "json", false, "Print JSON")
	return cmd
}

// printNoProjectStatus reports "no project" for `project status`, in
// whichever format was asked for — shared by projectRoot's ambient
// ErrNoProject and StatusAt's --project-directory ErrNoProject, so the two
// ways of finding out there is no project here read identically.
func printNoProjectStatus(jsonOutput bool) error {
	if jsonOutput {
		return printJSON(&project.ProjectStatus{})
	}
	utils.PrintMessage("No project here.")
	utils.PrintHint("Run %s to create one.", utils.StyleAction("condatainer project lock"))
	return nil
}

// printUnpinnedHelperOverlays prints the "used by helpers but not pinned"
// section §1/§2 add to both `project status` and `project lock`'s report —
// its own heading, kept out of the #DEP:-derived unpinned list so the two
// are never read as carrying the same weight.
func printUnpinnedHelperOverlays(names []string) {
	if len(names) == 0 {
		return
	}
	utils.PrintMessage("Overlays used by helpers but not pinned:")
	for _, name := range names {
		utils.PrintMessage("  %s", name)
	}
	utils.PrintHint("Run %s to pin one, if it should be locked.",
		utils.StyleAction("condatainer project pin <name> <identity>"))
}

// printManualPinUsage prints, next to the manual pins report, which
// helper(s) are recorded using each one — a fact, never a suggestion to
// unpin. Shared between `project status` and `project lock`'s report.
func printManualPinUsage(usage map[string][]string) {
	if len(usage) == 0 {
		return
	}
	keys := make([]string, 0, len(usage))
	for key := range usage {
		keys = append(keys, key)
	}
	sort.Strings(keys)
	utils.PrintMessage("Manual pin usage:")
	for _, key := range keys {
		names := usage[key]
		if len(names) == 0 {
			utils.PrintMessage("  %s: no recorded helper usage", key)
			continue
		}
		utils.PrintMessage("  %s: used by %s", key, strings.Join(names, ", "))
	}
}

func newProjectPinCmd() *cobra.Command {
	cmd := &cobra.Command{
		Use:   "pin <request> [identity]",
		Short: "Pin one declaration to an exact artifact",
		Long: `Pins one declaration to an exact artifact and copies what is needed to
rebuild it into cnt-lock/.

- name/version    pin it to one installed overlay
- a project .sqf  re-read that file and record its current identity

An identity is a full scheme@sha256:<hex> key or any unambiguous prefix of it.
Give one to choose between installed copies of the same name; a project path
takes none, since it already names the file it means.`,
		Example: `  condatainer project pin star/2.7.11b a31f902c12ab
  condatainer project pin overlays/combined.sqf`,
		Args:         cobra.RangeArgs(1, 2),
		SilenceUsage: true,
		RunE: func(cmd *cobra.Command, args []string) error {
			root, err := projectRoot(true)
			if err != nil {
				return err
			}
			current, err := lock.Load(root)
			if err != nil {
				return err
			}
			var identity string
			if len(args) > 1 {
				identity = args[1]
			}
			pinned, err := lock.Pin(root, args[0], identity, lock.PinOptions{})
			if err != nil {
				return err
			}
			pinned.Manual = !declaredByAScript(root, pinned.Request)
			upstream := orchestrate.RecordUpstream(cmd.Context(), root, current, pinned)
			if err := lock.Apply(root, current, pinned); err != nil {
				return err
			}

			utils.PrintSuccess("Pinned %s", utils.StyleName(pinned.Name))
			utils.PrintMessage("  identity %s", pinned.Identity.Digest())
			if pinned.Manual {
				utils.PrintMessage("  no script declares it, so it is recorded as a manual pin")
			}
			utils.PrintMessage("  read from %s", utils.StylePath(pinned.Path))
			for _, artifact := range pinned.Vendored {
				line := "  vendored " + artifact
				if remote, ok := upstream[artifact]; ok {
					line += " (fetchable from " + remote + ")"
				}
				utils.PrintMessage("%s", line)
			}
			printUnpublishedNote(current, []*lock.Pinned{pinned})
			return nil
		},
	}
	return cmd
}

func newProjectSelectDistroCmd() *cobra.Command {
	var auto bool
	cmd := &cobra.Command{
		Use:   "select-distro [distro]",
		Short: "Choose which distro's base this project restores to",
		Long: `Pins <distro>/base as this project's root, overriding what
'condatainer project lock' would otherwise derive from the configured
default_distro.

--auto clears the override and re-derives from the configured default_distro
immediately, rather than waiting for the next 'condatainer project lock'.`,
		Example: `  condatainer project select-distro rocky9
  condatainer project select-distro --auto`,
		Args:         cobra.MaximumNArgs(1),
		SilenceUsage: true,
		RunE: func(cmd *cobra.Command, args []string) error {
			if auto == (len(args) == 1) {
				return fmt.Errorf("give exactly one of a distro or --auto")
			}
			root, err := projectRoot(true)
			if err != nil {
				return err
			}
			current, err := lock.Load(root)
			if err != nil {
				return err
			}
			var pinned *lock.Pinned
			if auto {
				pinned, err = lock.AutoDistro(root, current, config.ResolvedDefaultDistro(), lock.PinOptions{})
			} else {
				pinned, err = lock.SelectDistro(root, current, args[0], lock.PinOptions{})
			}
			if err != nil {
				return err
			}
			upstream := orchestrate.RecordUpstream(cmd.Context(), root, current, pinned)
			if err := lock.Publish(root, current); err != nil {
				return err
			}

			if auto {
				utils.PrintSuccess("Root reset to the configured default: %s", utils.StyleName(pinned.Name))
			} else {
				utils.PrintSuccess("Root pinned to %s", utils.StyleName(pinned.Name))
			}
			utils.PrintMessage("  identity %s", pinned.Identity.Digest())
			utils.PrintMessage("  read from %s", utils.StylePath(pinned.Path))
			for _, artifact := range pinned.Vendored {
				line := "  vendored " + artifact
				if remote, ok := upstream[artifact]; ok {
					line += " (fetchable from " + remote + ")"
				}
				utils.PrintMessage("%s", line)
			}
			return nil
		},
	}
	cmd.Flags().BoolVar(&auto, "auto", false, "Clear a manual override and re-derive from the configured default_distro")
	return cmd
}

func newProjectUnpinCmd() *cobra.Command {
	cmd := &cobra.Command{
		Use:   "unpin <request>",
		Short: "Remove a manual pin",
		Long: `Removes a pin the project recorded itself, and the provenance nothing else
reaches once it is gone.

Only a manual pin — one no script declares. A pin a declaration produced is
removed by removing the declaration and running 'condatainer project lock',
which drops every pin nothing asks for any more.

The image itself is never deleted. A project path keeps its file; a named pin
keeps its overlay in the images root.`,
		Example: `  condatainer project unpin overlays/env.sqf
  condatainer project unpin ubuntu24/build-essential`,
		Args:         cobra.ExactArgs(1),
		SilenceUsage: true,
		RunE: func(cmd *cobra.Command, args []string) error {
			root, err := projectRoot(true)
			if err != nil {
				return err
			}
			l, err := lock.Load(root)
			if err != nil {
				return err
			}
			request, err := lock.CanonicalRequest(args[0])
			if err != nil {
				return err
			}
			entry, pinned := l.Pins[request]
			if !pinned {
				return fmt.Errorf("nothing pins %s", request)
			}
			if !entry.Manual {
				return fmt.Errorf("%s is pinned from a project declaration; remove the declaration and run `condatainer project lock`", request)
			}

			// Reachability before and after says what the removal orphans.
			// Publish prunes exactly that, so this reports rather than decides.
			before, _ := lock.Verify(root, l)
			delete(l.Pins, request)
			after, _ := lock.Verify(root, l)
			if err := lock.Publish(root, l); err != nil {
				return err
			}

			utils.PrintSuccess("Unpinned %s", utils.StyleName(request))
			for _, artifact := range orphaned(before, after) {
				utils.PrintMessage("  pruned %s", artifact)
			}
			if destination, ok := strings.CutPrefix(request, lock.PathPrefix); ok {
				utils.PrintMessage("  %s is left where it is", utils.StylePath(destination))
			}
			return nil
		},
	}
	return cmd
}

// orphaned lists the readable artifacts a removal left unreachable, sorted.
// Unreadable ones are left out because Prune leaves them on disk.
func orphaned(before, after *lock.Verified) []string {
	var out []string
	for artifact := range before.Reachable {
		if after.Reachable[artifact] {
			continue
		}
		if _, readable := before.Entries[artifact]; readable {
			out = append(out, artifact)
		}
	}
	sort.Strings(out)
	return out
}

func newProjectListCmd() *cobra.Command {
	var jsonOutput bool
	cmd := &cobra.Command{
		Use:   "list",
		Short: "List this project's pins",
		Long: `Lists every pin in cnt-lock/lock.json and the artifact it names.

Works in a fresh clone that has restored nothing. To ask whether an artifact is
actually available, use 'condatainer project restore --dry-run'.

A pin marked 'manual' is one the project recorded itself: no script declares it,
'project lock' will not sweep it, and 'project unpin' removes it.`,
		Args:         cobra.NoArgs,
		SilenceUsage: true,
		RunE: func(_ *cobra.Command, _ []string) error {
			root, err := projectRoot(!jsonOutput)
			if err != nil {
				return err
			}
			l, err := lock.Load(root)
			if err != nil {
				return err
			}
			// Problems are for `project validate` to report; here an entry that
			// does not verify is one line that says so, so a broken pin is still
			// listed with the key that addresses it.
			verified, _ := lock.Verify(root, l)
			pins := pinReports(l, verified)

			if jsonOutput {
				return printJSON(struct {
					Root string      `json:"root"`
					Pins []pinReport `json:"pins"`
				}{Root: root, Pins: pins})
			}
			if len(pins) == 0 {
				utils.PrintMessage("This project pins nothing.")
				utils.PrintHint("Run %s to pin what its scripts declare.",
					utils.StyleAction("condatainer project lock"))
				return nil
			}
			// Rows go to stdout so the listing can be piped; the summary and
			// any warning stay on the console, which is stderr.
			width := 0
			for _, pin := range pins {
				width = max(width, len(pinLabel(pin)))
			}
			var manual int
			for _, pin := range pins {
				label := pinLabel(pin)
				// Padded from the unstyled label: styling adds escape bytes that
				// a width verb would count as columns.
				row := "  " + utils.StyleName(label) + strings.Repeat(" ", width-len(label))
				switch {
				case pin.Problem != "":
					row += "  " + utils.StyleWarning("unreadable")
				case pin.Manual:
					manual++
					row += "  " + short(pin.Identity) + "  " + utils.StyleWarning("manual")
				default:
					row += "  " + short(pin.Identity)
				}
				fmt.Println(row)
				if pin.Problem != "" {
					utils.PrintWarning("  %s", pin.Problem)
				}
			}
			if manual > 0 {
				utils.PrintMessage("%d pin(s), %d manual", len(pins), manual)
				return nil
			}
			utils.PrintMessage("%d pin(s)", len(pins))
			return nil
		},
	}
	cmd.Flags().BoolVar(&jsonOutput, "json", false, "Print JSON")
	return cmd
}

// pinLabel is the pin key, and the artifact name after it when the key does not
// already carry one — a path addresses a file, which says nothing about what is
// in it, while a name request is the name.
func pinLabel(pin pinReport) string {
	if pin.Request == lock.BaseKey {
		if pin.Name == "" {
			return "base"
		}
		return "base (" + pin.Name + ")"
	}
	if pin.Name == "" || pin.Name == pin.Request {
		return pin.Request
	}
	return pin.Request + " (" + pin.Name + ")"
}

// pinReport is one pin as the listing renders it. Name and Identity come from
// the vendored records rather than from the key, which addresses an artifact
// without describing it.
type pinReport struct {
	Request  string `json:"request"`
	Artifact string `json:"artifact"`
	Name     string `json:"name,omitempty"`
	Identity string `json:"identity,omitempty"`
	Manual   bool   `json:"manual,omitempty"`
	// Problem is why the artifact could not be read, empty when it verified.
	Problem string `json:"problem,omitempty"`
}

func pinReports(l *lock.Lock, verified *lock.Verified) []pinReport {
	out := make([]pinReport, 0, len(l.Pins))
	for _, request := range l.Requests() {
		pin := l.Pins[request]
		report := pinReport{Request: request, Artifact: pin.Artifact, Manual: pin.Manual}
		if entry, ok := verified.Entries[pin.Artifact]; ok {
			report.Name, report.Identity = entry.Manifest.Name, entry.Identity.Digest()
		} else {
			report.Problem = pin.Artifact + " is missing or does not verify"
		}
		out = append(out, report)
	}
	return out
}

// declaredByAScript reports whether anything in the project asks for this
// request, which decides whether the pin is manual.
//
// Derived once, at pin time, and then recorded: `project lock` sweeps every pin
// the scan no longer produces, and a pin nothing ever declared has to survive
// that. Artifacts reach a project without a #DEP: routinely — a helper names
// them in #REQUIRED_OVERLAYS:, and a frozen environment is named by nobody.
//
// A scan that cannot run answers "declared", which keeps a pin sweepable rather
// than silently permanent: the cautious answer is the one a later lock can still
// correct.
func declaredByAScript(root, request string) bool {
	scanned, err := lock.Scan(root, lock.ScanOptions{})
	if err != nil {
		return true
	}
	for _, candidate := range scanned.Requests {
		if candidate.Key == request {
			return true
		}
	}
	return false
}

// printUnpublishedNote prints the same note `project lock`'s old
// noteUnpublished did, once, when pinned left any frozen-environment
// artifact without a recorded remote — the one case a pin cannot be rebuilt
// from what it vendored, so a registry copy is the only route another
// checkout has to obtain it.
func printUnpublishedNote(l *lock.Lock, pinned []*lock.Pinned) {
	if len(project.UnpublishedFrozenEnv(l, pinned)) == 0 {
		return
	}
	utils.PrintNote("A frozen environment cannot be rebuilt; run %s so another checkout can restore this project",
		utils.StyleAction("condatainer project push"))
}

func newProjectRegistryCmd() *cobra.Command {
	cmd := &cobra.Command{
		Use:   "registry",
		Short: "Manage where this project publishes its artifacts",
		Long: `Shows where 'condatainer project push' publishes, as recorded in
cnt-lock/lock.json.

Record one with 'set', forget it with 'unset'.`,
		Args:         cobra.NoArgs,
		SilenceUsage: true,
		RunE: func(_ *cobra.Command, _ []string) error {
			root, err := projectRoot(false)
			if err != nil {
				return err
			}
			l, err := lock.Load(root)
			if err != nil {
				return err
			}
			if l.OCI.Empty() {
				utils.PrintMessage("No registry is recorded. Set one with:")
				utils.PrintMessage("  condatainer project registry set <registry>/<owner>/<repo>")
				return nil
			}
			utils.PrintMessage("push     %s", utils.StylePath(l.OCI.Push))
			utils.PrintMessage("audience %s", audienceOrDefault(l.OCI.Audience))
			if l.Source != "" {
				utils.PrintMessage("source   %s", l.Source)
			}
			return nil
		},
	}
	cmd.AddCommand(newProjectRegistrySetCmd(), newProjectRegistryUnsetCmd())
	return cmd
}

func newProjectRegistrySetCmd() *cobra.Command {
	var audience, source string
	cmd := &cobra.Command{
		Use:   "set <registry>/<owner>/<repo>",
		Short: "Record where this project publishes its artifacts",
		Long: `Records where 'condatainer project push' publishes, in cnt-lock/lock.json.

Give a repository with no tag or digest. Every artifact goes into that one
repository, with its name carried in the tag, so a single package covers the
whole project.

- --audience says who can pull from the registry. A public one refuses an app
  whose recipe does not declare that it may be redistributed; a restricted one
  takes it.
- --source is the code repository the published packages link back to.`,
		Example:      `  condatainer project registry set ghcr.io/my-lab/rnaseq-2026/cnt --audience restricted`,
		Args:         cobra.ExactArgs(1),
		SilenceUsage: true,
		RunE: func(cmd *cobra.Command, args []string) error {
			root, err := projectRoot(true)
			if err != nil {
				return err
			}
			l, err := lock.Load(root)
			if err != nil {
				return err
			}
			previous := l.OCI
			l.OCI = lock.OCI{
				Push:     registry.TrimBaseScheme(args[0]),
				Audience: strings.ToLower(strings.TrimSpace(audience)),
			}
			if cmd.Flags().Changed("source") {
				l.Source = strings.TrimSpace(source)
			} else if l.Source == "" {
				// Derived once and recorded, never re-derived at push time: two
				// collaborators with different remote spellings would otherwise
				// publish two different annotations for one project.
				l.Source = publish.OriginURL(cmd.Context(), root)
			}
			if err := lock.Publish(root, l); err != nil {
				l.OCI = previous
				return err
			}

			utils.PrintSuccess("Publishing to %s", utils.StylePath(l.OCI.Push))
			utils.PrintMessage("  audience %s", audienceOrDefault(l.OCI.Audience))
			if l.Source != "" {
				utils.PrintMessage("  source   %s", l.Source)
			} else {
				utils.PrintMessage("  no repository URL was derived from origin; pass --source to record one")
			}
			return nil
		},
	}
	cmd.Flags().StringVar(&audience, "audience", "public", "Who can pull from this registry: public or restricted")
	cmd.Flags().StringVar(&source, "source", "", "Code repository the packages link back to (default: the origin remote)")
	return cmd
}

func newProjectRegistryUnsetCmd() *cobra.Command {
	return &cobra.Command{
		Use:   "unset",
		Short: "Forget where this project publishes",
		Long: `Forgets the recorded destination, so 'condatainer project push' has nowhere
to publish until one is set again.

Where the artifacts were already published is untouched, so a restore can still
download them.`,
		Args:         cobra.NoArgs,
		SilenceUsage: true,
		RunE: func(_ *cobra.Command, _ []string) error {
			root, err := projectRoot(true)
			if err != nil {
				return err
			}
			l, err := lock.Load(root)
			if err != nil {
				return err
			}
			if l.OCI.Empty() {
				utils.PrintMessage("No registry was recorded.")
				return nil
			}
			was := l.OCI.Push
			l.OCI = lock.OCI{}
			if err := lock.Publish(root, l); err != nil {
				return err
			}
			utils.PrintSuccess("No longer publishing to %s", utils.StylePath(was))
			return nil
		},
	}
}

// audienceOrDefault renders what an unrecorded audience means, rather than
// printing nothing and leaving the reader to guess which way it falls.
func audienceOrDefault(v string) string {
	if strings.TrimSpace(v) == "" {
		return "public (default)"
	}
	return v
}

func newProjectPushCmd() *cobra.Command {
	var all, closure, dryRun, jsonOutput bool
	var repository string
	cmd := &cobra.Command{
		Use:   "push",
		Short: "Publish this project's artifacts and record where they landed",
		Long: `Publishes the locked artifacts to the project's registry and records where
each landed, so a later restore downloads them instead of rebuilding.

- The project must pass 'condatainer project validate'; nothing is uploaded
  otherwise. --dry-run reports what is wrong and still shows the plan.
- Each artifact is tagged by its name and by its identity.
- Publishes only what is already built; 'project restore' produces the rest.
- Skips pins a recipe source already serves, since a restore tries those first.
  --all uploads them too.`,
		Example: `  condatainer project push --dry-run
  condatainer project push
  condatainer project push --all --closure`,
		Args:         cobra.NoArgs,
		SilenceUsage: true,
		RunE: func(cmd *cobra.Command, _ []string) error {
			root, err := projectRoot(true)
			if err != nil {
				return err
			}
			l, err := lock.Load(root)
			if err != nil {
				return err
			}
			verified, problems := lock.Verify(root, l)
			scanned, err := lock.Scan(root, lock.ScanOptions{})
			if err != nil {
				return err
			}
			// Push holds a project to what `project validate` asks, not to the
			// lock's internal consistency alone: what is published has to be what
			// the checkout describes, and a declaration the lock cannot reproduce
			// is missing from it. A dry run reports these and still shows the plan.
			invalid := projectProblems(l, scanned, problems)
			if len(invalid) > 0 && !dryRun {
				return fmt.Errorf("the project is not valid, so nothing was published:\n  %s",
					strings.Join(invalid, "\n  "))
			}

			opts := publish.Options{All: all, Closure: closure, Repository: repository}
			plan, err := publish.Build(cmd.Context(), root, l, verified, opts)
			if err != nil {
				return err
			}
			cat, err := config.OpenCatalog(cmd.Context())
			if err != nil {
				logging.FromContext(cmd.Context()).Debug("no catalog, so upstream copies were not checked", "err", err)
			}
			publish.Refine(cmd.Context(), root, plan, cat, opts)

			if dryRun {
				plan.Problems = append(plan.Problems, invalid...)
				return reportPushPlan(plan, jsonOutput)
			}
			if !plan.Complete() {
				return fmt.Errorf("nothing was published:\n  %s", strings.Join(plan.Problems, "\n  "))
			}
			report, runErr := publish.Run(cmd.Context(), root, plan, opts)
			if err := reportPush(report, jsonOutput); err != nil {
				return err
			}
			return runErr
		},
	}
	cmd.Flags().BoolVar(&all, "all", false, "Publish every pin, including ones a recipe source already serves")
	cmd.Flags().BoolVar(&closure, "closure", false, "Also publish build dependencies")
	cmd.Flags().StringVar(&repository, "registry", "", "Publish to this repository instead of the recorded one")
	cmd.Flags().BoolVar(&dryRun, "dry-run", false, "Report what would be published and upload nothing")
	cmd.Flags().BoolVar(&jsonOutput, "json", false, "Print JSON")
	return cmd
}

// projectProblems is everything `project validate` reports, as lines: a lock
// inconsistent with what it vendors, a declaration nothing can pin, and a
// declaration with no pin. Push and validate share it so a project cannot be
// publishable and invalid at once.
func projectProblems(l *lock.Lock, scanned *lock.ScanResult, problems []lock.Problem) []string {
	out := problemStrings(problems)
	if reason := missingBaseProblem(l); reason != "" {
		out = append(out, reason)
	}
	for _, finding := range scanned.Findings {
		out = append(out, fmt.Sprintf("%s:%d: %s", finding.Script, finding.Line, finding.Reason))
	}
	for _, request := range requestsNeedingPin(l, scanned) {
		out = append(out, fmt.Sprintf("%s is declared but not pinned (%s)",
			request.Key, strings.Join(request.Scripts, ", ")))
	}
	return out
}

// missingBaseProblem is the one line `project validate`/`push` add when no
// root is pinned, or "" when one is. `project lock` pins one unconditionally,
// so the only way to reach this is a lock written before this feature
// existed, or hand-edited.
func missingBaseProblem(l *lock.Lock) string {
	if _, ok := l.Pins[lock.BaseKey]; ok {
		return ""
	}
	return "no root is pinned; run `condatainer project lock`"
}

func reportPushPlan(plan *publish.Plan, jsonOutput bool) error {
	if jsonOutput {
		return printJSON(plan)
	}
	utils.PrintMessage("Publishing to %s (%s)", utils.StylePath(plan.Repository), audienceOrDefault(plan.Audience))
	if plan.Source != "" {
		utils.PrintMessage("  packages link to %s", plan.Source)
	}
	for _, step := range plan.Steps {
		switch step.Disposition {
		case publish.Upload:
			utils.PrintMessage("  upload   %s  %s", utils.StyleName(step.Name), strings.Join(step.Tags, " "))
		case publish.Present:
			utils.PrintMessage("  present  %s  already published", utils.StyleName(step.Name))
		case publish.Served:
			utils.PrintMessage("  upstream %s  %s", utils.StyleName(step.Name), step.Remote.Repository)
		case publish.Refused:
			utils.PrintWarning("  refused  %s  %s", utils.StyleName(step.Name), step.Reason)
		}
		// Reported, never judged: no allowlist of channels could be maintained
		// honestly, and interpreting a few hundred licence strings is what
		// #LICENSE: exists not to do. The operator sees what went into the solve.
		if len(step.Channels) > 0 {
			utils.PrintMessage("           channels: %s", strings.Join(step.Channels, ", "))
		}
	}
	for _, name := range plan.Ambiguous {
		utils.PrintWarning("  two pins are named %s, so neither takes the plain tag", utils.StyleName(name))
	}
	utils.PrintMessage("%d to upload", plan.Uploads())
	for _, problem := range plan.Problems {
		utils.PrintError("%s", problem)
	}
	if !plan.Complete() {
		return fmt.Errorf("this project cannot be published as it stands")
	}
	return nil
}

func reportPush(report *publish.Report, jsonOutput bool) error {
	if jsonOutput {
		return printJSON(report)
	}
	for _, step := range report.Published {
		utils.PrintMessage("  %-8s %s", step.Disposition, utils.StyleName(step.Name))
	}
	for _, name := range report.AmbiguousNames {
		utils.PrintWarning("  two pins are named %s, so neither takes the plain tag", utils.StyleName(name))
	}
	for _, failure := range report.Failures {
		utils.PrintError("  %s", failure)
	}
	if len(report.Failures) == 0 {
		utils.PrintSuccess("Published %d artifact(s) to %s", len(report.Published), utils.StylePath(report.Repository))
	}
	return nil
}

func problemStrings(problems []lock.Problem) []string {
	out := make([]string, 0, len(problems))
	for _, problem := range problems {
		out = append(out, problem.String())
	}
	return out
}

func newProjectValidateCmd() *cobra.Command {
	var jsonOutput bool
	cmd := &cobra.Command{
		Use:   "validate",
		Short: "Check the lock is complete and consistent with the scripts",
		Long: `Checks the lock and the files it carries:

- every declaration in the scripts has a pin;
- every pin points at a recorded artifact;
- every recorded artifact still regenerates the keys it claims;
- every dependency resolves to another recorded artifact.

It looks at the lock only, never at what is installed here. To ask whether a
script can actually run on this machine, use 'condatainer check <script>'.`,
		Args:         cobra.NoArgs,
		SilenceUsage: true,
		RunE: func(_ *cobra.Command, _ []string) error {
			root, err := projectRoot(!jsonOutput)
			if err != nil {
				return err
			}
			current, err := lock.Load(root)
			if err != nil {
				return err
			}
			scanned, err := lock.Scan(root, lock.ScanOptions{})
			if err != nil {
				return err
			}
			_, problems := lock.Verify(root, current)
			needPin := requestsNeedingPin(current, scanned)

			total := len(problems) + len(scanned.Findings) + len(needPin)
			if missingBaseProblem(current) != "" {
				total++
			}
			if jsonOutput {
				if err := printJSON(validateReport(root, current, problems, needPin, scanned)); err != nil {
					return err
				}
				// The exit status is the answer; the format only changes how it
				// is explained.
				if total > 0 {
					return fmt.Errorf("project is not valid: %d problem(s)", total)
				}
				return nil
			}
			for _, problem := range projectProblems(current, scanned, problems) {
				utils.PrintError("%s", problem)
			}
			if total > 0 {
				utils.PrintHint("Run %s to choose an artifact for each unpinned request.",
					utils.StyleAction("condatainer project pin <request> <identity>"))
				return fmt.Errorf("project is not valid: %d problem(s)", total)
			}
			utils.PrintSuccess("Project is valid: %d pin(s), %d script(s)",
				len(current.Pins), len(scanned.Scripts))
			return nil
		},
	}
	cmd.Flags().BoolVar(&jsonOutput, "json", false, "Print JSON")
	return cmd
}

// requestsNeedingPin reports declarations with no pin, without mutating
// the lock the way Reconcile does. An unpinnable request is not needPin:
// there is nowhere for restore to put an answer, and it is already a scan
// finding rather than a missing pin.
func requestsNeedingPin(l *lock.Lock, scanned *lock.ScanResult) []lock.Request {
	var out []lock.Request
	for _, request := range scanned.Requests {
		if !request.Kind.Pinnable() {
			continue
		}
		if _, ok := l.Pins[request.Key]; !ok {
			out = append(out, request)
		}
	}
	return out
}

type requestReport struct {
	Request  string   `json:"request"`
	Kind     string   `json:"kind"`
	Scripts  []string `json:"scripts"`
	Artifact string   `json:"artifact,omitempty"`
	// Pinnable is false for a declaration nothing can pin, which is a finding.
	Pinnable bool `json:"pinnable"`
}

func validateReport(root string, l *lock.Lock, problems []lock.Problem, needPin []lock.Request, scanned *lock.ScanResult) any {
	report := struct {
		Root       string          `json:"root"`
		Valid      bool            `json:"valid"`
		Problems   []string        `json:"problems,omitempty"`
		Findings   []lock.Finding  `json:"findings,omitempty"`
		Unselected []requestReport `json:"needPin,omitempty"`
		Scripts    []string        `json:"scripts"`
	}{Root: root, Scripts: scanned.Scripts, Findings: scanned.Findings}
	for _, problem := range problems {
		report.Problems = append(report.Problems, problem.String())
	}
	if reason := missingBaseProblem(l); reason != "" {
		report.Problems = append(report.Problems, reason)
	}
	for _, request := range needPin {
		report.Unselected = append(report.Unselected, requestReport{
			Request: request.Key, Kind: string(request.Kind), Scripts: request.Scripts,
			Pinnable: true})
	}
	report.Valid = len(report.Problems) == 0 && len(report.Findings) == 0 && len(report.Unselected) == 0
	return report
}

// reportLockState prints what a reconcile left behind and fails while anything
// is needPin, so a partial lock is published but never reported as complete.
func reportLockState(result *orchestrate.LockResult, jsonOutput bool) error {
	root, l, scanned, pinned, failed := result.Root, result.Current, result.Scanned, result.Pinned, result.Failed
	if jsonOutput {
		report := struct {
			Root                   string              `json:"root"`
			Base                   string              `json:"base,omitempty"`
			Requests               []requestReport     `json:"requests"`
			Failed                 []string            `json:"failed,omitempty"`
			Findings               []lock.Finding      `json:"findings,omitempty"`
			UnpinnedHelperOverlays []string            `json:"unpinned_helper_overlays,omitempty"`
			ManualPinUsage         map[string][]string `json:"manual_pin_usage,omitempty"`
		}{Root: root, Findings: scanned.Findings,
			UnpinnedHelperOverlays: result.UnpinnedHelperOverlays, ManualPinUsage: result.ManualPinUsage}
		if pin, ok := l.Pins[lock.BaseKey]; ok {
			report.Base = pin.Artifact
		}
		for _, request := range scanned.Requests {
			entry := requestReport{Request: request.Key, Kind: string(request.Kind),
				Scripts: request.Scripts, Pinnable: request.Kind.Pinnable()}
			if pin, ok := l.Pins[request.Key]; ok {
				entry.Artifact = pin.Artifact
			}
			report.Requests = append(report.Requests, entry)
		}
		report.Failed = errorStrings(failed)
		if err := printJSON(report); err != nil {
			return err
		}
		return pinFailure(failed)
	}

	if pin, ok := l.Pins[lock.BaseKey]; ok {
		utils.PrintMessage("  %s → %s", utils.StyleName("base"), pin.Artifact)
	}
	for _, request := range scanned.Requests {
		switch pin, ok := l.Pins[request.Key]; {
		case ok:
			utils.PrintMessage("  %s → %s", utils.StyleName(request.Key), pin.Artifact)
		case !request.Kind.Pinnable():
			utils.PrintMessage("  %s → %s", utils.StyleName(request.Key), utils.StyleWarning("cannot be pinned"))
		}
	}
	// Reported after the map, because an alternative only makes sense once the
	// reader can see which one was taken.
	alternatives := false
	for _, entry := range pinned {
		for _, other := range entry.Others {
			alternatives = true
			utils.PrintMessage("    also installed: %s %s",
				utils.StyleWarning(other.Identity.Digest()), utils.StylePath(other.Path))
		}
	}
	if alternatives {
		utils.PrintHint("Pin one of those instead with %s.",
			utils.StyleAction("condatainer project pin <request> <identity>"))
	}
	for _, finding := range scanned.Findings {
		utils.PrintWarning("%s:%d: %s", finding.Script, finding.Line, finding.Reason)
	}
	printUnpublishedNote(l, pinned)
	printUnpinnedHelperOverlays(result.UnpinnedHelperOverlays)
	printManualPinUsage(result.ManualPinUsage)
	if len(failed) == 0 {
		utils.PrintSuccess("Every declaration is pinned (%d).", len(l.Pins))
		return nil
	}
	for _, err := range failed {
		utils.PrintError("%v", err)
	}
	return pinFailure(failed)
}

// pinFailure is the exit status: what could be pinned already was and is
// published, so this says the project is incomplete, never that nothing ran.
func pinFailure(failed []error) error {
	if len(failed) == 0 {
		return nil
	}
	return fmt.Errorf("%d declaration(s) could not be pinned", len(failed))
}

func errorStrings(errs []error) []string {
	out := make([]string, 0, len(errs))
	for _, err := range errs {
		out = append(out, err.Error())
	}
	return out
}

func printJSON(value any) error {
	encoder := json.NewEncoder(os.Stdout)
	encoder.SetIndent("", "  ")
	return encoder.Encode(value)
}

func newProjectRestoreCmd() *cobra.Command {
	var (
		jsonOutput    bool
		dryRun        bool
		noPrebuilt    bool
		keepBuildDeps bool
		matchMode     string
		replace       bool
		only          string
	)
	cmd := &cobra.Command{
		Use:   "restore",
		Short: "Make every locked artifact available on this machine",
		Long: `Reuses, downloads, or rebuilds each locked artifact until the project can run.

- A named pin is installed under its name; a path pin at the path it declares.
- Downloads from a recorded registry whenever one is available.
- A rebuild that carries scheduler directives is submitted as a job.
- Build dependencies are discarded at the end unless --keep-build-deps.

A path pin is the only file a restore can overwrite, so one already holding
something the lock does not name is refused unless you pass --replace.`,
		Args:         cobra.NoArgs,
		SilenceUsage: true,
		RunE: func(cmd *cobra.Command, _ []string) error {
			match, err := parseMatch(matchMode)
			if err != nil {
				return err
			}
			root, err := projectRoot(!jsonOutput)
			if err != nil {
				return err
			}
			current, err := lock.Load(root)
			if err != nil {
				return err
			}
			opts := restore.Options{
				Match: match, SkipPrebuilt: noPrebuilt, KeepBuildDeps: keepBuildDeps,
				Replace: replace, Only: only, SubmitJobs: config.Global.SubmitJob,
			}

			if dryRun {
				return reportPlan(restore.Compute(root, current, opts), jsonOutput)
			}
			report, runErr := restore.Run(cmd.Context(), root, current, opts)
			if err := reportRestore(report, jsonOutput); err != nil {
				return err
			}
			if errors.Is(runErr, restore.ErrJobsSubmitted) {
				utils.PrintNote("%d scheduler job(s) submitted. exiting with code %d",
					len(report.Submitted), ExitCodeJobsSubmitted)
				os.Exit(ExitCodeJobsSubmitted)
			}
			return runErr
		},
	}
	cmd.Flags().BoolVar(&jsonOutput, "json", false, "Print JSON")
	cmd.Flags().BoolVar(&dryRun, "dry-run", false, "Report what would happen and change nothing")
	cmd.Flags().BoolVar(&noPrebuilt, "no-prebuilt", false,
		"Build missing artifacts locally where possible instead of downloading a prebuilt")
	cmd.Flags().BoolVar(&keepBuildDeps, "keep-build-deps", false,
		"Keep build dependencies instead of discarding them at the end")
	cmd.Flags().BoolVar(&replace, "replace", false,
		"Overwrite a project path holding something the lock does not name")
	cmd.Flags().StringVar(&matchMode, "match", string(restore.MatchEquivalent),
		"Key a restored artifact must agree with: equivalent or identity")
	// Set by the job a submitted rebuild runs, so it produces exactly what it
	// was sent for. Nothing a person types.
	cmd.Flags().StringVar(&only, "only", "", "Restore one recorded artifact and its dependencies")
	_ = cmd.Flags().MarkHidden("only")
	return cmd
}

// parseMatch rejects an unknown mode rather than falling back to the default:
// a misspelled --match identity would silently restore under the looser rule.
func parseMatch(mode string) (restore.Match, error) {
	switch restore.Match(mode) {
	case restore.MatchEquivalent:
		return restore.MatchEquivalent, nil
	case restore.MatchIdentity:
		return restore.MatchIdentity, nil
	}
	return "", fmt.Errorf("unknown --match %q: use %q or %q",
		mode, restore.MatchEquivalent, restore.MatchIdentity)
}

// reportPlan prints a --dry-run plan and fails if it could not run.
func reportPlan(plan *restore.Plan, jsonOutput bool) error {
	if jsonOutput {
		if err := printJSON(plan); err != nil {
			return err
		}
	} else {
		for _, step := range plan.Steps {
			utils.PrintMessage("  %-11s %s → %s", step.Action, utils.StyleName(step.Name), planDestination(step))
			if step.Replaces != "" {
				utils.PrintWarning("    replaces %s already there", short(step.Replaces))
			}
		}
		for _, problem := range plan.Problems {
			utils.PrintError("%s", problem)
		}
		if plan.Complete() {
			utils.PrintSuccess("%d artifact(s) to acquire, %d already available.",
				plan.Work(), len(plan.Steps)-plan.Work())
		}
	}
	if !plan.Complete() {
		return fmt.Errorf("restore would not complete: %d problem(s)", len(plan.Problems))
	}
	return nil
}

// reportRestore prints what a restore did. The run's own error is the exit
// status; this only explains it.
func reportRestore(report *restore.Report, jsonOutput bool) error {
	if jsonOutput {
		return printJSON(report)
	}
	for _, result := range report.Results {
		line := fmt.Sprintf("  %-8s %s → %s", result.Outcome, utils.StyleName(result.Name), result.Path)
		if result.Transient {
			line += " (build dep, discarded)"
		}
		if result.Found != "" {
			line += utils.StyleWarning(fmt.Sprintf(" (equivalent, not %s)", short(result.Identity)))
		}
		utils.PrintMessage("%s", line)
	}
	for _, problem := range report.Problems {
		utils.PrintError("%s", problem)
	}
	for _, failure := range report.Failures {
		utils.PrintError("%s: %s", utils.StyleName(failure.Name), failure.Reason)
		for _, diff := range failure.Diffs {
			utils.PrintMessage("      %s", diff)
		}
	}
	// Blocked artifacts follow the failures and name the one cause, so a broken
	// dependency reads as one problem rather than as one per artifact above it.
	for _, blocked := range report.Blocked {
		utils.PrintWarning("%s: not attempted, %s failed", utils.StyleName(blocked.Name), blocked.Cause)
	}
	for _, submitted := range report.Submitted {
		verb := "submitted as"
		if submitted.Queued {
			verb = "already queued as"
		}
		line := fmt.Sprintf("  %-8s %s %s job %s", "pending", utils.StyleName(submitted.Name), verb, submitted.JobID)
		if len(submitted.DependsOn) > 0 {
			line += fmt.Sprintf(" (after %s)", strings.Join(submitted.DependsOn, ", "))
		}
		utils.PrintMessage("%s", line)
	}
	if report.Complete() {
		utils.PrintSuccess("Project restored: %d artifact(s).", len(report.Results))
	}
	return nil
}

// planDestination says where a step's result will land, which differs by what
// asked for it: a project path, an images root, or nowhere at all.
func planDestination(step restore.Step) string {
	switch {
	case step.Destination != "":
		return step.Destination
	case !step.Direct:
		return "build dep, discarded after"
	default:
		return "store or flat"
	}
}

// short renders a digest for a one-line summary.
func short(digest string) string {
	if len(digest) > 19 {
		return digest[:19]
	}
	return digest
}
