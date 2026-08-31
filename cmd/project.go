package cmd

import (
	"context"
	"encoding/json"
	"errors"
	"fmt"
	"os"
	"strings"

	"github.com/Justype/condatainer/internal/config"
	"github.com/Justype/condatainer/internal/logging"
	"github.com/Justype/condatainer/internal/project/lock"
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
	projectCmd.AddCommand(newProjectLockCmd(), newProjectPinCmd(), newProjectValidateCmd(),
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
			current, err := lock.Load(root)
			if err != nil {
				return err
			}
			scanned, err := lock.Scan(root, lock.ScanOptions{})
			if err != nil {
				return err
			}
			unpinned := lock.Reconcile(root, current, scanned)
			if err := lock.Publish(root, current); err != nil {
				return err
			}
			// Pinning is what makes this a lock rather than a scan: a
			// declaration names what is needed, and the lock has to say which
			// exact build answers it.
			pinned, failed := lock.PinAll(root, current, unpinned, lock.PinOptions{})
			for _, p := range pinned {
				recordUpstream(cmd.Context(), root, current, p)
			}
			if len(pinned) > 0 {
				if err := lock.Publish(root, current); err != nil {
					return err
				}
			}
			return reportLockState(root, current, scanned, pinned, failed, jsonOutput)
		},
	}
	cmd.Flags().BoolVar(&jsonOutput, "json", false, "Print JSON")
	return cmd
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
			upstream := recordUpstream(cmd.Context(), root, current, pinned)
			if err := lock.Apply(root, current, pinned); err != nil {
				return err
			}

			utils.PrintSuccess("Pinned %s", utils.StyleName(pinned.Name))
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
	return cmd
}

// recordUpstream adds a fetch location for every artifact this pin
// vendored that its own recipe collection already publishes at the exact same
// identity, and reports which, for display.
//
// Free in both senses: the bytes are already there, so restore downloads instead
// of rebuilding, and nothing was uploaded to make it so.
//
// It never fails the pin. No network, no configured source, no declared
// endpoint and no match all record nothing — locking has to work offline, and an
// absent remote costs a rebuild rather than an error.
func recordUpstream(ctx context.Context, root string, l *lock.Lock, pinned *lock.Pinned) map[string]string {
	cat, err := config.OpenCatalog(ctx)
	if err != nil {
		logging.FromContext(ctx).Debug("no catalog, so no upstream locations were recorded", "err", err)
		return nil
	}
	found := publish.Upstream(ctx, root, pinned.Vendored, cat)
	shown := make(map[string]string, len(found))
	for artifact, remotes := range found {
		for _, remote := range remotes {
			if err := l.AddRemote(artifact, remote); err != nil {
				logging.FromContext(ctx).Debug("could not record an upstream location", "artifact", artifact, "err", err)
				continue
			}
			shown[artifact] = remote.Repository
		}
	}
	return shown
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
			if len(problems) > 0 {
				// Publication never fills gaps from the live catalog: what is
				// published has to be what the checkout already describes.
				return fmt.Errorf("the project is not valid, so nothing was published:\n  %s",
					strings.Join(problemStrings(problems), "\n  "))
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
	utils.PrintMessage("%d to upload", plan.Uploads())
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
			if jsonOutput {
				if err := printJSON(validateReport(root, problems, needPin, scanned)); err != nil {
					return err
				}
				// The exit status is the answer; the format only changes how it
				// is explained.
				if total > 0 {
					return fmt.Errorf("project is not valid: %d problem(s)", total)
				}
				return nil
			}
			for _, problem := range problems {
				utils.PrintError("%s", problem.String())
			}
			for _, finding := range scanned.Findings {
				utils.PrintError("%s:%d: %s", finding.Script, finding.Line, finding.Reason)
			}
			for _, request := range needPin {
				utils.PrintError("%s is declared but not pinned (%s)",
					utils.StyleName(request.Key), strings.Join(request.Scripts, ", "))
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
// there is nowhere for restore to put an answer, and one left undeclared is
// already a scan finding rather than a missing pin.
func requestsNeedingPin(l *lock.Lock, scanned *lock.ScanResult) []lock.Request {
	var out []lock.Request
	for _, request := range scanned.Requests {
		if !request.Kind.Pinnable() || request.Unpinned {
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
	Unpinned bool     `json:"unpinned,omitempty"`
	Reason   string   `json:"reason,omitempty"`
}

func validateReport(root string, problems []lock.Problem, needPin []lock.Request, scanned *lock.ScanResult) any {
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
	for _, request := range needPin {
		report.Unselected = append(report.Unselected, requestReport{
			Request: request.Key, Kind: string(request.Kind), Scripts: request.Scripts})
	}
	report.Valid = len(report.Problems) == 0 && len(report.Findings) == 0 && len(report.Unselected) == 0
	return report
}

// reportLockState prints what a reconcile left behind and fails while anything
// is needPin, so a partial lock is published but never reported as complete.
func reportLockState(root string, l *lock.Lock, scanned *lock.ScanResult,
	pinned []*lock.Pinned, failed []error, jsonOutput bool) error {
	if jsonOutput {
		report := struct {
			Root     string          `json:"root"`
			Requests []requestReport `json:"requests"`
			Failed   []string        `json:"failed,omitempty"`
			Findings []lock.Finding  `json:"findings,omitempty"`
		}{Root: root, Findings: scanned.Findings}
		for _, request := range scanned.Requests {
			entry := requestReport{Request: request.Key, Kind: string(request.Kind),
				Scripts: request.Scripts, Unpinned: request.Unpinned, Reason: request.Reason}
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

	for _, request := range scanned.Requests {
		switch pin, ok := l.Pins[request.Key]; {
		case ok:
			utils.PrintMessage("  %s → %s", utils.StyleName(request.Key), pin.Artifact)
		case request.Unpinned:
			utils.PrintMessage("  %s → %s", utils.StyleName(request.Key), utils.StyleWarning("unpinned"))
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
	cmd.Flags().BoolVar(&noPrebuilt, "no-prebuilt", false, "Build missing artifacts locally instead of downloading a prebuilt")
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
			utils.PrintMessage("  %-7s %s → %s", step.Action, utils.StyleName(step.Name), planDestination(step))
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
