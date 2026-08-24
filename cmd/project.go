package cmd

import (
	"encoding/json"
	"errors"
	"fmt"
	"os"
	"strings"

	"github.com/Justype/condatainer/internal/config"
	"github.com/Justype/condatainer/internal/project/lock"
	"github.com/Justype/condatainer/internal/project/restore"
	"github.com/Justype/condatainer/internal/utils"
	"github.com/spf13/cobra"
)

// projectDir is the --project flag, shared by every subcommand.
var projectDir string

var projectCmd = &cobra.Command{
	Use:   "project",
	Short: "Pin a project's dependencies to exact artifact identities",
	Long: `Record which exact artifact satisfies each #DEP: a project declares.

cnt-lock/ belongs in Git: it holds the selection map plus the manifests and
rebuild sources each selected artifact was built from. It holds no payload and
nothing machine-local, so a checkout is a complete rebuild specification even on
a machine that has never run CondaTainer.`,
}

func init() {
	rootCmd.AddCommand(projectCmd)
	projectCmd.PersistentFlags().StringVar(&projectDir, "project", "", "project root (default: the current directory)")
	projectCmd.AddCommand(newProjectLockCmd(), newProjectValidateCmd(), newProjectRestoreCmd())
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
			err, "cnt-lock/", "mkdir cnt-lock", "--project DIR")
	}
	if err != nil {
		return "", err
	}
	if announce {
		utils.PrintMessage("Project: %s", utils.StylePath(root))
	}
	return root, nil
}

func newProjectLockCmd() *cobra.Command {
	var jsonOutput bool
	cmd := &cobra.Command{
		Use:   "lock",
		Short: "Rescan declarations and reconcile the lock",
		Long: `Rescan every script for #DEP: declarations and reconcile them with the lock.

Selections nothing declares any more are dropped, and so are selections whose
artifact no longer verifies. Nothing is selected automatically: choosing an
artifact is an explicit act, so unselected requests are reported for
'condatainer project lock select' to resolve.`,
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
			unselected := lock.Reconcile(root, current, scanned)
			if err := lock.Publish(root, current); err != nil {
				return err
			}
			return reportLockState(root, current, scanned, unselected, jsonOutput)
		},
	}
	cmd.Flags().BoolVar(&jsonOutput, "json", false, "print JSON")
	cmd.AddCommand(newProjectSelectCmd())
	return cmd
}

func newProjectSelectCmd() *cobra.Command {
	cmd := &cobra.Command{
		Use:   "select <request> <identity|path>",
		Short: "Select one exact artifact and vendor its source closure",
		Long: `Resolve a declaration to one exact local artifact and vendor its sources.

The target is an identity — scheme@sha256:<hex>, a full digest, or an
unambiguous prefix — or a path to an immutable .sqf, which may live outside the
project. What gets recorded is the identity; the path is only how the artifact
was found.

Every copy in every readable image root is a candidate, not just the nearest
one, so a selection can name a copy that ordinary name resolution would hide.`,
		Example: `  condatainer project lock select star/2.7.11b a31f902c12ab
  condatainer project lock select star/2.7.11b /shared/overlays/star.sqf`,
		Args:         cobra.ExactArgs(2),
		SilenceUsage: true,
		RunE: func(_ *cobra.Command, args []string) error {
			root, err := projectRoot(true)
			if err != nil {
				return err
			}
			current, err := lock.Load(root)
			if err != nil {
				return err
			}
			selected, err := lock.Select(root, args[0], args[1], lock.SelectOptions{})
			if err != nil {
				return err
			}
			if err := lock.Apply(root, current, selected); err != nil {
				return err
			}

			utils.PrintSuccess("Selected %s", utils.StyleName(selected.Name))
			utils.PrintMessage("  identity %s", selected.Identity.Digest())
			utils.PrintMessage("  read from %s", utils.StylePath(selected.Path))
			for _, artifact := range selected.Vendored {
				utils.PrintMessage("  vendored %s", artifact)
			}
			return nil
		},
	}
	return cmd
}

func newProjectValidateCmd() *cobra.Command {
	var jsonOutput bool
	cmd := &cobra.Command{
		Use:   "validate",
		Short: "Validate the lock and its vendored sources, using only the checkout",
		Long: `Check that a checkout is a complete, internally consistent rebuild specification.

Nothing is read outside the checkout — no installed overlay, no store, no
catalog, no configuration, no network — and no payload is needed, so this
succeeds or fails identically in CI on a machine with no images at all.`,
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
			unselected := unselectedRequests(current, scanned)

			total := len(problems) + len(scanned.Findings) + len(unselected)
			if jsonOutput {
				if err := printJSON(validateReport(root, problems, unselected, scanned)); err != nil {
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
			for _, request := range unselected {
				utils.PrintError("%s is declared but not selected (%s)",
					utils.StyleName(request.Key), strings.Join(request.Scripts, ", "))
			}
			if total > 0 {
				utils.PrintHint("Run %s to choose an artifact for each unselected request.",
					utils.StyleAction("condatainer project lock select <request> <identity>"))
				return fmt.Errorf("project is not valid: %d problem(s)", total)
			}
			utils.PrintSuccess("Project is valid: %d selection(s), %d script(s)",
				len(current.Selections), len(scanned.Scripts))
			return nil
		},
	}
	cmd.Flags().BoolVar(&jsonOutput, "json", false, "print JSON")
	return cmd
}

// unselectedRequests reports declarations with no selection, without mutating
// the lock the way Reconcile does. An unpinnable request is not unselected:
// there is nowhere for restore to put an answer, and one left undeclared is
// already a scan finding rather than a missing selection.
func unselectedRequests(l *lock.Lock, scanned *lock.ScanResult) []lock.Request {
	var out []lock.Request
	for _, request := range scanned.Requests {
		if !request.Kind.Pinnable() || request.Unpinned {
			continue
		}
		if _, ok := l.Selections[request.Key]; !ok {
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

func validateReport(root string, problems []lock.Problem, unselected []lock.Request, scanned *lock.ScanResult) any {
	report := struct {
		Root       string          `json:"root"`
		Valid      bool            `json:"valid"`
		Problems   []string        `json:"problems,omitempty"`
		Findings   []lock.Finding  `json:"findings,omitempty"`
		Unselected []requestReport `json:"unselected,omitempty"`
		Scripts    []string        `json:"scripts"`
	}{Root: root, Scripts: scanned.Scripts, Findings: scanned.Findings}
	for _, problem := range problems {
		report.Problems = append(report.Problems, problem.String())
	}
	for _, request := range unselected {
		report.Unselected = append(report.Unselected, requestReport{
			Request: request.Key, Kind: string(request.Kind), Scripts: request.Scripts})
	}
	report.Valid = len(report.Problems) == 0 && len(report.Findings) == 0 && len(report.Unselected) == 0
	return report
}

// reportLockState prints what a reconcile left behind and fails while anything
// is unselected, so a partial lock is published but never reported as complete.
func reportLockState(root string, l *lock.Lock, scanned *lock.ScanResult, unselected []lock.Request, jsonOutput bool) error {
	if jsonOutput {
		report := struct {
			Root       string          `json:"root"`
			Selections []requestReport `json:"selections"`
			Unselected []requestReport `json:"unselected,omitempty"`
			Findings   []lock.Finding  `json:"findings,omitempty"`
		}{Root: root, Findings: scanned.Findings}
		for _, request := range scanned.Requests {
			entry := requestReport{Request: request.Key, Kind: string(request.Kind),
				Scripts: request.Scripts, Unpinned: request.Unpinned, Reason: request.Reason}
			if selection, ok := l.Selections[request.Key]; ok {
				entry.Artifact = selection.Artifact
			}
			report.Selections = append(report.Selections, entry)
		}
		for _, request := range unselected {
			report.Unselected = append(report.Unselected, requestReport{
				Request: request.Key, Kind: string(request.Kind), Scripts: request.Scripts})
		}
		if err := printJSON(report); err != nil {
			return err
		}
		if len(unselected) > 0 {
			return fmt.Errorf("%d request(s) are not selected", len(unselected))
		}
		return nil
	}

	for _, request := range scanned.Requests {
		switch selection, ok := l.Selections[request.Key]; {
		case ok:
			utils.PrintMessage("  %s → %s", utils.StyleName(request.Key), selection.Artifact)
		case request.Unpinned:
			utils.PrintMessage("  %s → %s", utils.StyleName(request.Key), utils.StyleWarning("unpinned"))
		}
	}
	for _, finding := range scanned.Findings {
		utils.PrintWarning("%s:%d: %s", finding.Script, finding.Line, finding.Reason)
	}
	if len(unselected) == 0 {
		utils.PrintSuccess("Every declaration is selected (%d).", len(l.Selections))
		return nil
	}
	utils.PrintMessage("Unselected:")
	for _, request := range unselected {
		utils.PrintMessage("  - %s (%s)", utils.StyleWarning(request.Key), strings.Join(request.Scripts, ", "))
	}
	utils.PrintHint("Choose one with %s.", utils.StyleAction("condatainer project lock select <request> <identity>"))
	return fmt.Errorf("%d request(s) are not selected", len(unselected))
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
		only          string
	)
	cmd := &cobra.Command{
		Use:   "restore",
		Short: "Make every locked artifact available on this machine",
		Long: `Reuse, fetch, or rebuild each locked artifact until the project can run.

Nothing is mounted and cnt-lock/ is never modified. A selected artifact lands
where the store's destination rule puts it — the flat name when it is free,
store/ when that name is already held at a different identity — and a
project-path selection is materialized at exactly the path it declares.

An artifact nothing selected is a build dependency: it exists only so its dependent
can be built, so it is produced in a temporary directory and removed when the
restore ends. One already installed is adopted in place and never copied.
--keep-build-deps installs newly produced ones through the ordinary rule instead.

--no-prebuilt builds every missing artifact from source rather than downloading
one the lock records. It is not an offline mode: every build needs the network —
a Conda replay downloads each pinned package, a recipe fetches its own sources —
which is what the proxy is for on a compute node without egress. An artifact
already present is adopted either way.

--match decides which key a restored artifact must agree with. The default,
equivalent, accepts anything that can substitute for what was locked, which is
what the equivalence key exists to decide; identity accepts only the exact build
the lock names, which additionally requires a Conda artifact to replay its
explicit.txt to the byte and a data artifact to be rebuilt in the same
environment.

A rebuild whose recipe carries scheduler directives is submitted rather than run
here, and so is anything waiting on it; dependencies become afterok edges. The
command then exits with the jobs-submitted code, having made nothing available
yet. Re-running the restore is how it resumes: a job still queued is reported
rather than submitted twice, and one that finished is adopted.

Restore is atomic per artifact, not across the project: a later failure leaves
earlier results in place, and re-running adopts them.`,
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
				Only: only, SubmitJobs: config.Global.SubmitJob,
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
	cmd.Flags().BoolVar(&jsonOutput, "json", false, "print JSON")
	cmd.Flags().BoolVar(&dryRun, "dry-run", false, "report what would happen and acquire nothing")
	cmd.Flags().BoolVar(&noPrebuilt, "no-prebuilt", false, "build missing artifacts from source instead of downloading a prebuilt")
	cmd.Flags().BoolVar(&keepBuildDeps, "keep-build-deps", false,
		"install build dependencies instead of discarding them when the restore ends")
	cmd.Flags().StringVar(&matchMode, "match", string(restore.MatchEquivalent),
		"key a restored artifact must agree with: equivalent or identity")
	// Set by the job a submitted rebuild runs, so it produces exactly what it
	// was sent for. Nothing a person types.
	cmd.Flags().StringVar(&only, "only", "", "restore one vendored artifact and its closure")
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
