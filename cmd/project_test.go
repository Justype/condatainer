package cmd

import (
	"encoding/json"
	"errors"
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/Justype/condatainer/catalog"
	"github.com/Justype/condatainer/internal/artifact/capsule"
	"github.com/Justype/condatainer/internal/artifact/key"
	"github.com/Justype/condatainer/internal/artifact/meta"
	"github.com/Justype/condatainer/internal/project/lock"
	"github.com/Justype/condatainer/internal/project/restore"
	"github.com/spf13/cobra"
	"github.com/spf13/pflag"
)

// newProject builds a project root with cnt-lock/ and chdirs into it, so the
// commands exercise their own upward root discovery.
func newProject(t *testing.T) string {
	t.Helper()
	root := t.TempDir()
	if err := os.MkdirAll(filepath.Join(root, lock.DirName), 0o775); err != nil {
		t.Fatal(err)
	}
	t.Chdir(root)
	projectDir = ""
	t.Cleanup(func() { projectDir = "" })
	return root
}

func writeScript(t *testing.T, root, rel, body string) {
	t.Helper()
	path := filepath.Join(root, rel)
	if err := os.MkdirAll(filepath.Dir(path), 0o775); err != nil {
		t.Fatal(err)
	}
	if err := os.WriteFile(path, []byte(body), 0o664); err != nil {
		t.Fatal(err)
	}
}

// vendorArtifact puts a verifiable artifact into cnt-lock/provenance/ and returns
// its relative path, standing in for a selection that already happened.
func vendorArtifact(t *testing.T, root, name, recipe string) string {
	t.Helper()
	manifest := meta.Manifest{
		SchemaVersion: meta.SchemaVersion,
		Name:          name,
		Type:          catalog.TypeApp,
		BuildType:     "script",
		// This host, in uname form: a fixture claiming another architecture
		// would be refused by restore's platform check before anything else.
		Platform: meta.NativePlatform(),
		Source:   meta.Source{Files: []string{meta.RecipeFileName}},
	}
	derived, err := key.Generate(manifest, key.Sources{meta.RecipeFileName: []byte(recipe)})
	if err != nil {
		t.Fatal(err)
	}
	manifest.Keys = derived.Keys()
	data, err := json.MarshalIndent(manifest, "", "  ")
	if err != nil {
		t.Fatal(err)
	}
	relative, err := lock.StageEntry(root, capsule.EntryName(name, manifest.Keys.Identity.Digest()),
		map[string][]byte{meta.FileName: data, meta.RecipeFileName: []byte(recipe)})
	if err != nil {
		t.Fatal(err)
	}
	return relative
}

// resetFlags returns every project flag to its default. Cobra builds the
// command tree once, so a flag set by one test would otherwise persist into the
// next and make it pass for the wrong reason.
func resetFlags(cmd *cobra.Command) {
	cmd.Flags().VisitAll(func(f *pflag.Flag) {
		_ = f.Value.Set(f.DefValue)
		f.Changed = false
	})
	for _, child := range cmd.Commands() {
		resetFlags(child)
	}
}

// run executes the project command tree with args and captures stdout.
func run(t *testing.T, args ...string) (string, error) {
	t.Helper()
	resetFlags(projectCmd)
	stdout := os.Stdout
	read, write, err := os.Pipe()
	if err != nil {
		t.Fatal(err)
	}
	os.Stdout = write

	rootCmd.SetArgs(args)
	runErr := rootCmd.Execute()

	write.Close()
	os.Stdout = stdout
	var out strings.Builder
	buf := make([]byte, 4096)
	for {
		n, readErr := read.Read(buf)
		out.Write(buf[:n])
		if readErr != nil {
			break
		}
	}
	read.Close()
	return out.String(), runErr
}

// A declared but needPin request leaves a valid partial lock and a nonzero
// exit: the lock is published, but never reported as complete.
func TestProjectLockPublishesPartiallyAndFails(t *testing.T) {
	root := newProject(t)
	writeScript(t, root, "run.sh", "#DEP: star/2.7.11b\nrun\n")

	if _, err := run(t, "project", "lock", "--json"); err == nil {
		t.Fatal("project lock succeeded with an needPin request")
	}
	if _, err := os.Stat(lock.FilePath(root)); err != nil {
		t.Fatalf("no lock was published: %v", err)
	}
	loaded, err := lock.Load(root)
	if err != nil {
		t.Fatalf("the published lock does not load: %v", err)
	}
	if len(loaded.Pins) != 0 {
		t.Fatalf("selections = %#v, want none invented", loaded.Pins)
	}
}

func TestProjectValidateSucceedsOnACompleteProject(t *testing.T) {
	root := newProject(t)
	writeScript(t, root, "run.sh", "#DEP: star/2.7.11b\nrun\n")
	artifact := vendorArtifact(t, root, "star/2.7.11b", "echo star\n")
	base := vendorArtifact(t, root, "ubuntu24/base", "echo base\n")

	l := lock.New()
	l.Pins["star/2.7.11b"] = lock.PinEntry{Artifact: artifact}
	l.Pins[lock.BaseKey] = lock.PinEntry{Artifact: base}
	if err := lock.Publish(root, l); err != nil {
		t.Fatal(err)
	}

	if _, err := run(t, "project", "validate"); err != nil {
		t.Fatalf("validate failed on a complete project: %v", err)
	}
}

// --installed asks what plain validate does not: whether the pins are on this
// machine. A complete lock passes the first and fails the second when nothing
// is installed.
func TestProjectValidateInstalledFailsWhenAPinIsNotInstalled(t *testing.T) {
	root := newProject(t)
	writeScript(t, root, "run.sh", "#DEP: star/2.7.11b\nrun\n")
	artifact := vendorArtifact(t, root, "star/2.7.11b", "echo star\n")
	base := vendorArtifact(t, root, "ubuntu24/base", "echo base\n")

	l := lock.New()
	l.Pins["star/2.7.11b"] = lock.PinEntry{Artifact: artifact}
	l.Pins[lock.BaseKey] = lock.PinEntry{Artifact: base}
	if err := lock.Publish(root, l); err != nil {
		t.Fatal(err)
	}

	if _, err := run(t, "project", "validate"); err != nil {
		t.Fatalf("plain validate failed: %v", err)
	}
	out, err := run(t, "project", "validate", "--installed", "--json")
	if err == nil {
		t.Fatalf("--installed passed with nothing installed:\n%s", out)
	}
	if !strings.Contains(out, "is not installed here") {
		t.Errorf("--installed does not say which pin is missing:\n%s", out)
	}
}

// A declaration the lock cannot reproduce fails validate, even when every other
// declaration is pinned. Nothing in a script silences it: freezing the overlay
// is what closes it, and the finding says so.
func TestProjectValidateFailsOnAnUnpinnableDeclaration(t *testing.T) {
	root := newProject(t)
	writeScript(t, root, "run.sh", "#DEP: star/2.7.11b\n#DEP: env.img\nrun\n")
	artifact := vendorArtifact(t, root, "star/2.7.11b", "echo star\n")

	l := lock.New()
	l.Pins["star/2.7.11b"] = lock.PinEntry{Artifact: artifact}
	if err := lock.Publish(root, l); err != nil {
		t.Fatal(err)
	}

	out, err := runConsole(t, "project", "validate")
	if err == nil {
		t.Fatal("validate passed a project whose environment is a writable overlay")
	}
	if !strings.Contains(out, "overlay freeze") {
		t.Errorf("validate does not point at freeze:\n%s", out)
	}
}

// Push holds a project to what validate asks, so a lock that cannot reproduce
// the project is never uploaded. The refusal is reached before any registry is
// contacted.
func TestProjectPushRefusesWhatValidateRefuses(t *testing.T) {
	root := newProject(t)
	writeScript(t, root, "run.sh", "#DEP: star/2.7.11b\n#DEP: env.img\nrun\n")
	artifact := vendorArtifact(t, root, "star/2.7.11b", "echo star\n")

	l := lock.New()
	l.Pins["star/2.7.11b"] = lock.PinEntry{Artifact: artifact}
	if err := lock.Publish(root, l); err != nil {
		t.Fatal(err)
	}

	_, err := runConsole(t, "project", "push")
	if err == nil {
		t.Fatal("push published a project validate rejects")
	}
	if !strings.Contains(err.Error(), "env.img") {
		t.Errorf("the refusal does not name the declaration: %v", err)
	}
}

// A #DEP: below the header is inert, and validate says so rather than passing.
func TestProjectValidateFailsOnAnInertDeclaration(t *testing.T) {
	root := newProject(t)
	writeScript(t, root, "run.sh", "#DEP: star/2.7.11b\nrun\n#DEP: cutadapt/5.0\n")
	artifact := vendorArtifact(t, root, "star/2.7.11b", "echo star\n")

	l := lock.New()
	l.Pins["star/2.7.11b"] = lock.PinEntry{Artifact: artifact}
	if err := lock.Publish(root, l); err != nil {
		t.Fatal(err)
	}

	if _, err := run(t, "project", "validate"); err == nil {
		t.Fatal("validate passed a project with an inert declaration")
	}
}

func TestProjectValidateJSONReportsWhatIsWrong(t *testing.T) {
	root := newProject(t)
	writeScript(t, root, "run.sh", "#DEP: star/2.7.11b\nrun\n")

	out, err := run(t, "project", "validate", "--json")
	if err == nil {
		t.Fatal("validate passed with nothing selected")
	}
	var report struct {
		Valid      bool `json:"valid"`
		Unselected []struct {
			Request string   `json:"request"`
			Scripts []string `json:"scripts"`
		} `json:"needPin"`
	}
	if err := json.Unmarshal([]byte(out), &report); err != nil {
		t.Fatalf("output is not JSON: %v\n%s", err, out)
	}
	if report.Valid {
		t.Error("report claims valid")
	}
	if len(report.Unselected) != 1 || report.Unselected[0].Request != "star/2.7.11b" {
		t.Fatalf("needPin = %#v", report.Unselected)
	}
	if len(report.Unselected[0].Scripts) != 1 || report.Unselected[0].Scripts[0] != "run.sh" {
		t.Errorf("scripts = %v, want the declaring script", report.Unselected[0].Scripts)
	}
}

// Reconcile drops a selection nothing declares any more. The lock is
// published right after Reconcile, before the base pin is even attempted, so
// this holds independent of whether a root can be derived.
func TestProjectLockDropsAStaleSelection(t *testing.T) {
	root := newProject(t)
	artifact := vendorArtifact(t, root, "star/2.7.11b", "echo star\n")
	l := lock.New()
	l.Pins["star/2.7.11b"] = lock.PinEntry{Artifact: artifact}
	if err := lock.Publish(root, l); err != nil {
		t.Fatal(err)
	}
	writeScript(t, root, "run.sh", "echo no declarations\n")

	// No default_distro is configured, so the run still fails deriving the
	// base pin — but that failure is reported alongside whatever Reconcile
	// already published, not instead of it.
	if _, err := run(t, "project", "lock"); err == nil {
		t.Fatal("project lock succeeded with no default_distro configured")
	}
	loaded, err := lock.Load(root)
	if err != nil {
		t.Fatal(err)
	}
	if len(loaded.Pins) != 0 {
		t.Fatalf("a selection nothing declares survived: %#v", loaded.Pins)
	}
	if _, err := os.Stat(filepath.Join(lock.Dir(root), artifact)); !os.IsNotExist(err) {
		t.Errorf("its artifact was not pruned: %v", err)
	}
}

func TestProjectFlagSelectsTheNamedRoot(t *testing.T) {
	root := newProject(t)
	writeScript(t, root, "run.sh", "#DEP: star/2.7.11b\nrun\n")
	artifact := vendorArtifact(t, root, "star/2.7.11b", "echo star\n")
	base := vendorArtifact(t, root, "ubuntu24/base", "echo base\n")
	l := lock.New()
	l.Pins["star/2.7.11b"] = lock.PinEntry{Artifact: artifact}
	l.Pins[lock.BaseKey] = lock.PinEntry{Artifact: base}
	if err := lock.Publish(root, l); err != nil {
		t.Fatal(err)
	}

	elsewhere := t.TempDir()
	t.Chdir(elsewhere)
	if _, err := run(t, "project", "validate", "--project", root); err != nil {
		t.Fatalf("--project did not select the named root: %v", err)
	}
}

// --dry-run reports the plan and acquires nothing. The project below has no
// images at all, so every step is work and none of it may happen.
func TestProjectRestoreDryRunAcquiresNothing(t *testing.T) {
	root := newProject(t)
	writeScript(t, root, "run.sh", "#DEP: star/2.7.11b\nrun\n")
	relative := vendorArtifact(t, root, "star/2.7.11b", "echo star\n")
	l := lock.New()
	l.Pins["star/2.7.11b"] = lock.PinEntry{Artifact: relative}
	if err := lock.Publish(root, l); err != nil {
		t.Fatal(err)
	}

	before := lockBytes(t, root)
	if out, err := run(t, "project", "restore", "--dry-run", "--no-prebuilt"); err != nil {
		t.Fatalf("dry run failed: %v\n%s", err, out)
	}
	if err := filepath.WalkDir(root, func(path string, entry os.DirEntry, err error) error {
		if err == nil && strings.HasSuffix(entry.Name(), ".sqf") {
			t.Errorf("a dry run produced %s", path)
		}
		return err
	}); err != nil {
		t.Fatal(err)
	}
	if after := lockBytes(t, root); after != before {
		t.Error("a dry run rewrote the lock")
	}
}

func lockBytes(t *testing.T, root string) string {
	t.Helper()
	data, err := os.ReadFile(lock.FilePath(root))
	if err != nil {
		t.Fatal(err)
	}
	return string(data)
}

// A misspelled mode is refused rather than silently restoring under the looser
// default, which is the whole point of asking for the stricter one.
func TestProjectRestoreRejectsAnUnknownMatchMode(t *testing.T) {
	root := newProject(t)
	writeScript(t, root, "run.sh", "#DEP: star/2.7.11b\nrun\n")

	out, err := run(t, "project", "restore", "--match", "exact")
	if err == nil {
		t.Fatalf("an unknown --match was accepted: %s", out)
	}
	if !strings.Contains(err.Error(), "identity") {
		t.Errorf("error does not name the valid modes: %v", err)
	}
}

func TestProjectRestoreDryRunJSONCarriesTheMatchMode(t *testing.T) {
	root := newProject(t)
	writeScript(t, root, "run.sh", "#DEP: star/2.7.11b\nrun\n")
	relative := vendorArtifact(t, root, "star/2.7.11b", "echo star\n")
	l := lock.New()
	l.Pins["star/2.7.11b"] = lock.PinEntry{Artifact: relative}
	if err := lock.Publish(root, l); err != nil {
		t.Fatal(err)
	}

	out, err := run(t, "project", "restore", "--dry-run", "--no-prebuilt", "--match", "identity", "--json")
	if err != nil {
		t.Fatalf("dry run failed: %v\n%s", err, out)
	}
	var report struct {
		Match string `json:"match"`
		Steps []struct {
			Name   string `json:"name"`
			Action string `json:"action"`
		} `json:"steps"`
	}
	if err := json.Unmarshal([]byte(out), &report); err != nil {
		t.Fatalf("output is not JSON: %v\n%s", err, out)
	}
	if report.Match != "identity" {
		t.Errorf("match = %q, want the requested mode", report.Match)
	}
	if len(report.Steps) != 1 || report.Steps[0].Action != "build" {
		t.Errorf("steps = %#v", report.Steps)
	}
}

// An unsound lock stops restore before it acquires anything, and says why.
func TestProjectRestoreFailsOnAnInvalidLock(t *testing.T) {
	root := newProject(t)
	writeScript(t, root, "run.sh", "#DEP: star/2.7.11b\nrun\n")
	l := lock.New()
	l.Pins["star/2.7.11b"] = lock.PinEntry{Artifact: "provenance/star--2.7.11b@000000000000"}
	if err := lock.Publish(root, l); err != nil {
		t.Fatal(err)
	}

	if out, err := run(t, "project", "restore", "--no-prebuilt"); err == nil {
		t.Fatalf("restore accepted an invalid lock: %s", out)
	}
}

// Scanning a directory for declarations is what makes it a project, so lock has
// nothing to refuse. Every other subcommand acts on selections that must already
// exist, and must not seed a project from a mistyped path.
func TestOnlyLockCreatesAProject(t *testing.T) {
	prev := projectDir
	t.Cleanup(func() { projectDir = prev })
	projectDir = ""

	chdir := func(t *testing.T) string {
		t.Helper()
		dir := t.TempDir()
		wd, err := os.Getwd()
		if err != nil {
			t.Fatal(err)
		}
		if err := os.Chdir(dir); err != nil {
			t.Fatal(err)
		}
		t.Cleanup(func() { _ = os.Chdir(wd) })
		// macOS hands out a symlinked temp path; Getwd resolves it.
		resolved, err := os.Getwd()
		if err != nil {
			t.Fatal(err)
		}
		return resolved
	}

	dir := chdir(t)
	root, err := projectRootOrInit(false)
	if err != nil || root != dir {
		t.Fatalf("projectRootOrInit = (%q, %v), want %q", root, err, dir)
	}

	if _, err := projectRoot(false); !errors.Is(err, lock.ErrNoProject) {
		t.Fatalf("projectRoot = %v, want ErrNoProject", err)
	}
	if _, err := os.Stat(filepath.Join(dir, lock.DirName)); !os.IsNotExist(err) {
		t.Error("resolving a root created cnt-lock/; only publishing should")
	}
}

// runConsole captures what utils.Print* writes, which is stderr — run captures
// stdout, where only --json output lands.
func runConsole(t *testing.T, args ...string) (string, error) {
	t.Helper()
	read, write, err := os.Pipe()
	if err != nil {
		t.Fatal(err)
	}
	stderr := os.Stderr
	os.Stderr = write

	_, runErr := run(t, args...)

	write.Close()
	os.Stderr = stderr
	var out strings.Builder
	buf := make([]byte, 4096)
	for {
		n, readErr := read.Read(buf)
		out.Write(buf[:n])
		if readErr != nil {
			break
		}
	}
	read.Close()
	return out.String(), runErr
}

// A manual pin is what `project lock` cannot sweep, so `project unpin` is the
// only way it goes — and what it alone vendored goes with it.
func TestProjectUnpinRemovesAManualPinAndItsProvenance(t *testing.T) {
	root := newProject(t)
	artifact := vendorArtifact(t, root, "ubuntu24/build-essential", "echo apt\n")
	l := lock.New()
	l.Pins["ubuntu24/build-essential"] = lock.PinEntry{Artifact: artifact, Manual: true}
	if err := lock.Publish(root, l); err != nil {
		t.Fatal(err)
	}

	out, err := runConsole(t, "project", "unpin", "ubuntu24/build-essential")
	if err != nil {
		t.Fatalf("unpin failed: %v", err)
	}
	if !strings.Contains(out, artifact) {
		t.Errorf("the pruned artifact was not reported:\n%s", out)
	}
	loaded, err := lock.Load(root)
	if err != nil {
		t.Fatal(err)
	}
	if len(loaded.Pins) != 0 {
		t.Fatalf("the pin survived: %#v", loaded.Pins)
	}
	if _, err := os.Stat(filepath.Join(lock.Dir(root), artifact)); !os.IsNotExist(err) {
		t.Errorf("its provenance was not pruned: %v", err)
	}
}

// A declared pin is removed by removing the declaration: unpinning one here
// would be undone by the next lock, possibly at a different identity.
func TestProjectUnpinRefusesADeclaredPin(t *testing.T) {
	root := newProject(t)
	writeScript(t, root, "run.sh", "#DEP: star/2.7.11b\nrun\n")
	artifact := vendorArtifact(t, root, "star/2.7.11b", "echo star\n")
	l := lock.New()
	l.Pins["star/2.7.11b"] = lock.PinEntry{Artifact: artifact}
	if err := lock.Publish(root, l); err != nil {
		t.Fatal(err)
	}

	_, err := run(t, "project", "unpin", "star/2.7.11b")
	if err == nil {
		t.Fatal("a declared pin was unpinned")
	}
	if !strings.Contains(err.Error(), "project lock") {
		t.Errorf("the refusal does not name the way to remove it: %v", err)
	}
	loaded, err := lock.Load(root)
	if err != nil {
		t.Fatal(err)
	}
	if _, ok := loaded.Pins["star/2.7.11b"]; !ok {
		t.Error("the refusal still removed the pin")
	}
}

func TestProjectUnpinRefusesWhatIsNotPinned(t *testing.T) {
	newProject(t)
	if _, err := run(t, "project", "unpin", "star/2.7.11b"); err == nil {
		t.Fatal("unpinning nothing succeeded")
	}
}

// The listing reads the lock and the records beside it, and is the only view
// that shows a manual pin — no scan produces one.
func TestProjectListShowsEveryPinAndMarksTheManualOne(t *testing.T) {
	root := newProject(t)
	declared := vendorArtifact(t, root, "star/2.7.11b", "echo star\n")
	manual := vendorArtifact(t, root, "ubuntu24/build-essential", "echo apt\n")
	l := lock.New()
	l.Pins["star/2.7.11b"] = lock.PinEntry{Artifact: declared}
	l.Pins["ubuntu24/build-essential"] = lock.PinEntry{Artifact: manual, Manual: true}
	if err := lock.Publish(root, l); err != nil {
		t.Fatal(err)
	}

	rows, err := run(t, "project", "list")
	if err != nil {
		t.Fatalf("list failed: %v", err)
	}
	for _, want := range []string{"star/2.7.11b", "ubuntu24/build-essential", "manual"} {
		if !strings.Contains(rows, want) {
			t.Errorf("stdout does not carry %q:\n%s", want, rows)
		}
	}
	// A named pin's key is its name, so the listing must not print it twice.
	if strings.Contains(rows, "star/2.7.11b (star/2.7.11b)") {
		t.Errorf("the name was repeated after the key it equals:\n%s", rows)
	}
	console, err := runConsole(t, "project", "list")
	if err != nil {
		t.Fatalf("list failed: %v", err)
	}
	if !strings.Contains(console, "2 pin(s), 1 manual") {
		t.Errorf("the summary is not on the console:\n%s", console)
	}
}

// JSON carries the flag as data, never as the styled word the text view prints.
func TestProjectListJSONReportsTheManualFlag(t *testing.T) {
	root := newProject(t)
	artifact := vendorArtifact(t, root, "ubuntu24/build-essential", "echo apt\n")
	l := lock.New()
	l.Pins["ubuntu24/build-essential"] = lock.PinEntry{Artifact: artifact, Manual: true}
	if err := lock.Publish(root, l); err != nil {
		t.Fatal(err)
	}

	out, err := run(t, "project", "list", "--json")
	if err != nil {
		t.Fatalf("list --json failed: %v", err)
	}
	var report struct {
		Pins []struct {
			Request  string `json:"request"`
			Name     string `json:"name"`
			Identity string `json:"identity"`
			Manual   bool   `json:"manual"`
		} `json:"pins"`
	}
	if err := json.Unmarshal([]byte(out), &report); err != nil {
		t.Fatalf("cannot decode: %v\n%s", err, out)
	}
	if len(report.Pins) != 1 {
		t.Fatalf("pins = %#v", report.Pins)
	}
	pin := report.Pins[0]
	if !pin.Manual || pin.Name != "ubuntu24/build-essential" || !strings.HasPrefix(pin.Identity, "sha256:") {
		t.Errorf("pin = %#v", pin)
	}
}

// A pin whose artifact is gone is still listed: the key is what addresses it,
// and hiding it would hide the thing to unpin.
func TestProjectListReportsAnUnreadableArtifact(t *testing.T) {
	root := newProject(t)
	artifact := vendorArtifact(t, root, "star/2.7.11b", "echo star\n")
	l := lock.New()
	l.Pins["star/2.7.11b"] = lock.PinEntry{Artifact: artifact, Manual: true}
	if err := lock.Publish(root, l); err != nil {
		t.Fatal(err)
	}
	if err := os.RemoveAll(filepath.Join(lock.Dir(root), artifact)); err != nil {
		t.Fatal(err)
	}

	rows, err := run(t, "project", "list")
	if err != nil {
		t.Fatalf("list failed: %v", err)
	}
	if !strings.Contains(rows, "star/2.7.11b") || !strings.Contains(rows, "unreadable") {
		t.Errorf("a broken pin was not listed as one:\n%s", rows)
	}
}

// A path key addresses a file and says nothing about what is in it, so the name
// is printed after it — the one case where the second column carries something.
func TestProjectListNamesWhatAPathPinHolds(t *testing.T) {
	root := newProject(t)
	artifact := vendorArtifact(t, root, "combined/1.0", "echo combined\n")
	l := lock.New()
	l.Pins["path:overlays/combined.sqf"] = lock.PinEntry{Artifact: artifact, Manual: true}
	if err := lock.Publish(root, l); err != nil {
		t.Fatal(err)
	}

	rows, err := run(t, "project", "list")
	if err != nil {
		t.Fatalf("list failed: %v", err)
	}
	if !strings.Contains(rows, "path:overlays/combined.sqf (combined/1.0)") {
		t.Errorf("a path pin does not name its artifact:\n%s", rows)
	}
}

// An equivalent artifact stands in, and the note says which inputs differ. The
// field names hold colons themselves, so only the colon-and-space ends one.
func TestEquivalentNoteNamesTheInputsThatDiffer(t *testing.T) {
	result := restore.Result{
		Identity: "sha256:41ab1c2d3e4f5a6b7c8d",
		Found:    "sha256:9f2b7e04a1c3d5e6f708",
		Diffs: []string{
			"src:gtf: 9f2b7e04a1c3 -> 3c81d5e2b7a4",
			"dep:samtools/1.23.1: 11aa22bb33cc -> 44dd55ee66ff",
		},
	}
	got := equivalentNote(result)
	for _, want := range []string{"equivalent, not sha256:41ab1c2d3e4f;", "differs: src:gtf, dep:samtools/1.23.1"} {
		if !strings.Contains(got, want) {
			t.Errorf("note = %q, want it to contain %q", got, want)
		}
	}
	// Digests are for --json, not a restore of many artifacts.
	if strings.Contains(got, "3c81d5e2b7a4") {
		t.Errorf("note = %q carries a digest", got)
	}
}

// With no diffs to name, the note is what it was before diffs existed.
func TestEquivalentNoteWithoutDiffsNamesNothing(t *testing.T) {
	got := equivalentNote(restore.Result{Identity: "sha256:41ab1c2d3e4f5a6b7c8d"})
	if want := " (equivalent, not sha256:41ab1c2d3e4f)"; got != want {
		t.Errorf("note = %q, want %q", got, want)
	}
}

// --verify reads a .sqf; anything else is refused before it reads a byte.
func TestInfoVerifyRefusesAnythingButASquashFS(t *testing.T) {
	img := filepath.Join(t.TempDir(), "env.img")
	if err := os.WriteFile(img, nil, 0o644); err != nil {
		t.Fatal(err)
	}
	out, err := runConsole(t, "info", img, "--verify")
	if err == nil {
		t.Fatalf("--verify accepted an .img:\n%s", out)
	}
	if !strings.Contains(err.Error(), "--verify reads a .sqf") {
		t.Errorf("error = %v", err)
	}
}

// A submitted job re-runs create from the name alone, so everything that changes
// what it builds or where it lands travels as a flag; submission and mode flags
// do not.
func TestCreateJobFlagsCarryWhatChangesTheBuild(t *testing.T) {
	saved := []any{createChannels, createSources, createLayer, createBlockSize, createDataBlockSize, createAlwaysSubmitData, createUpdate}
	t.Cleanup(func() {
		createChannels, createSources = saved[0].([]string), saved[1].([]string)
		createLayer, createBlockSize, createDataBlockSize = saved[2].(string), saved[3].(string), saved[4].(string)
		createAlwaysSubmitData, createUpdate = saved[5].(bool), saved[6].(bool)
	})
	createChannels = []string{"conda-forge", "bioconda"}
	createSources = []string{"lab"}
	createLayer, createBlockSize, createDataBlockSize = "u", "256k", "1m"
	createAlwaysSubmitData, createUpdate = true, true

	got := strings.Join(createJobFlags(), " ")
	want := "--channel conda-forge --channel bioconda --source lab --layer u --block-size 256k --data-block-size 1m"
	if got != want {
		t.Errorf("flags = %q, want %q", got, want)
	}
}
