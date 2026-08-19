package cmd

import (
	"encoding/json"
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/Justype/condatainer/catalog"
	"github.com/Justype/condatainer/internal/artifact/capsule"
	"github.com/Justype/condatainer/internal/artifact/key"
	"github.com/Justype/condatainer/internal/artifact/meta"
	"github.com/Justype/condatainer/internal/project/lock"
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

// vendorArtifact puts a verifiable artifact into cnt-lock/artifacts/ and returns
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
	relative, err := lock.StageArtifact(root, capsule.EntryName(name, manifest.Keys.Identity.Digest()),
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

// A declared but unselected request leaves a valid partial lock and a nonzero
// exit: the lock is published, but never reported as complete.
func TestProjectLockPublishesPartiallyAndFails(t *testing.T) {
	root := newProject(t)
	writeScript(t, root, "run.sh", "#DEP: star/2.7.11b\nrun\n")

	if _, err := run(t, "project", "lock", "--json"); err == nil {
		t.Fatal("project lock succeeded with an unselected request")
	}
	if _, err := os.Stat(lock.FilePath(root)); err != nil {
		t.Fatalf("no lock was published: %v", err)
	}
	loaded, err := lock.Load(root)
	if err != nil {
		t.Fatalf("the published lock does not load: %v", err)
	}
	if len(loaded.Selections) != 0 {
		t.Fatalf("selections = %#v, want none invented", loaded.Selections)
	}
}

func TestProjectValidateSucceedsOnACompleteProject(t *testing.T) {
	root := newProject(t)
	writeScript(t, root, "run.sh", "#DEP: star/2.7.11b\nrun\n")
	artifact := vendorArtifact(t, root, "star/2.7.11b", "echo star\n")

	l := lock.New()
	l.Selections["star/2.7.11b"] = lock.Selection{Artifact: artifact}
	if err := lock.Publish(root, l); err != nil {
		t.Fatal(err)
	}

	if _, err := run(t, "project", "validate"); err != nil {
		t.Fatalf("validate failed on a complete project: %v", err)
	}
}

// An unpinnable declaration carrying the marker is not an unselected request.
func TestProjectValidateAcceptsAnUnpinnedDeclaration(t *testing.T) {
	root := newProject(t)
	writeScript(t, root, "run.sh", "#DEP: star/2.7.11b\n#DEP: env.img  ## unpinned — scratch\nrun\n")
	artifact := vendorArtifact(t, root, "star/2.7.11b", "echo star\n")

	l := lock.New()
	l.Selections["star/2.7.11b"] = lock.Selection{Artifact: artifact}
	if err := lock.Publish(root, l); err != nil {
		t.Fatal(err)
	}

	if _, err := run(t, "project", "validate"); err != nil {
		t.Fatalf("validate rejected a declared-unpinnable dependency: %v", err)
	}
}

// A #DEP: below the header is inert, and validate says so rather than passing.
func TestProjectValidateFailsOnAnInertDeclaration(t *testing.T) {
	root := newProject(t)
	writeScript(t, root, "run.sh", "#DEP: star/2.7.11b\nrun\n#DEP: cutadapt/5.0\n")
	artifact := vendorArtifact(t, root, "star/2.7.11b", "echo star\n")

	l := lock.New()
	l.Selections["star/2.7.11b"] = lock.Selection{Artifact: artifact}
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
		} `json:"unselected"`
	}
	if err := json.Unmarshal([]byte(out), &report); err != nil {
		t.Fatalf("output is not JSON: %v\n%s", err, out)
	}
	if report.Valid {
		t.Error("report claims valid")
	}
	if len(report.Unselected) != 1 || report.Unselected[0].Request != "star/2.7.11b" {
		t.Fatalf("unselected = %#v", report.Unselected)
	}
	if len(report.Unselected[0].Scripts) != 1 || report.Unselected[0].Scripts[0] != "run.sh" {
		t.Errorf("scripts = %v, want the declaring script", report.Unselected[0].Scripts)
	}
}

// Reconcile drops a selection nothing declares any more.
func TestProjectLockDropsAStaleSelection(t *testing.T) {
	root := newProject(t)
	artifact := vendorArtifact(t, root, "star/2.7.11b", "echo star\n")
	l := lock.New()
	l.Selections["star/2.7.11b"] = lock.Selection{Artifact: artifact}
	if err := lock.Publish(root, l); err != nil {
		t.Fatal(err)
	}
	writeScript(t, root, "run.sh", "echo no declarations\n")

	if _, err := run(t, "project", "lock"); err != nil {
		t.Fatalf("lock failed: %v", err)
	}
	loaded, err := lock.Load(root)
	if err != nil {
		t.Fatal(err)
	}
	if len(loaded.Selections) != 0 {
		t.Fatalf("a selection nothing declares survived: %#v", loaded.Selections)
	}
	if _, err := os.Stat(filepath.Join(lock.Dir(root), artifact)); !os.IsNotExist(err) {
		t.Errorf("its artifact was not pruned: %v", err)
	}
}

func TestProjectFlagSelectsTheNamedRoot(t *testing.T) {
	root := newProject(t)
	writeScript(t, root, "run.sh", "#DEP: star/2.7.11b\nrun\n")
	artifact := vendorArtifact(t, root, "star/2.7.11b", "echo star\n")
	l := lock.New()
	l.Selections["star/2.7.11b"] = lock.Selection{Artifact: artifact}
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
	l.Selections["star/2.7.11b"] = lock.Selection{Artifact: relative}
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
	l.Selections["star/2.7.11b"] = lock.Selection{Artifact: relative}
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
	l.Selections["star/2.7.11b"] = lock.Selection{Artifact: "artifacts/star--2.7.11b@000000000000"}
	if err := lock.Publish(root, l); err != nil {
		t.Fatal(err)
	}

	if out, err := run(t, "project", "restore", "--no-prebuilt"); err == nil {
		t.Fatalf("restore accepted an invalid lock: %s", out)
	}
}
