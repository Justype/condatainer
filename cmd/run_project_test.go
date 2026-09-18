package cmd

import (
	"context"
	"errors"
	"path/filepath"
	"strings"
	"testing"

	"github.com/Justype/condatainer/internal/project/lock"
	"github.com/Justype/condatainer/internal/scheduler"
)

func specs() *scheduler.ScriptSpecs { return &scheduler.ScriptSpecs{} }

// A script outside any project keeps today's behaviour untouched.
func TestProjectRunContextIgnoresANonProjectScript(t *testing.T) {
	dir := t.TempDir()
	script := filepath.Join(dir, "run.sh")
	writeScript(t, dir, "run.sh", "#DEP: star/2.7.11b\nrun\n")
	t.Chdir(dir)

	got, err := projectRunContext(context.Background(), script, specs())
	if err != nil {
		t.Fatal(err)
	}
	if got != nil {
		t.Fatalf("context = %#v, want none outside a project", got)
	}
}

// The project is the current directory, so naming a project's script from
// outside resolves by name exactly as an unlocked script does. This is what
// replaced --no-project: `cd` elsewhere is the opt-out.
func TestProjectRunContextIgnoresAScriptNamedFromOutside(t *testing.T) {
	root := newProject(t)
	writeScript(t, root, "run.sh", "#DEP: star/2.7.11b\nrun\n")
	t.Chdir(t.TempDir())

	got, err := projectRunContext(context.Background(), filepath.Join(root, "run.sh"), specs())
	if err != nil {
		t.Fatal(err)
	}
	if got != nil {
		t.Fatalf("context = %#v, want none from outside the root", got)
	}
}

// Only the directory holding cnt-lock/ is the root; a subdirectory is not.
func TestProjectRunContextFindsAProjectFromASubdirectory(t *testing.T) {
	root := newProject(t)
	writeScript(t, root, "scripts/run.sh", "#DEP: env.img\nrun\n")
	if err := lock.Publish(root, lock.New()); err != nil {
		t.Fatal(err)
	}
	t.Chdir(filepath.Join(root, "scripts"))

	got, err := projectRunContext(context.Background(), "run.sh", specs())
	if err != nil {
		t.Fatal(err)
	}
	if got == nil {
		t.Fatal("a subdirectory of a project was not recognized as standing in it")
	}
	if got.Root != root {
		t.Fatalf("context.Root = %q, want %q", got.Root, root)
	}
}

// Inside a project an unsatisfiable declaration is fatal, and the message names
// the one remedy rather than leaving the reader to guess.
func TestProjectRunContextFailsWithTheRestoreRemedy(t *testing.T) {
	root := newProject(t)
	writeScript(t, root, "run.sh", "#DEP: star/2.7.11b\nrun\n")
	t.Chdir(root)

	_, err := projectRunContext(context.Background(), filepath.Join(root, "run.sh"), specs())
	if err == nil {
		t.Fatal("a project with nothing selected was allowed to run")
	}
	if !strings.Contains(err.Error(), "project restore") {
		t.Errorf("error does not name the remedy: %v", err)
	}
	if !strings.Contains(err.Error(), "star/2.7.11b") {
		t.Errorf("error does not name the declaration: %v", err)
	}
}

// An unpinnable declaration is mounted as the literal path it names, so a
// project carrying only those runs with no selections at all.
func TestProjectRunContextMountsUnpinnableDeclarations(t *testing.T) {
	root := newProject(t)
	writeScript(t, root, "run.sh", "#DEP: env.img\nrun\n")
	if err := lock.Publish(root, lock.New()); err != nil {
		t.Fatal(err)
	}
	t.Chdir(root)

	got, err := projectRunContext(context.Background(), filepath.Join(root, "run.sh"), specs())
	if err != nil {
		t.Fatal(err)
	}
	if got == nil {
		t.Fatal("no project context inside a project")
	}
	want := filepath.Join(root, "env.img")
	if len(got.Overlays) != 1 || got.Overlays[0] != want {
		t.Fatalf("overlays = %v, want %q anchored on the root", got.Overlays, want)
	}
}

// Nothing sets a working directory today, so a job takes the scheduler's —
// $HOME under PBS. A project states the root instead.
func TestProjectRunContextStatesTheRootAsTheWorkingDirectory(t *testing.T) {
	root := newProject(t)
	writeScript(t, root, "run.sh", "run\n")
	if err := lock.Publish(root, lock.New()); err != nil {
		t.Fatal(err)
	}

	t.Chdir(root)
	script := filepath.Join(root, "run.sh")
	scriptSpecs := specs()
	if _, err := projectRunContext(context.Background(), script, scriptSpecs); err != nil {
		t.Fatal(err)
	}
	if scriptSpecs.Control.WorkDir != mustEvalPath(t, root) && scriptSpecs.Control.WorkDir != root {
		t.Fatalf("WorkDir = %q, want the project root %q", scriptSpecs.Control.WorkDir, root)
	}
}

// A declared working directory elsewhere is refused, not overridden: it was
// written on purpose and only its author can resolve the conflict.
func TestProjectRunContextRefusesAForeignWorkingDirectory(t *testing.T) {
	root := newProject(t)
	writeScript(t, root, "run.sh", "run\n")
	if err := lock.Publish(root, lock.New()); err != nil {
		t.Fatal(err)
	}

	t.Chdir(root)
	elsewhere := t.TempDir()
	scriptSpecs := specs()
	scriptSpecs.Control.WorkDir = elsewhere

	_, err := projectRunContext(context.Background(), filepath.Join(root, "run.sh"), scriptSpecs)
	if err == nil {
		t.Fatal("a working directory outside the project was accepted")
	}
	if !strings.Contains(err.Error(), elsewhere) {
		t.Errorf("error does not name the declared directory: %v", err)
	}
}

// The script may live anywhere; the root is what declarations anchor on. Run
// from the root, `run scripts/align.sh` is an ordinary project run.
func TestProjectRunContextAnchorsASubdirectoryScriptOnTheRoot(t *testing.T) {
	root := newProject(t)
	writeScript(t, root, "scripts/align.sh", "#DEP: overlays/tool.img\nrun\n")
	if err := lock.Publish(root, lock.New()); err != nil {
		t.Fatal(err)
	}
	t.Chdir(root)

	got, err := projectRunContext(context.Background(), filepath.Join(root, "scripts", "align.sh"), specs())
	if err != nil {
		t.Fatal(err)
	}
	if got == nil {
		t.Fatal("a subdirectory script did not find the project above it")
	}
	want := filepath.Join(root, "overlays", "tool.img")
	if len(got.Overlays) != 1 || got.Overlays[0] != want {
		t.Fatalf("overlays = %v, want %q — the root, not the script directory", got.Overlays, want)
	}
}

func mustEvalPath(t *testing.T, path string) string {
	t.Helper()
	resolved, err := filepath.EvalSymlinks(path)
	if err != nil {
		return path
	}
	return resolved
}

func TestResolveDepsRefusesAVersionConstraint(t *testing.T) {
	dir := t.TempDir()
	writeScript(t, dir, "run.sh", "#DEP: samtools>=1.20\nrun\n")

	if _, err := resolveDeps(context.Background(), filepath.Join(dir, "run.sh"), "run.sh"); !errors.Is(err, errRunAborted) {
		t.Fatalf("err = %v, want the run to abort on a constrained #DEP:", err)
	}
}
