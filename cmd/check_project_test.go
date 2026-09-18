package cmd

import (
	"context"
	"path/filepath"
	"strings"
	"testing"

	"github.com/Justype/condatainer/internal/project/lock"
)

// A directory with no cnt-lock/ keeps today's behaviour, -a included.
func TestProjectCheckIgnoresNonProjectScripts(t *testing.T) {
	dir := t.TempDir()
	writeScript(t, dir, "run.sh", "#DEP: star/2.7.11b\nrun\n")
	t.Chdir(dir)

	handled, err := projectCheck(context.Background(), []string{filepath.Join(dir, "run.sh")}, nil)
	if err != nil {
		t.Fatal(err)
	}
	if handled {
		t.Fatal("a non-project directory was handled by the project guard")
	}
}

// The project is the current directory and nothing above it. Naming a project's
// script from outside is an ordinary by-name check — this is what replaced
// --no-project, so an ad-hoc invocation needs no flag.
func TestProjectCheckIgnoresAProjectScriptNamedFromOutside(t *testing.T) {
	root := newProject(t)
	writeScript(t, root, "run.sh", "#DEP: star/2.7.11b\nrun\n")
	t.Chdir(t.TempDir())

	handled, err := projectCheck(context.Background(), []string{filepath.Join(root, "run.sh")}, nil)
	if err != nil {
		t.Fatal(err)
	}
	if handled {
		t.Fatal("a script named from outside its project went through the project guard")
	}
}

// Standing in a subdirectory of a project is standing in the project: the
// search for cnt-lock/ walks upward, so only the top of the tree marks a
// boundary, not the directory a script happens to sit in.
func TestProjectCheckFindsAProjectFromASubdirectory(t *testing.T) {
	root := newProject(t)
	writeScript(t, root, "scripts/run.sh", "#DEP: env.img\nrun\n")
	if err := lock.Publish(root, lock.New()); err != nil {
		t.Fatal(err)
	}
	t.Chdir(filepath.Join(root, "scripts"))

	handled, err := projectCheck(context.Background(), []string{"run.sh"}, nil)
	if !handled {
		t.Fatal("a subdirectory of a project was not recognized as standing in it")
	}
	if err != nil {
		t.Fatalf("a resolvable script was reported unrunnable: %v", err)
	}
}

// -a installs by name into a shared images directory, which is one checkout's
// dependency replacing what every other user of that directory sees.
func TestProjectCheckRefusesAutoInstallInsideAProject(t *testing.T) {
	root := newProject(t)
	writeScript(t, root, "run.sh", "#DEP: star/2.7.11b\nrun\n")
	t.Chdir(root)

	checkAutoInstall = true
	t.Cleanup(func() { checkAutoInstall = false })

	handled, err := projectCheck(context.Background(), []string{"run.sh"}, nil)
	if !handled || err == nil {
		t.Fatalf("-a was allowed inside a project (handled=%v, err=%v)", handled, err)
	}
	if !strings.Contains(err.Error(), "project restore") {
		t.Errorf("refusal does not name the remedy: %v", err)
	}
}

// The refusal holds even with nothing selected: the answer is to lock them, not
// to install them by name.
func TestProjectCheckRefusesAutoInstallWithAnEmptyLock(t *testing.T) {
	root := newProject(t)
	writeScript(t, root, "run.sh", "#DEP: star/2.7.11b\nrun\n")
	if err := lock.Publish(root, lock.New()); err != nil {
		t.Fatal(err)
	}
	t.Chdir(root)

	checkAutoInstall = true
	t.Cleanup(func() { checkAutoInstall = false })

	if _, err := projectCheck(context.Background(), []string{"run.sh"}, nil); err == nil {
		t.Fatal("-a was allowed against an empty lock")
	}
}

// check merges declarations and answers once; in a project the answer is per
// script, and the merged question is `project restore --dry-run`.
func TestProjectCheckRefusesSeveralScriptsInOneProject(t *testing.T) {
	root := newProject(t)
	writeScript(t, root, "a.sh", "#DEP: star/2.7.11b\nrun\n")
	writeScript(t, root, "b.sh", "#DEP: cutadapt/5.0\nrun\n")
	t.Chdir(root)

	_, err := projectCheck(context.Background(), []string{"a.sh", "b.sh"}, nil)
	if err == nil {
		t.Fatal("two scripts in one project were checked together")
	}
	if !strings.Contains(err.Error(), "--dry-run") {
		t.Errorf("refusal does not point at the whole-project question: %v", err)
	}
}

// A `<name>` addresses a recipe and belongs to no project, so it cannot be
// merged with a project script.
func TestProjectCheckRefusesMixingAProjectScriptWithANameArgument(t *testing.T) {
	root := newProject(t)
	writeScript(t, root, "run.sh", "#DEP: star/2.7.11b\nrun\n")
	t.Chdir(root)

	_, err := projectCheck(context.Background(), []string{"run.sh"}, []string{"samtools/1.22"})
	if err == nil {
		t.Fatal("a project script was merged with a name argument")
	}
	if !strings.Contains(err.Error(), "belongs to no project") {
		t.Errorf("refusal does not explain why: %v", err)
	}
}

// A bare `<name>` is a recipe wherever it is typed, so standing in a project
// must not capture it.
func TestProjectCheckIgnoresANameArgumentInsideAProject(t *testing.T) {
	root := newProject(t)
	t.Chdir(root)

	handled, err := projectCheck(context.Background(), nil, []string{"samtools/1.22"})
	if err != nil {
		t.Fatal(err)
	}
	if handled {
		t.Fatal("a `<name>` was captured by the project guard")
	}
}

// A single script whose declarations all resolve reports that it can run.
func TestProjectCheckReportsARunnableScript(t *testing.T) {
	root := newProject(t)
	writeScript(t, root, "run.sh", "#DEP: env.img\nrun\n")
	if err := lock.Publish(root, lock.New()); err != nil {
		t.Fatal(err)
	}
	t.Chdir(root)

	handled, err := projectCheck(context.Background(), []string{"run.sh"}, nil)
	if !handled {
		t.Fatal("a project script was not handled")
	}
	if err != nil {
		t.Fatalf("a resolvable script was reported unrunnable: %v", err)
	}
}

// An needPin declaration fails, and the hint names restore rather than -a.
func TestProjectCheckFailsOnAnUnselectedDeclaration(t *testing.T) {
	root := newProject(t)
	writeScript(t, root, "run.sh", "#DEP: star/2.7.11b\nrun\n")
	if err := lock.Publish(root, lock.New()); err != nil {
		t.Fatal(err)
	}
	t.Chdir(root)

	_, err := projectCheck(context.Background(), []string{"run.sh"}, nil)
	if err == nil {
		t.Fatal("an needPin declaration was reported runnable")
	}
}
