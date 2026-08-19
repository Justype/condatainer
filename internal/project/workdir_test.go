package project

import (
	"os"
	"path/filepath"
	"strings"
	"testing"
)

// Nothing sets a working directory today, so a job inherits the scheduler's —
// the submission directory under SLURM, `$HOME` under PBS. A project states the
// root instead of taking an answer that varies by scheduler.
func TestWorkDirDefaultsToTheProjectRoot(t *testing.T) {
	root := t.TempDir()
	got, err := WorkDir(root, "")
	if err != nil {
		t.Fatal(err)
	}
	if got != root {
		t.Fatalf("WorkDir = %q, want the project root %q", got, root)
	}
}

func TestWorkDirAcceptsTheRootDeclaredExplicitly(t *testing.T) {
	root := t.TempDir()
	for _, declared := range []string{root, root + "/", filepath.Join(root, "scripts", "..")} {
		got, err := WorkDir(root, declared)
		if err != nil {
			t.Errorf("WorkDir(%q) = %v, want the root accepted", declared, err)
			continue
		}
		if got != root {
			t.Errorf("WorkDir(%q) = %q, want %q", declared, got, root)
		}
	}
}

// A #SBATCH --chdir elsewhere splits the two anchors: overlays would resolve
// against the root while the script's own relative paths resolve against the
// declared directory.
func TestWorkDirRefusesADirectoryThatIsNotTheRoot(t *testing.T) {
	root := t.TempDir()
	elsewhere := t.TempDir()

	_, err := WorkDir(root, elsewhere)
	if err == nil {
		t.Fatal("a working directory outside the project was accepted")
	}
	if !strings.Contains(err.Error(), root) || !strings.Contains(err.Error(), elsewhere) {
		t.Errorf("error names neither directory: %v", err)
	}
}

// A subdirectory is not the root either: the script's relative paths would
// resolve one level down from everything the lock describes.
func TestWorkDirRefusesASubdirectoryOfTheRoot(t *testing.T) {
	root := t.TempDir()
	scripts := filepath.Join(root, "scripts")
	if err := os.MkdirAll(scripts, 0o775); err != nil {
		t.Fatal(err)
	}
	if _, err := WorkDir(root, scripts); err == nil {
		t.Fatal("a subdirectory was accepted as the working directory")
	}
}

// A symlinked root is the same root. Refusing it would fail a job for the shape
// of the path rather than for where it points.
func TestWorkDirFollowsSymlinksToTheSameDirectory(t *testing.T) {
	root := t.TempDir()
	link := filepath.Join(t.TempDir(), "link")
	if err := os.Symlink(root, link); err != nil {
		t.Skipf("symlinks unavailable: %v", err)
	}
	if _, err := WorkDir(root, link); err != nil {
		t.Fatalf("a symlink to the root was refused: %v", err)
	}
}
