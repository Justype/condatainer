package cmd

import (
	"os"
	"path/filepath"
	"testing"
)

// Nothing to relocate to is a no-op: the process stays where it is.
func TestApplyProjectRelocationNoopWithNoDir(t *testing.T) {
	before, err := os.Getwd()
	if err != nil {
		t.Fatal(err)
	}

	if err := applyProjectRelocation(""); err != nil {
		t.Fatal(err)
	}
	after, err := os.Getwd()
	if err != nil {
		t.Fatal(err)
	}
	if after != before {
		t.Fatalf("cwd changed from %q to %q with no --project", before, after)
	}
}

// --project changes into the named directory.
func TestApplyProjectRelocationChangesDirectory(t *testing.T) {
	dir := t.TempDir()
	t.Chdir(t.TempDir())

	if err := applyProjectRelocation(dir); err != nil {
		t.Fatal(err)
	}
	got, err := os.Getwd()
	if err != nil {
		t.Fatal(err)
	}
	if got != dir {
		// macOS temp dirs resolve through a symlink, so compare the
		// symlink-resolved form rather than the raw string.
		resolved, _ := filepath.EvalSymlinks(dir)
		if got != resolved {
			t.Fatalf("cwd = %q, want %q", got, dir)
		}
	}
}

// --project and --no-project contradict each other.
func TestApplyProjectRelocationRefusesWithNoProject(t *testing.T) {
	noProjectRequested = true
	t.Cleanup(func() { noProjectRequested = false })

	if err := applyProjectRelocation(t.TempDir()); err == nil {
		t.Fatal("--project and --no-project were both accepted")
	}
}

// announceProject fires once per invocation: a second call is a no-op, so
// resolving overlays and the root in the same run does not print twice.
func TestAnnounceProjectFiresOnce(t *testing.T) {
	projectAnnounced = false
	t.Cleanup(func() { projectAnnounced = false })

	announceProject("/some/project")
	if !projectAnnounced {
		t.Fatal("announceProject did not mark itself announced")
	}
	announceProject("/some/project")
	if !projectAnnounced {
		t.Fatal("a second call cleared the announced flag")
	}
}
