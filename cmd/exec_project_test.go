package cmd

import (
	"path/filepath"
	"strings"
	"testing"

	"github.com/Justype/condatainer/internal/project/lock"
)

// Outside a project, `-o` means exactly what it always meant.
func TestProjectOverlaysPassesThroughOutsideAProject(t *testing.T) {
	t.Chdir(t.TempDir())

	got, err := projectOverlays([]string{"star/2.7.11b", "env.img"})
	if err != nil {
		t.Fatal(err)
	}
	if len(got) != 2 || got[0] != "star/2.7.11b" || got[1] != "env.img" {
		t.Fatalf("overlays = %v, want them untouched", got)
	}
}

// Standing in a project, a name that is not locked is refused rather than
// resolved to whatever currently answers to it.
func TestProjectOverlaysRefusesAnUnselectedName(t *testing.T) {
	root := newProject(t)
	if err := lock.Publish(root, lock.New()); err != nil {
		t.Fatal(err)
	}

	_, err := projectOverlays([]string{"star/2.7.11b"})
	if err == nil {
		t.Fatal("an unlocked name was accepted inside a project")
	}
	if !strings.Contains(err.Error(), "star/2.7.11b") {
		t.Errorf("error does not name the overlay: %v", err)
	}
	if !strings.Contains(err.Error(), "project pin") {
		t.Errorf("error does not name the remedy: %v", err)
	}
}

// A writable .img has no identity to pin, so it is mounted as written — and
// anchored on the root, which is where a project's relative paths resolve.
func TestProjectOverlaysMountsAWritableImageLiterally(t *testing.T) {
	root := newProject(t)
	if err := lock.Publish(root, lock.New()); err != nil {
		t.Fatal(err)
	}

	got, err := projectOverlays([]string{"env.img"})
	if err != nil {
		t.Fatal(err)
	}
	want := filepath.Join(root, "env.img")
	if len(got) != 1 || got[0] != want {
		t.Fatalf("overlays = %v, want %q", got, want)
	}
}

// The suffix is a mount mode, not part of what is being addressed, and it has
// to reach the runtime.
func TestProjectOverlaysKeepsTheMountMode(t *testing.T) {
	root := newProject(t)
	if err := lock.Publish(root, lock.New()); err != nil {
		t.Fatal(err)
	}

	got, err := projectOverlays([]string{"env.img:rw"})
	if err != nil {
		t.Fatal(err)
	}
	want := filepath.Join(root, "env.img") + ":rw"
	if len(got) != 1 || got[0] != want {
		t.Fatalf("overlays = %v, want %q", got, want)
	}
}

// A range is a build-recipe feature. Inside a project the lock says which
// version, so "any version in this range" is not a request.
func TestProjectOverlaysRefusesAVersionConstraint(t *testing.T) {
	root := newProject(t)
	if err := lock.Publish(root, lock.New()); err != nil {
		t.Fatal(err)
	}

	_, err := projectOverlays([]string{"star/2.7.11b>=2.7.0"})
	if err == nil {
		t.Fatal("a constrained overlay was accepted inside a project")
	}
	if !strings.Contains(err.Error(), "constraint") {
		t.Errorf("error does not explain the refusal: %v", err)
	}
}

// A project-relative .sqf is a restore output, so it answers to the lock like
// any other selection rather than being mounted on sight.
func TestProjectOverlaysRefusesAnUnselectedProjectPath(t *testing.T) {
	root := newProject(t)
	if err := lock.Publish(root, lock.New()); err != nil {
		t.Fatal(err)
	}
	writeScript(t, root, "overlays/tool.sqf", "not really an overlay\n")

	if _, err := projectOverlays([]string{"overlays/tool.sqf"}); err == nil {
		t.Fatal("an unlocked project path was mounted on sight")
	}
}

// Nothing to resolve is not a reason to look for a project.
func TestProjectOverlaysIgnoresAnEmptyList(t *testing.T) {
	newProject(t)

	got, err := projectOverlays(nil)
	if err != nil {
		t.Fatal(err)
	}
	if len(got) != 0 {
		t.Fatalf("overlays = %v, want none", got)
	}
}

func TestSplitOverlayModeSeparatesTheMountMode(t *testing.T) {
	for _, tc := range []struct{ in, value, suffix string }{
		{"env.img:rw", "env.img", ":rw"},
		{"tool.sqf:ro", "tool.sqf", ":ro"},
		{"star/2.7.11b", "star/2.7.11b", ""},
	} {
		value, suffix := splitOverlayMode(tc.in)
		if value != tc.value || suffix != tc.suffix {
			t.Errorf("splitOverlayMode(%q) = %q, %q; want %q, %q", tc.in, value, suffix, tc.value, tc.suffix)
		}
	}
}
