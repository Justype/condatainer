package project

import (
	"context"
	"strings"
	"testing"

	"github.com/Justype/condatainer/internal/project/lock"
)

// Outside a project there is nothing to stand in, and the caller's ordinary
// fallback takes over.
func TestStandingAtNilOutsideAProject(t *testing.T) {
	root := t.TempDir()

	standing, err := StandingAt(root)
	if err != nil {
		t.Fatal(err)
	}
	if standing != nil {
		t.Fatalf("standing = %+v outside a project", standing)
	}
}

// Nothing to resolve is not a reason to look for a project.
func TestStandingResolveNamesEmptyWithNoNames(t *testing.T) {
	standing := &Standing{Root: projectRoot(t), Lock: lock.New()}

	if paths, err := standing.ResolveNames(context.Background(), nil); paths != nil || err != nil {
		t.Fatalf("paths = %v, err = %v, want nil and no error with no names", paths, err)
	}
}

// A required overlay the project has not pinned is a refusal naming
// `project restore`, never a fall back to whatever currently answers to the
// name.
func TestStandingResolveNamesRefusesAnUnpinnedName(t *testing.T) {
	root := projectRoot(t)
	if err := lock.Publish(root, lock.New()); err != nil {
		t.Fatal(err)
	}
	standing, err := StandingAt(root)
	if err != nil || standing == nil {
		t.Fatalf("StandingAt: standing = %v, err = %v", standing, err)
	}

	_, err = standing.ResolveNames(context.Background(), []string{"ubuntu24/build-essential"})
	if err == nil {
		t.Fatal("an unpinned required overlay was accepted")
	}
	if !strings.Contains(err.Error(), "project restore") {
		t.Errorf("error does not name the remedy: %v", err)
	}
}

// A pinned name whose artifact is not vendored is the same refusal, not a
// crash or a silent skip.
func TestStandingResolveNamesRefusesAMissingArtifact(t *testing.T) {
	root := projectRoot(t)
	l := lock.New()
	l.Pins["ubuntu24/build-essential"] = lock.PinEntry{Artifact: "provenance/ghost@000000000000"}
	if err := lock.Publish(root, l); err != nil {
		t.Fatal(err)
	}
	standing, err := StandingAt(root)
	if err != nil || standing == nil {
		t.Fatalf("StandingAt: standing = %v, err = %v", standing, err)
	}

	if _, err := standing.ResolveNames(context.Background(), []string{"ubuntu24/build-essential"}); err == nil {
		t.Fatal("a pin with no vendored artifact was accepted")
	}
}

// A lock with no base pin at all — one an older or hand-built lock could
// still have — falls through rather than erroring, matching a caller with
// its own ordinary default.
func TestStandingBaseEmptyWithNoBasePin(t *testing.T) {
	standing := &Standing{Root: projectRoot(t), Lock: lock.New()}

	path, err := standing.Base(context.Background())
	if path != "" || err != nil {
		t.Fatalf("path = %q, err = %v, want empty and no error with no base pin", path, err)
	}
}

// An unresolved base pin refuses naming `project restore`, never a silent
// fall back to this machine's configured default_distro.
func TestStandingBaseRefusesAnUnresolvedBase(t *testing.T) {
	root := projectRoot(t)
	l := lock.New()
	l.Pins[lock.BaseKey] = lock.PinEntry{Artifact: "provenance/ghost@000000000000"}
	if err := lock.Publish(root, l); err != nil {
		t.Fatal(err)
	}
	standing, err := StandingAt(root)
	if err != nil || standing == nil {
		t.Fatalf("StandingAt: standing = %v, err = %v", standing, err)
	}

	if _, err := standing.Base(context.Background()); err == nil {
		t.Fatal("an unresolved base pin was accepted")
	}
}

// With no base pin, SelectedDistro answers "" rather than erroring — the
// safe default for a completion or display path.
func TestStandingSelectedDistroEmptyWithNoBasePin(t *testing.T) {
	standing := &Standing{Root: projectRoot(t), Lock: lock.New()}

	if got := standing.SelectedDistro(); got != "" {
		t.Fatalf("SelectedDistro() = %q, want empty with no base pin", got)
	}
}
