package project

import (
	"os"
	"path/filepath"
	"testing"

	"github.com/Justype/condatainer/internal/helperhistory"
	"github.com/Justype/condatainer/internal/project/lock"
)

func TestClassifyOverlay(t *testing.T) {
	cases := []struct {
		name     string
		overlay  string
		location string
		wantKey  string
		wantOK   bool
	}{
		{"catalog name", "rstudio-server/4.4.3", ".", "rstudio-server/4.4.3", true},
		{"project-relative path at root", "overlays/combined.sqf", ".", lock.PathPrefix + "overlays/combined.sqf", true},
		{"project-relative path under a location", "combined.sqf", "steps1", lock.PathPrefix + "steps1/combined.sqf", true},
		{"writable img is never pinnable", "env.img", ".", "", false},
		{"absolute path is external", "/scratch/test.sqf", ".", "", false},
		{"escapes the root once joined", "../../scratch/test.sqf", ".", "", false},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			key, ok := classifyOverlay(c.overlay, c.location)
			if ok != c.wantOK {
				t.Fatalf("ok = %v, want %v", ok, c.wantOK)
			}
			if ok && key != c.wantKey {
				t.Fatalf("key = %q, want %q", key, c.wantKey)
			}
		})
	}
}

// UnpinnedHelperOverlays surfaces a recorded overlay with no matching pin,
// and drops it the moment a pin is added — never writing one itself.
func TestUnpinnedHelperOverlays(t *testing.T) {
	root := projectRoot(t)
	if err := helperhistory.RecordUsed(root, "rstudio-server", ".", []string{"rstudio-server/4.4.3"}); err != nil {
		t.Fatal(err)
	}

	l := lock.New()
	unpinned, err := UnpinnedHelperOverlays(root, l)
	if err != nil {
		t.Fatal(err)
	}
	if len(unpinned) != 1 || unpinned[0] != "rstudio-server/4.4.3" {
		t.Fatalf("unpinned = %v, want [rstudio-server/4.4.3]", unpinned)
	}

	l.Pins["rstudio-server/4.4.3"] = lock.PinEntry{Artifact: "provenance/ghost@000000000000"}
	unpinned, err = UnpinnedHelperOverlays(root, l)
	if err != nil {
		t.Fatal(err)
	}
	if len(unpinned) != 0 {
		t.Fatalf("unpinned = %v, want none once pinned", unpinned)
	}
}

// An overlay a `-o` path recorded but that escapes the project root is
// never suggested — nothing could pin it either direction.
func TestUnpinnedHelperOverlaysExcludesExternalPaths(t *testing.T) {
	root := projectRoot(t)
	if err := helperhistory.RecordUsed(root, "rstudio-server", ".", []string{"../../scratch/test.sqf"}); err != nil {
		t.Fatal(err)
	}
	unpinned, err := UnpinnedHelperOverlays(root, lock.New())
	if err != nil {
		t.Fatal(err)
	}
	if len(unpinned) != 0 {
		t.Fatalf("unpinned = %v, want none for an external path", unpinned)
	}
}

// ManualPinUsage names every helper recorded using a manual pin, and reports
// a fact ("no recorded helper usage") rather than nothing for one with none
// — an absent list is not the same as "not manual".
func TestManualPinUsage(t *testing.T) {
	root := projectRoot(t)
	if err := helperhistory.RecordUsed(root, "rstudio-server", ".", []string{"build-essential"}); err != nil {
		t.Fatal(err)
	}
	if err := helperhistory.RecordUsed(root, "jupyterlab", ".", []string{"build-essential"}); err != nil {
		t.Fatal(err)
	}

	l := lock.New()
	l.Pins["build-essential"] = lock.PinEntry{Artifact: "provenance/a@000000000000", Manual: true}
	l.Pins["unrelated/1.0"] = lock.PinEntry{Artifact: "provenance/b@000000000000", Manual: true}
	// A #DEP:-derived (non-manual) pin never appears here at all.
	l.Pins["scripted/1.0"] = lock.PinEntry{Artifact: "provenance/c@000000000000"}

	usage, err := ManualPinUsage(root, l)
	if err != nil {
		t.Fatal(err)
	}
	if _, ok := usage["scripted/1.0"]; ok {
		t.Fatalf("a non-manual pin must not appear in ManualPinUsage: %v", usage)
	}
	names := usage["build-essential"]
	if len(names) != 2 || names[0] != "jupyterlab" || names[1] != "rstudio-server" {
		t.Fatalf("build-essential usage = %v, want [jupyterlab rstudio-server]", names)
	}
	if got, ok := usage["unrelated/1.0"]; !ok || len(got) != 0 {
		t.Fatalf("unrelated/1.0 usage = %v, ok=%v, want present and empty (a fact, not an omission)", got, ok)
	}
}

// StatusAt never walks up to an ancestor project — the CLI's --project
// convention, unlike Status's ambient StandingAt walk-up.
func TestStatusAtDoesNotWalkUpToAnAncestor(t *testing.T) {
	root := projectRoot(t)
	if err := lock.Publish(root, lock.New()); err != nil {
		t.Fatal(err)
	}
	sub := filepath.Join(root, "sub")
	if err := os.MkdirAll(sub, 0o775); err != nil {
		t.Fatal(err)
	}

	if _, err := StatusAt(sub); err == nil {
		t.Fatal("StatusAt(sub) found the parent project; it must refuse instead")
	}

	// Status, by contrast, is the ambient/dashboard entry point and is
	// expected to walk up and find it.
	status, err := Status(sub)
	if err != nil {
		t.Fatal(err)
	}
	if !status.HasProject || status.Root != root {
		t.Fatalf("Status(sub) = %+v, want the ancestor project at %s", status, root)
	}
}

// StatusAt on the project's own root succeeds and reports it.
func TestStatusAtOnTheRootItself(t *testing.T) {
	root := projectRoot(t)
	if err := lock.Publish(root, lock.New()); err != nil {
		t.Fatal(err)
	}
	status, err := StatusAt(root)
	if err != nil {
		t.Fatal(err)
	}
	if !status.HasProject || status.Root != root {
		t.Fatalf("StatusAt(root) = %+v, want HasProject at %s", status, root)
	}
}
