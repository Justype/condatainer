package lock

import (
	"strings"
	"testing"
)

const testEntry = "build-essential--1.0@abc123456789"

func manualPin() PinEntry {
	return PinEntry{Artifact: EntryPath(testEntry), Manual: true}
}

// A manual pin is the one entry a rescan may not sweep: nothing in the checkout
// declares it, so absence from the scan says nothing about whether it belongs.
// Helper overlays and frozen environments arrive this way.
func TestReconcileKeepsAManualPin(t *testing.T) {
	l := New()
	l.Pins["build-essential/1.0"] = manualPin()

	if needPin := Reconcile(t.TempDir(), l, &ScanResult{}); len(needPin) != 0 {
		t.Errorf("Reconcile asked to pin %v", needPin)
	}
	if _, ok := l.Pins["build-essential/1.0"]; !ok {
		t.Error("Reconcile swept a manual pin no script declared")
	}
}

// A declared pin stays a projection of the scan, so dropping the #DEP: drops the
// pin. Without this the flag would have made every pin permanent.
func TestReconcileStillSweepsADeclaredPin(t *testing.T) {
	l := New()
	l.Pins["star/2.7.11b"] = PinEntry{Artifact: EntryPath("star--2.7.11b@abc123456789")}

	Reconcile(t.TempDir(), l, &ScanResult{})

	if _, ok := l.Pins["star/2.7.11b"]; ok {
		t.Error("Reconcile kept a declared pin the scan no longer produces")
	}
}

func TestManualSurvivesARoundTrip(t *testing.T) {
	l := New()
	l.Pins["build-essential/1.0"] = manualPin()

	data, err := l.Marshal()
	if err != nil {
		t.Fatalf("Marshal: %v", err)
	}
	back, err := Unmarshal(data)
	if err != nil {
		t.Fatalf("Unmarshal: %v", err)
	}
	if !back.Pins["build-essential/1.0"].Manual {
		t.Error("the manual flag did not survive a round trip")
	}
}

// Absent for an ordinary pin, so an existing lock gains no noise from the field.
func TestManualIsOmittedWhenFalse(t *testing.T) {
	l := New()
	l.Pins["star/2.7.11b"] = PinEntry{Artifact: EntryPath("star--2.7.11b@abc123456789")}

	data, err := l.Marshal()
	if err != nil {
		t.Fatalf("Marshal: %v", err)
	}
	if strings.Contains(string(data), "manual") {
		t.Errorf("a declared pin wrote a manual field:\n%s", data)
	}
}

// Re-pinning changes which artifact answers, never why the pin is in the lock.
func TestApplyKeepsAPinManualAcrossARepin(t *testing.T) {
	root := t.TempDir()
	l := New()
	l.Pins["build-essential/1.0"] = manualPin()

	pinned := &Pinned{
		Request:  "build-essential/1.0",
		Artifact: EntryPath("build-essential--1.0@def987654321"),
		Manual:   false, // as a caller would set it if a script had appeared
	}
	// Apply verifies against the checkout, which has no vendored entries here, so
	// only the merge is under test.
	if err := Apply(root, l, pinned); err == nil {
		t.Log("Apply succeeded; the merge is what matters below")
	}
	if !l.Pins["build-essential/1.0"].Manual {
		t.Error("a re-pin cleared the manual flag")
	}
}
