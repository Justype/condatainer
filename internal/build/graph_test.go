package build

import (
	"os"
	"path/filepath"
	"strings"
	"testing"
)

// The job re-runs create on a node, so every choice the user made on the command
// line that changes what is built has to be in that command.
func TestSchedulerCreateCommandCarriesTheBuildFlags(t *testing.T) {
	got := buildSchedulerCreateCommand([]string{"star/2.7.11b"}, nil, true, true, true, nil)
	for _, want := range []string{"--update", "--store", "--no-prebuilt", "star/2.7.11b"} {
		if !strings.Contains(got, want) {
			t.Errorf("command %q lacks %q", got, want)
		}
	}
	if plain := buildSchedulerCreateCommand([]string{"star/2.7.11b"}, nil, false, false, false, nil); plain != "condatainer create star/2.7.11b" {
		t.Errorf("a plain build renders %q", plain)
	}
}

// A flag value is quoted, so a channel or source name cannot break the job script.
func TestSchedulerCreateCommandQuotesFlagValues(t *testing.T) {
	got := buildSchedulerCreateCommand([]string{"star/2.7.11b"}, []string{"--channel", "my chan; rm -rf /"}, false, false, false, nil)
	want := "condatainer create --channel 'my chan; rm -rf /' star/2.7.11b"
	if got != want {
		t.Errorf("command = %q, want %q", got, want)
	}
}

// A store build is filed beside the installed name, so an installed name is no
// reason for the graph to skip it — that is what `create --store` is for.
func TestGraphDoesNotSkipAnInstalledNameForAStoreBuild(t *testing.T) {
	target := filepath.Join(t.TempDir(), "demo--1.sqf")
	if err := os.WriteFile(target, nil, 0o644); err != nil {
		t.Fatal(err)
	}
	plain := &BuildObject{spec: Spec{Image: ImageSpec{Name: "demo/1"}}, tgt: targetFor(target)}
	store := &BuildObject{spec: Spec{Image: ImageSpec{Name: "demo/1"}}, tgt: targetFor(target), storeOverflow: true}
	bg := &BuildGraph{}

	if !bg.installedAlready(plain) {
		t.Error("an installed name is not skipped")
	}
	if bg.installedAlready(store) {
		t.Error("an installed name is skipped for a store build")
	}
}
