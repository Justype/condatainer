package build

import (
	"os"
	"path/filepath"
	"testing"

	"github.com/Justype/condatainer/catalog"
	"github.com/Justype/condatainer/internal/config"
)

// always_submit_data sends a data build to the scheduler without directives; an
// app build is not sent, and a build made from packages or a file never is,
// because a job re-runs `create <name>` and could not reproduce it.
func TestAlwaysSubmitDataAppliesToNamedDataBuildsOnly(t *testing.T) {
	prevData, prevSubmit := config.Global.Build.AlwaysSubmitData, config.Global.SubmitJob
	config.Global.Build.AlwaysSubmitData, config.Global.SubmitJob = true, true
	t.Cleanup(func() { config.Global.Build.AlwaysSubmitData, config.Global.SubmitJob = prevData, prevSubmit })

	named := func(typ catalog.Type) *BuildObject {
		return &BuildObject{submitJob: true, spec: Spec{Image: ImageSpec{Name: "demo/1", Type: typ}}}
	}
	if !named(catalog.TypeData).RequiresScheduler() {
		t.Error("a named data build is not submitted under always_submit_data")
	}
	if named(catalog.TypeApp).RequiresScheduler() {
		t.Error("an app build is submitted under always_submit_data")
	}

}

// An external shell script with directives is submitted once the arguments that
// rebuild it are known, and never a definition: a job could not reproduce one.
func TestExternalScriptIsSubmittableOnlyWithItsJobArgs(t *testing.T) {
	prev := config.Global.SubmitJob
	config.Global.SubmitJob = true
	t.Cleanup(func() { config.Global.SubmitJob = prev })

	dir := t.TempDir()
	src := filepath.Join(dir, "demo.sh")
	if err := os.WriteFile(src, []byte("#!/usr/bin/env bash\n#DESC:demo\n#SBATCH --cpus-per-task=2\necho hi\n"), 0o644); err != nil {
		t.Fatal(err)
	}
	b, err := FromExternalSource(t.Context(), filepath.Join(dir, "demo"), src, false, dir, false)
	if err != nil {
		t.Fatal(err)
	}
	if b.RequiresScheduler() {
		t.Error("an external script is submitted before anything says how to rebuild it")
	}

	b.SetJobArgs([]string{"--prefix", filepath.Join(dir, "demo"), "--file", src})
	if !b.RequiresScheduler() {
		t.Error("an external script with directives is not submitted")
	}
	want := "condatainer create --prefix " + filepath.Join(dir, "demo") + " --file " + src
	if got := buildSchedulerCreateCommand(b.jobTarget(), nil, false, false, false, nil); got != want {
		t.Errorf("command = %q, want %q", got, want)
	}

	def := &BuildObject{buildType: BuildTypeDef}
	def.SetJobArgs([]string{"--file", "x.def"})
	if def.submitJob {
		t.Error("a definition was made submittable")
	}
}
