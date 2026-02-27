package scheduler

import (
	"os"
	"path/filepath"
	"strings"
	"testing"
	"time"
)

func TestSlurmArrayJobScriptGeneration(t *testing.T) {
	tmpDir := t.TempDir()

	inputFile := filepath.Join(tmpDir, "input.txt")
	os.WriteFile(inputFile, []byte("sampleA\nsampleB extra_arg\nsampleC\n"), 0644)

	sched := newTestSlurmScheduler()

	jobSpec := &JobSpec{
		Name:    "testjob_20260226",
		Command: "condatainer run myscript.sh",
		Specs: &ScriptSpecs{
			Spec:    &ResourceSpec{CpusPerTask: 4, MemPerNodeMB: 8192, Time: time.Hour, Nodes: 1, TasksPerNode: 1},
			Control: RuntimeConfig{JobName: "testjob"},
		},
		Metadata: map[string]string{},
		Array: &ArraySpec{
			InputFile: inputFile,
			Count:     3,
			Limit:     2,
		},
	}

	scriptPath, err := sched.CreateScriptWithSpec(jobSpec, tmpDir)
	if err != nil {
		t.Fatalf("CreateScriptWithSpec failed: %v", err)
	}
	data, _ := os.ReadFile(scriptPath)
	script := string(data)

	checks := []struct{ label, want string }{
		{"array directive",         "#SBATCH --array=1-3%2"},
		{"output to /dev/null",     "#SBATCH --output=/dev/null"},
		{"error to /dev/null",      "#SBATCH --error=/dev/null"},
		{"ARRAY_IDX assignment",    "_ARRAY_IDX=$SLURM_ARRAY_TASK_ID"},
		{"sed extraction",          "sed -n \"${_ARRAY_IDX}p\" " + inputFile},
		{"ARRAY_ARGS export",     "export ARRAY_ARGS"},
		{"first-token extraction",  "_ARRAY_FIRST=${ARRAY_ARGS%% *}"},
		{"combined exec redirect",  "exec &>"},
		{"array job ID in header",  "SLURM_ARRAY_JOB_ID"},
		{"array index in metadata", "SLURM_ARRAY_TASK_ID"},
	}

	for _, c := range checks {
		if !strings.Contains(script, c.want) {
			t.Errorf("[%s] missing %q\nScript:\n%s", c.label, c.want, script)
		}
	}

	// Verify non-array behavior is unchanged
	jobSpecPlain := &JobSpec{
		Name:    "plainjob",
		Command: "condatainer run plain.sh",
		Specs: &ScriptSpecs{
			Spec:    &ResourceSpec{CpusPerTask: 2, MemPerNodeMB: 4096, Time: 30 * time.Minute, Nodes: 1, TasksPerNode: 1},
			Control: RuntimeConfig{JobName: "plainjob"},
		},
		Metadata: map[string]string{},
	}
	plainPath, err := sched.CreateScriptWithSpec(jobSpecPlain, tmpDir)
	if err != nil {
		t.Fatalf("plain CreateScriptWithSpec failed: %v", err)
	}
	plainData, _ := os.ReadFile(plainPath)
	plain := string(plainData)
	if strings.Contains(plain, "#SBATCH --array") {
		t.Error("plain job should not contain --array directive")
	}
	if strings.Contains(plain, "/dev/null") {
		t.Error("plain job should not have /dev/null output")
	}
	if !strings.Contains(plain, "plainjob.log") {
		t.Error("plain job should have named log file")
	}
}
