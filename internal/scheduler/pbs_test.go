package scheduler

import (
	"fmt"
	"os"
	"path/filepath"
	"regexp"
	"strings"
	"testing"
)

func newTestPbsScheduler() *PbsScheduler {
	return &PbsScheduler{
		qsubBin:     "/usr/bin/qsub",
		directiveRe: regexp.MustCompile(`^#PBS\s+(.+)`),
		jobIDRe:     regexp.MustCompile(`^(\d+\..*|^\d+)$`),
	}
}

// Test that scheduler uses outputDir for log files
func TestPbsCreateScriptUsesOutputDir(t *testing.T) {
	tmpOutputDir := t.TempDir()

	pbs := newTestPbsScheduler()

	jobSpec := &JobSpec{
		Name:    "pbstest/job",
		Command: "echo 'hello'",
		Specs: &ScriptSpecs{
			RawFlags: []string{},
			Ncpus:    1,
		},
	}

	scriptPath, err := pbs.CreateScriptWithSpec(jobSpec, tmpOutputDir)
	if err != nil {
		t.Fatalf("CreateScriptWithSpec failed: %v", err)
	}

	// Ensure the generated script uses outputDir for logs
	logPath := filepath.Join(tmpOutputDir, "pbstest--job.log")
	content, err := os.ReadFile(scriptPath)
	if err != nil {
		t.Fatalf("Failed to read generated script: %v", err)
	}
	if !strings.Contains(string(content), fmt.Sprintf("-o %s", logPath)) {
		t.Errorf("Generated script does not contain expected output path %s\nScript:\n%s", logPath, string(content))
	}
}
