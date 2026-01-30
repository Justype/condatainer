package scheduler

import (
	"fmt"
	"os"
	"path/filepath"
	"regexp"
	"strings"
	"testing"

	"github.com/Justype/condatainer/internal/config"
)

func newTestPbsScheduler() *PbsScheduler {
	return &PbsScheduler{
		qsubBin:     "/usr/bin/qsub",
		directiveRe: regexp.MustCompile(`^#PBS\s+(.+)`),
		jobIDRe:     regexp.MustCompile(`^(\d+\..*|^\d+)$`),
	}
}

func TestPbsCreateScriptUsesConfigLogsDir(t *testing.T) {
	tmpScriptDir := t.TempDir()
	tmpLogsDir := t.TempDir()

	// Save and restore original value
	orig := config.Global.LogsDir
	defer func() { config.Global.LogsDir = orig }()

	config.Global.LogsDir = tmpLogsDir

	pbs := newTestPbsScheduler()

	jobSpec := &JobSpec{
		Name:    "pbstest/job",
		Command: "echo 'hello'",
		Specs: &ScriptSpecs{
			RawFlags: []string{},
			Ncpus:    1,
		},
	}

	scriptPath, err := pbs.CreateScriptWithSpec(jobSpec, tmpScriptDir)
	if err != nil {
		t.Fatalf("CreateScriptWithSpec failed: %v", err)
	}

	// Ensure the log directory was created in config.Global.LogsDir
	if info, err := os.Stat(tmpLogsDir); err != nil || !info.IsDir() {
		t.Fatalf("Expected logs directory %s to exist and be a directory", tmpLogsDir)
	}

	// Ensure the generated script uses the absolute path
	logPath := filepath.Join(tmpLogsDir, "pbstest--job.log")
	content, err := os.ReadFile(scriptPath)
	if err != nil {
		t.Fatalf("Failed to read generated script: %v", err)
	}
	if !strings.Contains(string(content), fmt.Sprintf("-o %s", logPath)) {
		t.Errorf("Generated script does not contain expected output path %s\nScript:\n%s", logPath, string(content))
	}
}
