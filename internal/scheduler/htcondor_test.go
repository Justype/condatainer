package scheduler

import (
	"fmt"
	"os"
	"path/filepath"
	"regexp"
	"strings"
	"testing"
	"time"
)

// newTestHTCondorScheduler creates an HTCondor scheduler instance for testing
// without requiring condor_submit to be installed
func newTestHTCondorScheduler() *HTCondorScheduler {
	return &HTCondorScheduler{
		condorSubmitBin: "/usr/bin/condor_submit", // fake path for testing
		condorQBin:      "/usr/bin/condor_q",      // fake path for testing
		condorStatusBin: "/usr/bin/condor_status", // fake path for testing
		jobIDRe:         regexp.MustCompile(`submitted to cluster (\d+)`),
	}
}

func TestHTCondorEmailParsing(t *testing.T) {
	tests := []struct {
		name         string
		scriptLines  []string
		wantBegin    bool
		wantEnd      bool
		wantFail     bool
		wantMailUser string
	}{
		{
			name: "Always notification",
			scriptLines: []string{
				"# HTCondor Submit File",
				"notification = Always",
				"notify_user = test@example.com",
				"queue",
			},
			wantBegin:    true,
			wantEnd:      true,
			wantFail:     true,
			wantMailUser: "test@example.com",
		},
		{
			name: "Complete notification",
			scriptLines: []string{
				"# HTCondor Submit File",
				"notification = Complete",
				"notify_user = user@domain.org",
				"queue",
			},
			wantBegin:    false,
			wantEnd:      true,
			wantFail:     false,
			wantMailUser: "user@domain.org",
		},
		{
			name: "Error notification",
			scriptLines: []string{
				"# HTCondor Submit File",
				"notification = Error",
				"queue",
			},
			wantBegin:    false,
			wantEnd:      false,
			wantFail:     true,
			wantMailUser: "",
		},
		{
			name: "Never notification",
			scriptLines: []string{
				"# HTCondor Submit File",
				"notification = Never",
				"queue",
			},
			wantBegin: false,
			wantEnd:   false,
			wantFail:  false,
		},
		{
			name: "Case insensitive notification",
			scriptLines: []string{
				"# HTCondor Submit File",
				"notification = always",
				"notify_user = Test@Example.COM",
				"queue",
			},
			wantBegin:    true,
			wantEnd:      true,
			wantFail:     true,
			wantMailUser: "Test@Example.COM",
		},
		{
			name: "No notification directives",
			scriptLines: []string{
				"# HTCondor Submit File",
				"request_cpus = 4",
				"queue",
			},
			wantBegin: false,
			wantEnd:   false,
			wantFail:  false,
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			tmpDir := t.TempDir()
			scriptPath := filepath.Join(tmpDir, "test.sub")
			content := strings.Join(tt.scriptLines, "\n")
			if err := os.WriteFile(scriptPath, []byte(content), 0644); err != nil {
				t.Fatalf("Failed to create test submit file: %v", err)
			}

			htcondor := newTestHTCondorScheduler()
			specs, err := htcondor.ReadScriptSpecs(scriptPath)
			if err != nil {
				t.Fatalf("Failed to parse submit file: %v", err)
			}

			if specs.Control.EmailOnBegin != tt.wantBegin {
				t.Errorf("EmailOnBegin = %v; want %v", specs.Control.EmailOnBegin, tt.wantBegin)
			}
			if specs.Control.EmailOnEnd != tt.wantEnd {
				t.Errorf("EmailOnEnd = %v; want %v", specs.Control.EmailOnEnd, tt.wantEnd)
			}
			if specs.Control.EmailOnFail != tt.wantFail {
				t.Errorf("EmailOnFail = %v; want %v", specs.Control.EmailOnFail, tt.wantFail)
			}
			if specs.Control.MailUser != tt.wantMailUser {
				t.Errorf("MailUser = %q; want %q", specs.Control.MailUser, tt.wantMailUser)
			}
		})
	}
}

func TestHTCondorEmailScriptGeneration(t *testing.T) {
	tests := []struct {
		name             string
		emailOnBegin     bool
		emailOnEnd       bool
		emailOnFail      bool
		mailUser         string
		wantNotification string // expected notification value (empty = Never)
		wantNotifyUser   string
	}{
		{
			name:             "All three notifications (Always)",
			emailOnBegin:     true,
			emailOnEnd:       true,
			emailOnFail:      true,
			mailUser:         "test@example.com",
			wantNotification: "Always",
			wantNotifyUser:   "test@example.com",
		},
		{
			name:             "Only END (Complete)",
			emailOnBegin:     false,
			emailOnEnd:       true,
			emailOnFail:      false,
			mailUser:         "user@domain.org",
			wantNotification: "Complete",
			wantNotifyUser:   "user@domain.org",
		},
		{
			name:             "Only FAIL (Error)",
			emailOnBegin:     false,
			emailOnEnd:       false,
			emailOnFail:      true,
			mailUser:         "admin@company.com",
			wantNotification: "Error",
			wantNotifyUser:   "admin@company.com",
		},
		{
			name:             "No email notifications (Never)",
			emailOnBegin:     false,
			emailOnEnd:       false,
			emailOnFail:      false,
			wantNotification: "Never",
			wantNotifyUser:   "",
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			tmpDir := t.TempDir()

			htcondor := newTestHTCondorScheduler()

			jobSpec := &JobSpec{
				Name:    "test_job",
				Command: "echo 'test'",
				Specs: &ScriptSpecs{
					Spec: &ResourceSpec{
						Nodes:        1,
						TasksPerNode: 1,
						CpusPerTask:  4,
						MemPerNodeMB: 8000,
						Time:         time.Hour,
					},
					Control: RuntimeConfig{
						JobName:      "test_job",
						EmailOnBegin: tt.emailOnBegin,
						EmailOnEnd:   tt.emailOnEnd,
						EmailOnFail:  tt.emailOnFail,
						MailUser:     tt.mailUser,
					},
					RemainingFlags: []string{},
				},
			}

			subPath, err := htcondor.CreateScriptWithSpec(jobSpec, tmpDir)
			if err != nil {
				t.Fatalf("Failed to create script: %v", err)
			}

			// Read the generated submit file
			content, err := os.ReadFile(subPath)
			if err != nil {
				t.Fatalf("Failed to read generated submit file: %v", err)
			}
			subContent := string(content)

			// Check notification directive
			expectedLine := "notification = " + tt.wantNotification
			if !strings.Contains(subContent, expectedLine) {
				t.Errorf("Submit file missing expected line: %q\nSubmit content:\n%s", expectedLine, subContent)
			}

			// Check notify_user directive
			if tt.wantNotifyUser != "" {
				expectedLine := "notify_user = " + tt.wantNotifyUser
				if !strings.Contains(subContent, expectedLine) {
					t.Errorf("Submit file missing expected line: %q\nSubmit content:\n%s", expectedLine, subContent)
				}
			} else {
				if strings.Contains(subContent, "notify_user =") {
					t.Errorf("Submit file should not contain notify_user directive\nSubmit content:\n%s", subContent)
				}
			}
		})
	}
}

func TestHTCondorEmailRoundTrip(t *testing.T) {
	tmpDir := t.TempDir()

	// Create a submit file with notification directives
	originalScript := `# HTCondor Submit File
universe = vanilla
executable = job.sh
request_cpus = 8
request_memory = 16384
+MaxRuntime = 7200
notification = Always
notify_user = roundtrip@example.com
queue
`
	scriptPath := filepath.Join(tmpDir, "original.sub")
	if err := os.WriteFile(scriptPath, []byte(originalScript), 0644); err != nil {
		t.Fatalf("Failed to create test script: %v", err)
	}

	htcondor := newTestHTCondorScheduler()

	// Parse the original script
	specs, err := htcondor.ReadScriptSpecs(scriptPath)
	if err != nil {
		t.Fatalf("Failed to parse script: %v", err)
	}

	// Verify parsing
	if !specs.Control.EmailOnBegin || !specs.Control.EmailOnEnd || !specs.Control.EmailOnFail {
		t.Errorf("Parsing failed: EmailOnBegin=%v, EmailOnEnd=%v, EmailOnFail=%v",
			specs.Control.EmailOnBegin, specs.Control.EmailOnEnd, specs.Control.EmailOnFail)
	}
	if specs.Control.MailUser != "roundtrip@example.com" {
		t.Errorf("MailUser = %q; want %q", specs.Control.MailUser, "roundtrip@example.com")
	}
	if specs.Spec.Time != 2*time.Hour {
		t.Errorf("Time = %v; want 2h", specs.Spec.Time)
	}

	// Generate a new submit file from parsed specs
	jobSpec := &JobSpec{
		Name:    "roundtrip_test",
		Command: "echo 'Running job'",
		Specs:   specs,
	}

	subPath, err := htcondor.CreateScriptWithSpec(jobSpec, tmpDir)
	if err != nil {
		t.Fatalf("Failed to create script: %v", err)
	}

	// Read the generated submit file
	content, err := os.ReadFile(subPath)
	if err != nil {
		t.Fatalf("Failed to read generated submit file: %v", err)
	}
	subContent := string(content)

	// Verify the generated submit file contains the notification directives
	if !strings.Contains(subContent, "notification = Always") {
		t.Errorf("Generated submit file missing notification directive:\n%s", subContent)
	}
	if !strings.Contains(subContent, "notify_user = roundtrip@example.com") {
		t.Errorf("Generated submit file missing notify_user directive:\n%s", subContent)
	}
	if !strings.Contains(subContent, "+MaxRuntime = 7200") {
		t.Errorf("Generated submit file missing +MaxRuntime directive:\n%s", subContent)
	}

	// Verify NO duplicates
	notificationCount := strings.Count(subContent, "notification =")
	notifyUserCount := strings.Count(subContent, "notify_user =")
	maxRuntimeCount := strings.Count(subContent, "+MaxRuntime =")
	if notificationCount != 1 {
		t.Errorf("Expected exactly 1 notification directive, found %d:\n%s", notificationCount, subContent)
	}
	if notifyUserCount != 1 {
		t.Errorf("Expected exactly 1 notify_user directive, found %d:\n%s", notifyUserCount, subContent)
	}
	if maxRuntimeCount != 1 {
		t.Errorf("Expected exactly 1 +MaxRuntime directive, found %d:\n%s", maxRuntimeCount, subContent)
	}
}

func TestHTCondorCreateScriptUsesOutputDir(t *testing.T) {
	tmpOutputDir := t.TempDir()

	htcondor := newTestHTCondorScheduler()

	jobSpec := &JobSpec{
		Name:    "test/job",
		Command: "echo 'hello'",
		Specs: &ScriptSpecs{
			Spec: &ResourceSpec{
				Nodes:       1,
				CpusPerTask: 1,
			},
			RemainingFlags: []string{},
		},
	}

	subPath, err := htcondor.CreateScriptWithSpec(jobSpec, tmpOutputDir)
	if err != nil {
		t.Fatalf("CreateScriptWithSpec failed: %v", err)
	}

	// Ensure the generated submit file uses outputDir for logs
	logPath := filepath.Join(tmpOutputDir, "test--job.log")
	content, err := os.ReadFile(subPath)
	if err != nil {
		t.Fatalf("Failed to read generated submit file: %v", err)
	}
	if !strings.Contains(string(content), fmt.Sprintf("output = %s", logPath)) {
		t.Errorf("Generated submit file does not contain expected output path %s\nSubmit:\n%s", logPath, string(content))
	}

	// Verify wrapper script was also created
	shPath := filepath.Join(tmpOutputDir, "test--job.sh")
	if _, err := os.Stat(shPath); os.IsNotExist(err) {
		t.Errorf("Wrapper script was not created at %s", shPath)
	}
}

func TestHTCondorResourceParsing(t *testing.T) {
	tests := []struct {
		name        string
		scriptLines []string
		wantCpus    int
		wantMemMB   int64
		wantGpus    int
		wantTime    time.Duration
	}{
		{
			name: "All resources",
			scriptLines: []string{
				"# HTCondor Submit File",
				"executable = job.sh",
				"request_cpus = 8",
				"request_memory = 16384",
				"request_gpus = 2",
				"+MaxRuntime = 3600",
				"queue",
			},
			wantCpus:  8,
			wantMemMB: 16384,
			wantGpus:  2,
			wantTime:  time.Hour,
		},
		{
			name: "Memory with GB suffix",
			scriptLines: []string{
				"# HTCondor Submit File",
				"request_cpus = 4",
				"request_memory = 8GB",
				"queue",
			},
			wantCpus:  4,
			wantMemMB: 8192,
			wantGpus:  0,
			wantTime:  4 * time.Hour, // default
		},
		{
			name: "Memory with MB suffix",
			scriptLines: []string{
				"# HTCondor Submit File",
				"request_memory = 4096MB",
				"queue",
			},
			wantCpus:  2, // default
			wantMemMB: 4096,
			wantGpus:  0,
			wantTime:  4 * time.Hour, // default
		},
		{
			name: "Plain memory (default MB)",
			scriptLines: []string{
				"# HTCondor Submit File",
				"request_memory = 2048",
				"queue",
			},
			wantCpus:  2, // default
			wantMemMB: 2048,
			wantGpus:  0,
			wantTime:  4 * time.Hour, // default
		},
		{
			name: "Only CPUs",
			scriptLines: []string{
				"# HTCondor Submit File",
				"request_cpus = 16",
				"queue",
			},
			wantCpus:  16,
			wantMemMB: 8192, // default
			wantGpus:  0,
			wantTime:  4 * time.Hour, // default
		},
		{
			name: "Only GPUs",
			scriptLines: []string{
				"# HTCondor Submit File",
				"request_gpus = 4",
				"queue",
			},
			wantCpus:  2,    // default
			wantMemMB: 8192, // default
			wantGpus:  4,
			wantTime:  4 * time.Hour, // default
		},
		{
			name: "Time in seconds",
			scriptLines: []string{
				"# HTCondor Submit File",
				"+MaxRuntime = 86400",
				"queue",
			},
			wantCpus:  2,    // default
			wantMemMB: 8192, // default
			wantGpus:  0,
			wantTime:  24 * time.Hour,
		},
		{
			name: "No resource directives (defaults)",
			scriptLines: []string{
				"# HTCondor Submit File",
				"# Just a comment",
			},
			wantCpus:  2,    // default
			wantMemMB: 8192, // default
			wantGpus:  0,
			wantTime:  4 * time.Hour, // default
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			tmpDir := t.TempDir()
			scriptPath := filepath.Join(tmpDir, "test.sub")
			content := strings.Join(tt.scriptLines, "\n")
			if err := os.WriteFile(scriptPath, []byte(content), 0644); err != nil {
				t.Fatalf("Failed to create test submit file: %v", err)
			}

			htcondor := newTestHTCondorScheduler()
			specs, err := htcondor.ReadScriptSpecs(scriptPath)
			if err != nil {
				t.Fatalf("Failed to parse submit file: %v", err)
			}

			if specs.Spec.CpusPerTask != tt.wantCpus {
				t.Errorf("Ncpus = %d; want %d", specs.Spec.CpusPerTask, tt.wantCpus)
			}
			if specs.Spec.MemPerNodeMB != tt.wantMemMB {
				t.Errorf("MemMB = %d; want %d", specs.Spec.MemPerNodeMB, tt.wantMemMB)
			}
			if tt.wantGpus > 0 {
				if specs.Spec.Gpu == nil {
					t.Errorf("Gpu is nil; want count %d", tt.wantGpus)
				} else if specs.Spec.Gpu.Count != tt.wantGpus {
					t.Errorf("Gpu.Count = %d; want %d", specs.Spec.Gpu.Count, tt.wantGpus)
				}
			} else {
				if specs.Spec.Gpu != nil {
					t.Errorf("Gpu = %+v; want nil", specs.Spec.Gpu)
				}
			}
			if specs.Spec.Time != tt.wantTime {
				t.Errorf("Time = %v; want %v", specs.Spec.Time, tt.wantTime)
			}
		})
	}
}

func TestHTCondorGetJobResources(t *testing.T) {
	sched := &HTCondorScheduler{}

	t.Run("not in job", func(t *testing.T) {
		clearJobEnvVars(t)
		if res := sched.GetJobResources(); res != nil {
			t.Fatalf("expected nil, got %+v", res)
		}
	})

	t.Run("full resources", func(t *testing.T) {
		clearJobEnvVars(t)
		t.Setenv("_CONDOR_JOB_AD", "/tmp/condor_job_ad")
		t.Setenv("_CONDOR_REQUEST_CPUS", "16")
		t.Setenv("_CONDOR_REQUEST_MEMORY", "8192")
		t.Setenv("CUDA_VISIBLE_DEVICES", "0,1")

		res := sched.GetJobResources()
		if res == nil {
			t.Fatal("expected non-nil")
		}
		if res.Ncpus == nil || *res.Ncpus != 16 {
			t.Errorf("Ncpus = %v; want 16", res.Ncpus)
		}
		if res.MemMB == nil || *res.MemMB != 8192 {
			t.Errorf("MemMB = %v; want 8192", res.MemMB)
		}
		if res.Ngpus == nil || *res.Ngpus != 2 {
			t.Errorf("Ngpus = %v; want 2", res.Ngpus)
		}
	})

	t.Run("partial data", func(t *testing.T) {
		clearJobEnvVars(t)
		t.Setenv("_CONDOR_JOB_AD", "/tmp/condor_job_ad")
		t.Setenv("_CONDOR_REQUEST_CPUS", "4")

		res := sched.GetJobResources()
		if res == nil {
			t.Fatal("expected non-nil")
		}
		if res.Ncpus == nil || *res.Ncpus != 4 {
			t.Errorf("Ncpus = %v; want 4", res.Ncpus)
		}
		if res.MemMB != nil {
			t.Errorf("MemMB should be nil, got %d", *res.MemMB)
		}
		if res.Ngpus != nil {
			t.Errorf("Ngpus should be nil, got %d", *res.Ngpus)
		}
	})

	t.Run("invalid values", func(t *testing.T) {
		clearJobEnvVars(t)
		t.Setenv("_CONDOR_JOB_AD", "/tmp/condor_job_ad")
		t.Setenv("_CONDOR_REQUEST_CPUS", "not-a-number")
		t.Setenv("_CONDOR_REQUEST_MEMORY", "-100")

		res := sched.GetJobResources()
		if res == nil {
			t.Fatal("expected non-nil")
		}
		if res.Ncpus != nil {
			t.Errorf("Ncpus should be nil for invalid value, got %d", *res.Ncpus)
		}
		if res.MemMB != nil {
			t.Errorf("MemMB should be nil for negative value, got %d", *res.MemMB)
		}
	})
}

func TestHTCondorSubmitFileFormat(t *testing.T) {
	tmpDir := t.TempDir()

	htcondor := newTestHTCondorScheduler()

	jobSpec := &JobSpec{
		Name:    "format_test",
		Command: "echo 'hello world'",
		Specs: &ScriptSpecs{
			Spec: &ResourceSpec{
				Nodes:        1,
				TasksPerNode: 1,
				CpusPerTask:  8,
				MemPerNodeMB: 16384,
				Time:         2 * time.Hour,
			},
			Control: RuntimeConfig{
				JobName: "format_test",
			},
			RemainingFlags: []string{
				"request_cpus = 8",
				"request_memory = 16384",
			},
		},
	}

	subPath, err := htcondor.CreateScriptWithSpec(jobSpec, tmpDir)
	if err != nil {
		t.Fatalf("Failed to create script: %v", err)
	}

	content, err := os.ReadFile(subPath)
	if err != nil {
		t.Fatalf("Failed to read submit file: %v", err)
	}
	subContent := string(content)

	// Verify submit file structure
	requiredLines := []string{
		"universe = vanilla",
		"executable =",
		"transfer_executable = false",
		"+MaxRuntime = 7200",
		"notification = Never",
		"queue",
	}

	for _, line := range requiredLines {
		if !strings.Contains(subContent, line) {
			t.Errorf("Submit file missing required line: %q\nContent:\n%s", line, subContent)
		}
	}

	// Verify the wrapper script exists and is executable
	shPath := filepath.Join(tmpDir, "format_test.sh")
	info, err := os.Stat(shPath)
	if err != nil {
		t.Fatalf("Wrapper script not found: %v", err)
	}
	if info.Mode()&0111 == 0 {
		t.Errorf("Wrapper script is not executable: mode %v", info.Mode())
	}

	// Read wrapper script and verify it contains the command
	shContent, err := os.ReadFile(shPath)
	if err != nil {
		t.Fatalf("Failed to read wrapper script: %v", err)
	}
	if !strings.Contains(string(shContent), "echo 'hello world'") {
		t.Errorf("Wrapper script missing command:\n%s", string(shContent))
	}
	if !strings.Contains(string(shContent), "#!/bin/bash") {
		t.Errorf("Wrapper script missing shebang:\n%s", string(shContent))
	}
}

func TestHTCondorMemoryParsing(t *testing.T) {
	tests := []struct {
		name    string
		input   string
		wantMB  int64
		wantErr bool
	}{
		{"Plain number (default MB)", "4096", 4096, false},
		{"MB suffix", "8192MB", 8192, false},
		{"GB suffix", "16GB", 16384, false},
		{"G suffix", "8G", 8192, false},
		{"KB suffix", "1048576KB", 1024, false},
		{"TB suffix", "1TB", 1048576, false},
		{"Empty string", "", 0, true},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			got, err := parseHTCondorMemory(tt.input)
			if (err != nil) != tt.wantErr {
				t.Errorf("parseHTCondorMemory(%q) error = %v, wantErr %v", tt.input, err, tt.wantErr)
				return
			}
			if !tt.wantErr && got != tt.wantMB {
				t.Errorf("parseHTCondorMemory(%q) = %d; want %d", tt.input, got, tt.wantMB)
			}
		})
	}
}

func TestHTCondorTimeParsing(t *testing.T) {
	tests := []struct {
		name     string
		seconds  string
		wantTime time.Duration
	}{
		{"1 hour", "3600", time.Hour},
		{"24 hours", "86400", 24 * time.Hour},
		{"30 minutes", "1800", 30 * time.Minute},
		{"1 second", "1", time.Second},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			tmpDir := t.TempDir()
			scriptPath := filepath.Join(tmpDir, "test.sub")
			content := fmt.Sprintf("# HTCondor Submit File\n+MaxRuntime = %s\nqueue\n", tt.seconds)
			if err := os.WriteFile(scriptPath, []byte(content), 0644); err != nil {
				t.Fatalf("Failed to create test submit file: %v", err)
			}

			htcondor := newTestHTCondorScheduler()
			specs, err := htcondor.ReadScriptSpecs(scriptPath)
			if err != nil {
				t.Fatalf("Failed to parse submit file: %v", err)
			}
			if specs.Spec.Time != tt.wantTime {
				t.Errorf("Time = %v; want %v", specs.Spec.Time, tt.wantTime)
			}
		})
	}
}

func TestHTCondorDefaultNodesTasks(t *testing.T) {
	tests := []struct {
		name  string
		lines []string
	}{
		{
			name: "with resources",
			lines: []string{
				"# HTCondor Submit File",
				"request_cpus = 8",
				"request_memory = 16384",
				"queue",
			},
		},
		{
			name: "no directives",
			lines: []string{
				"# Empty submit file",
			},
		},
		{
			name: "full job",
			lines: []string{
				"# HTCondor Submit File",
				"executable = job.sh",
				"request_cpus = 16",
				"request_memory = 32768",
				"request_gpus = 2",
				"+MaxRuntime = 3600",
				"queue",
			},
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			tmpDir := t.TempDir()
			scriptPath := filepath.Join(tmpDir, "test.sub")
			content := strings.Join(tt.lines, "\n")
			if err := os.WriteFile(scriptPath, []byte(content), 0644); err != nil {
				t.Fatalf("Failed to create test submit file: %v", err)
			}

			htcondor := newTestHTCondorScheduler()
			specs, err := htcondor.ReadScriptSpecs(scriptPath)
			if err != nil {
				t.Fatalf("Failed to parse submit file: %v", err)
			}

			// HTCondor is inherently single-node, so Nodes and Ntasks should always be 1
			if specs.Spec.Nodes != 1 {
				t.Errorf("Nodes = %d; want 1", specs.Spec.Nodes)
			}
			if specs.Spec.TasksPerNode != 1 {
				t.Errorf("Ntasks = %d; want 1", specs.Spec.TasksPerNode)
			}
		})
	}
}

func TestHTCondorIsAvailable(t *testing.T) {
	t.Run("not in job", func(t *testing.T) {
		clearJobEnvVars(t)
		htcondor := newTestHTCondorScheduler()
		if !htcondor.IsAvailable() {
			t.Error("Expected IsAvailable to return true when not in a job")
		}
	})

	t.Run("inside job", func(t *testing.T) {
		clearJobEnvVars(t)
		t.Setenv("_CONDOR_JOB_AD", "/tmp/condor_job_ad")
		htcondor := newTestHTCondorScheduler()
		if htcondor.IsAvailable() {
			t.Error("Expected IsAvailable to return false when inside a job")
		}
	})

	t.Run("no binary", func(t *testing.T) {
		clearJobEnvVars(t)
		htcondor := &HTCondorScheduler{}
		if htcondor.IsAvailable() {
			t.Error("Expected IsAvailable to return false when no binary is set")
		}
	})
}

func TestHTCondorGetInfo(t *testing.T) {
	clearJobEnvVars(t)
	htcondor := newTestHTCondorScheduler()

	info := htcondor.GetInfo()
	if info.Type != "HTCondor" {
		t.Errorf("Type = %q; want %q", info.Type, "HTCondor")
	}
	if info.Binary != "/usr/bin/condor_submit" {
		t.Errorf("Binary = %q; want %q", info.Binary, "/usr/bin/condor_submit")
	}
	if info.InJob {
		t.Error("InJob should be false")
	}
	if !info.Available {
		t.Error("Available should be true")
	}
}

func TestTryParseHTCondorScript(t *testing.T) {
	tmpDir := t.TempDir()

	script := `# HTCondor Submit File
universe = vanilla
executable = job.sh
request_cpus = 4
request_memory = 8192
request_gpus = 1
+MaxRuntime = 7200
notification = Complete
notify_user = user@example.com
queue
`
	scriptPath := filepath.Join(tmpDir, "test.sub")
	if err := os.WriteFile(scriptPath, []byte(script), 0644); err != nil {
		t.Fatalf("Failed to create test submit file: %v", err)
	}

	specs, err := TryParseHTCondorScript(scriptPath)
	if err != nil {
		t.Fatalf("TryParseHTCondorScript failed: %v", err)
	}

	if specs.ScriptPath != "job.sh" {
		t.Errorf("ScriptPath = %q; want %q", specs.ScriptPath, "job.sh")
	}
	if specs.Spec.CpusPerTask != 4 {
		t.Errorf("Ncpus = %d; want 4", specs.Spec.CpusPerTask)
	}
	if specs.Spec.MemPerNodeMB != 8192 {
		t.Errorf("MemMB = %d; want 8192", specs.Spec.MemPerNodeMB)
	}
	if specs.Spec.Gpu == nil || specs.Spec.Gpu.Count != 1 {
		t.Errorf("Gpu = %+v; want count 1", specs.Spec.Gpu)
	}
	if specs.Spec.Time != 2*time.Hour {
		t.Errorf("Time = %v; want 2h", specs.Spec.Time)
	}
	if !specs.Control.EmailOnEnd {
		t.Error("EmailOnEnd should be true for Complete notification")
	}
	if specs.Control.MailUser != "user@example.com" {
		t.Errorf("MailUser = %q; want %q", specs.Control.MailUser, "user@example.com")
	}
	// RemainingFlags should be empty - all flags above are recognized and parsed into typed fields
	if len(specs.RemainingFlags) != 0 {
		t.Errorf("RemainingFlags count = %d; want 0 (all flags were recognized)", len(specs.RemainingFlags))
	}
}

// TestHTCondorFalsePositiveBashScript verifies that a plain bash build script containing
// bash variable assignments (e.g. "gencode_version=M6") is NOT misidentified as an
// HTCondor submit file. HasDirectives must be false so ParseScriptAny does not return
// the script as HTCondor type.
func TestHTCondorFalsePositiveBashScript(t *testing.T) {
	tmpDir := t.TempDir()

	// Simulate a typical CondaTainer build script with bash key=value assignments
	// but no valid HTCondor directives.
	script := `#!/usr/bin/bash
#DEP:samtools/1.22.1
#WHATIS:GRCm39 GENCODE M6 transcript FASTA
#ENV:TRANSCRIPT_FASTA=$app_root/gencode.vM6.transcripts.fa

gencode_version=M6
reference_index=/opt/references/index

install_app() {
    cd "$target_dir"
    wget -nv -O "gencode.v${gencode_version}.transcripts.fa.gz" "https://example.com"
}
`
	// Write with a non-.sub extension (as build scripts always are)
	scriptPath := filepath.Join(tmpDir, "transcript-gencode-M6")
	if err := os.WriteFile(scriptPath, []byte(script), 0644); err != nil {
		t.Fatalf("Failed to create test script: %v", err)
	}

	sched := newTestHTCondorScheduler()
	specs, err := sched.ReadScriptSpecs(scriptPath)
	if err != nil {
		t.Fatalf("ReadScriptSpecs failed: %v", err)
	}

	if specs.HasDirectives {
		t.Errorf("HasDirectives = true for plain bash script; bash variable assignments must not be treated as HTCondor directives (RawFlags: %v)", specs.RawFlags)
	}
	if len(specs.RawFlags) != 0 {
		t.Errorf("RawFlags = %v; want empty for plain bash script", specs.RawFlags)
	}
}
