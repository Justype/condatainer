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
		condorSubmitBin:  "/usr/bin/condor_submit",  // fake path for testing
		condorQBin:       "/usr/bin/condor_q",       // fake path for testing
		condorStatusBin:  "/usr/bin/condor_status",  // fake path for testing
		directiveRe:      regexp.MustCompile(`^\s*#CONDOR\s+(.+)$`),
		jobIDRe:          regexp.MustCompile(`submitted to cluster (\d+)`),
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
				"#!/bin/bash",
				"#CONDOR notification = Always",
				"#CONDOR notify_user = test@example.com",
			},
			wantBegin:    true,
			wantEnd:      true,
			wantFail:     true,
			wantMailUser: "test@example.com",
		},
		{
			name: "Complete notification",
			scriptLines: []string{
				"#!/bin/bash",
				"#CONDOR notification = Complete",
				"#CONDOR notify_user = user@domain.org",
			},
			wantBegin:    false,
			wantEnd:      true,
			wantFail:     false,
			wantMailUser: "user@domain.org",
		},
		{
			name: "Error notification",
			scriptLines: []string{
				"#!/bin/bash",
				"#CONDOR notification = Error",
			},
			wantBegin:    false,
			wantEnd:      false,
			wantFail:     true,
			wantMailUser: "",
		},
		{
			name: "Never notification",
			scriptLines: []string{
				"#!/bin/bash",
				"#CONDOR notification = Never",
			},
			wantBegin: false,
			wantEnd:   false,
			wantFail:  false,
		},
		{
			name: "Case insensitive notification",
			scriptLines: []string{
				"#!/bin/bash",
				"#CONDOR notification = always",
				"#CONDOR notify_user = Test@Example.COM",
			},
			wantBegin:    true,
			wantEnd:      true,
			wantFail:     true,
			wantMailUser: "Test@Example.COM",
		},
		{
			name: "No notification directives",
			scriptLines: []string{
				"#!/bin/bash",
				"#CONDOR request_cpus = 4",
			},
			wantBegin: false,
			wantEnd:   false,
			wantFail:  false,
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			tmpDir := t.TempDir()
			scriptPath := filepath.Join(tmpDir, "test.sh")
			content := strings.Join(tt.scriptLines, "\n")
			if err := os.WriteFile(scriptPath, []byte(content), 0644); err != nil {
				t.Fatalf("Failed to create test script: %v", err)
			}

			htcondor := newTestHTCondorScheduler()
			specs, err := htcondor.ReadScriptSpecs(scriptPath)
			if err != nil {
				t.Fatalf("Failed to parse script: %v", err)
			}

			if specs.EmailOnBegin != tt.wantBegin {
				t.Errorf("EmailOnBegin = %v; want %v", specs.EmailOnBegin, tt.wantBegin)
			}
			if specs.EmailOnEnd != tt.wantEnd {
				t.Errorf("EmailOnEnd = %v; want %v", specs.EmailOnEnd, tt.wantEnd)
			}
			if specs.EmailOnFail != tt.wantFail {
				t.Errorf("EmailOnFail = %v; want %v", specs.EmailOnFail, tt.wantFail)
			}
			if specs.MailUser != tt.wantMailUser {
				t.Errorf("MailUser = %q; want %q", specs.MailUser, tt.wantMailUser)
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
					JobName:      "test_job",
					Ncpus:        4,
					MemMB:        8000,
					Time:         time.Hour,
					EmailOnBegin: tt.emailOnBegin,
					EmailOnEnd:   tt.emailOnEnd,
					EmailOnFail:  tt.emailOnFail,
					MailUser:     tt.mailUser,
					RawFlags:     []string{},
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

	// Create a script with notification directives
	originalScript := `#!/bin/bash
#CONDOR request_cpus = 8
#CONDOR request_memory = 16384
#CONDOR +MaxRuntime = 7200
#CONDOR notification = Always
#CONDOR notify_user = roundtrip@example.com

echo "Running job"
`
	scriptPath := filepath.Join(tmpDir, "original.sh")
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
	if !specs.EmailOnBegin || !specs.EmailOnEnd || !specs.EmailOnFail {
		t.Errorf("Parsing failed: EmailOnBegin=%v, EmailOnEnd=%v, EmailOnFail=%v",
			specs.EmailOnBegin, specs.EmailOnEnd, specs.EmailOnFail)
	}
	if specs.MailUser != "roundtrip@example.com" {
		t.Errorf("MailUser = %q; want %q", specs.MailUser, "roundtrip@example.com")
	}
	if specs.Time != 2*time.Hour {
		t.Errorf("Time = %v; want 2h", specs.Time)
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
			RawFlags: []string{},
			Ncpus:    1,
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
				"#!/bin/bash",
				"#CONDOR request_cpus = 8",
				"#CONDOR request_memory = 16384",
				"#CONDOR request_gpus = 2",
				"#CONDOR +MaxRuntime = 3600",
			},
			wantCpus:  8,
			wantMemMB: 16384,
			wantGpus:  2,
			wantTime:  time.Hour,
		},
		{
			name: "Memory with GB suffix",
			scriptLines: []string{
				"#!/bin/bash",
				"#CONDOR request_cpus = 4",
				"#CONDOR request_memory = 8GB",
			},
			wantCpus:  4,
			wantMemMB: 8192,
			wantGpus:  0,
			wantTime:  0,
		},
		{
			name: "Memory with MB suffix",
			scriptLines: []string{
				"#!/bin/bash",
				"#CONDOR request_memory = 4096MB",
			},
			wantCpus:  4, // default
			wantMemMB: 4096,
			wantGpus:  0,
			wantTime:  0,
		},
		{
			name: "Plain memory (default MB)",
			scriptLines: []string{
				"#!/bin/bash",
				"#CONDOR request_memory = 2048",
			},
			wantCpus:  4, // default
			wantMemMB: 2048,
			wantGpus:  0,
			wantTime:  0,
		},
		{
			name: "Only CPUs",
			scriptLines: []string{
				"#!/bin/bash",
				"#CONDOR request_cpus = 16",
			},
			wantCpus:  16,
			wantMemMB: 0,
			wantGpus:  0,
			wantTime:  0,
		},
		{
			name: "Only GPUs",
			scriptLines: []string{
				"#!/bin/bash",
				"#CONDOR request_gpus = 4",
			},
			wantCpus:  4, // default
			wantMemMB: 0,
			wantGpus:  4,
			wantTime:  0,
		},
		{
			name: "Time in seconds",
			scriptLines: []string{
				"#!/bin/bash",
				"#CONDOR +MaxRuntime = 86400",
			},
			wantCpus:  4, // default
			wantMemMB: 0,
			wantGpus:  0,
			wantTime:  24 * time.Hour,
		},
		{
			name: "No resource directives (defaults)",
			scriptLines: []string{
				"#!/bin/bash",
				"echo hello",
			},
			wantCpus:  4, // default
			wantMemMB: 0,
			wantGpus:  0,
			wantTime:  0,
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			tmpDir := t.TempDir()
			scriptPath := filepath.Join(tmpDir, "test.sh")
			content := strings.Join(tt.scriptLines, "\n")
			if err := os.WriteFile(scriptPath, []byte(content), 0644); err != nil {
				t.Fatalf("Failed to create test script: %v", err)
			}

			htcondor := newTestHTCondorScheduler()
			specs, err := htcondor.ReadScriptSpecs(scriptPath)
			if err != nil {
				t.Fatalf("Failed to parse script: %v", err)
			}

			if specs.Ncpus != tt.wantCpus {
				t.Errorf("Ncpus = %d; want %d", specs.Ncpus, tt.wantCpus)
			}
			if specs.MemMB != tt.wantMemMB {
				t.Errorf("MemMB = %d; want %d", specs.MemMB, tt.wantMemMB)
			}
			if tt.wantGpus > 0 {
				if specs.Gpu == nil {
					t.Errorf("Gpu is nil; want count %d", tt.wantGpus)
				} else if specs.Gpu.Count != tt.wantGpus {
					t.Errorf("Gpu.Count = %d; want %d", specs.Gpu.Count, tt.wantGpus)
				}
			} else {
				if specs.Gpu != nil {
					t.Errorf("Gpu = %+v; want nil", specs.Gpu)
				}
			}
			if specs.Time != tt.wantTime {
				t.Errorf("Time = %v; want %v", specs.Time, tt.wantTime)
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
			JobName:  "format_test",
			Ncpus:    8,
			MemMB:    16384,
			Time:     2 * time.Hour,
			RawFlags: []string{
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
			scriptPath := filepath.Join(tmpDir, "test.sh")
			content := fmt.Sprintf("#!/bin/bash\n#CONDOR +MaxRuntime = %s\n", tt.seconds)
			if err := os.WriteFile(scriptPath, []byte(content), 0644); err != nil {
				t.Fatalf("Failed to create test script: %v", err)
			}

			htcondor := newTestHTCondorScheduler()
			specs, err := htcondor.ReadScriptSpecs(scriptPath)
			if err != nil {
				t.Fatalf("Failed to parse script: %v", err)
			}
			if specs.Time != tt.wantTime {
				t.Errorf("Time = %v; want %v", specs.Time, tt.wantTime)
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

	script := `#!/bin/bash
#CONDOR request_cpus = 4
#CONDOR request_memory = 8192
#CONDOR request_gpus = 1
#CONDOR +MaxRuntime = 7200
#CONDOR notification = Complete
#CONDOR notify_user = user@example.com

echo "Running job"
`
	scriptPath := filepath.Join(tmpDir, "test.sh")
	if err := os.WriteFile(scriptPath, []byte(script), 0644); err != nil {
		t.Fatalf("Failed to create test script: %v", err)
	}

	specs, err := TryParseHTCondorScript(scriptPath)
	if err != nil {
		t.Fatalf("TryParseHTCondorScript failed: %v", err)
	}

	if specs.Ncpus != 4 {
		t.Errorf("Ncpus = %d; want 4", specs.Ncpus)
	}
	if specs.MemMB != 8192 {
		t.Errorf("MemMB = %d; want 8192", specs.MemMB)
	}
	if specs.Gpu == nil || specs.Gpu.Count != 1 {
		t.Errorf("Gpu = %+v; want count 1", specs.Gpu)
	}
	if specs.Time != 2*time.Hour {
		t.Errorf("Time = %v; want 2h", specs.Time)
	}
	if !specs.EmailOnEnd {
		t.Error("EmailOnEnd should be true for Complete notification")
	}
	if specs.MailUser != "user@example.com" {
		t.Errorf("MailUser = %q; want %q", specs.MailUser, "user@example.com")
	}
	if len(specs.RawFlags) != 6 {
		t.Errorf("RawFlags count = %d; want 6", len(specs.RawFlags))
	}
}
