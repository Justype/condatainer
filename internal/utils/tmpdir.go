package utils

import (
	"os"
	"path/filepath"
)

// GetTmpDir returns the best available tmp directory for condatainer builds.
// Priority:
//  1. CNT_TMPDIR (explicit override)
//  2. Scheduler-specific local node storage: SLURM_TMPDIR, PBS_TMPDIR, LSF_TMPDIR, _CONDOR_SCRATCH_DIR
//  3. Common tmp env vars: TMPDIR, TEMP, TMP
//  4. /tmp (always available)
//
// Always appends cnt-$USER to avoid collisions between users.
func GetTmpDir() string {
	user := os.Getenv("USER")
	if user == "" {
		user = "condatainer"
	}
	sub := "cnt-" + user

	// Priority 1: explicit override
	if v := os.Getenv("CNT_TMPDIR"); v != "" {
		return filepath.Join(v, sub)
	}

	// Priority 2: scheduler-assigned local node storage (fast per-job scratch)
	for _, env := range []string{"SLURM_TMPDIR", "PBS_TMPDIR", "LSF_TMPDIR", "_CONDOR_SCRATCH_DIR"} {
		if v := os.Getenv(env); v != "" {
			return filepath.Join(v, sub)
		}
	}

	// Priority 3: common tmp env vars
	for _, env := range []string{"TMPDIR", "TEMP", "TMP"} {
		if v := os.Getenv(env); v != "" {
			return filepath.Join(v, sub)
		}
	}

	return filepath.Join("/tmp", sub)
}
