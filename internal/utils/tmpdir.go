package utils

import (
	"os"
	"path/filepath"
)

// EnsureTmpSubdir creates the leaf directory dir with permissions based on the parent:
//   - Parent world-writable (o+w, e.g. /tmp with mode 1777): use 0700 so other users
//     on the same node cannot read or list the contents.
//   - Parent private (scheduler job dir, user-controlled path): use PermDir (0775).
//
// Only the leaf is created; the parent must already exist.
// An already-existing directory is accepted silently (safe for concurrent builds).
func EnsureTmpSubdir(dir string) error {
	perm := os.FileMode(PermDir)
	if info, err := os.Stat(filepath.Dir(dir)); err == nil {
		if info.Mode()&0002 != 0 { // parent has o+w → shared/world-writable
			perm = 0700
		}
	}
	if err := os.Mkdir(dir, perm); err != nil && !os.IsExist(err) {
		return err
	}
	return nil
}

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
