package utils

import (
	"os"
	"path/filepath"
)

// EnsureTmpSubdir creates the leaf directory dir with permissions based on the parent:
//   - Parent world-writable (o+w, e.g. /tmp with mode 1777): use 0700 so other users
//     on the same node cannot read or list the contents; not group-shared.
//   - Otherwise (private path, or a group-shared scratch): use PermDir and
//     ShareWithParentGroup, so a group-writable (2775) parent propagates group-write
//     down into the workspace built underneath.
//
// Only the leaf is created; the parent must already exist.
// An already-existing directory is accepted silently (safe for concurrent builds).
func EnsureTmpSubdir(dir string) error {
	if info, err := os.Stat(filepath.Dir(dir)); err == nil && info.Mode()&0002 != 0 {
		// Parent world-writable: keep the scratch private, don't share.
		if err := os.Mkdir(dir, 0700); err != nil && !os.IsExist(err) {
			return err
		}
		return nil
	}
	if err := os.Mkdir(dir, PermDir); err != nil && !os.IsExist(err) {
		return err
	}
	ShareWithParentGroup(dir)
	return nil
}

// GetTmpDir returns the best available tmp directory for condatainer builds.
// Priority:
//  1. CNT_TMPDIR (explicit override)
//  2. Scheduler-specific local node storage: SLURM_TMPDIR, PBS_TMPDIR, LSF_TMPDIR, _CONDOR_SCRATCH_DIR
//  3. TMPDIR (POSIX standard)
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

	// Priority 3: POSIX standard tmp dir
	if v := os.Getenv("TMPDIR"); v != "" {
		return filepath.Join(v, sub)
	}

	return filepath.Join("/tmp", sub)
}
