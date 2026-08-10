package utils

import (
	"bytes"
	"context"
	"errors"
	"fmt"
	"os"
	"os/exec"
	"path/filepath"
	"strings"
	"syscall"
)

// MoveFile moves src to dst, falling back to a copy across filesystems, and
// reports whether a copy was performed. sparse=true preserves holes.
//
// Both branches share the result with the destination's group: which one ran is
// an accident of how the machine is mounted and must not decide who can write.
func MoveFile(ctx context.Context, src, dst string, sparse bool) (copied bool, err error) {
	if err := MkdirAllShared(filepath.Dir(dst)); err != nil {
		return false, fmt.Errorf("create destination directory: %w", err)
	}

	if err := os.Rename(src, dst); err == nil {
		// A rename carries the source inode's mode across untouched.
		ShareWithParentGroup(dst)
		return false, nil
	} else if !errors.Is(err, syscall.EXDEV) {
		return false, fmt.Errorf("rename %s → %s: %w", src, dst, err)
	}

	// Separate mounts for scratch and images is the normal HPC layout, not an
	// edge case.
	if err := crossFsCopy(ctx, src, dst, sparse); err != nil {
		return false, err
	}
	os.Remove(src) //nolint:errcheck
	ShareWithParentGroup(dst)
	return true, nil
}

// crossFsCopy copies src to dst with the system cp, so a long copy stays
// interruptible rather than leaving a truncated image behind.
func crossFsCopy(ctx context.Context, src, dst string, sparse bool) error {
	sparseFlag := "--sparse=never"
	if sparse {
		sparseFlag = "--sparse=always"
	}
	var outBuf bytes.Buffer
	cmd := exec.Command("cp", sparseFlag, src, dst)
	cmd.Stdout = &outBuf
	cmd.Stderr = &outBuf
	if err := cmd.Start(); err != nil {
		return fmt.Errorf("cp %s → %s: %w", src, dst, err)
	}

	waitErr := make(chan error, 1)
	go func() { waitErr <- cmd.Wait() }()

	select {
	case err := <-waitErr:
		if err != nil {
			os.Remove(dst) //nolint:errcheck
			return fmt.Errorf("cp %s → %s: %w\n%s", src, dst, err, strings.TrimSpace(outBuf.String()))
		}
	case <-ctx.Done():
		cmd.Process.Kill() //nolint:errcheck
		os.Remove(dst)     //nolint:errcheck
		return ctx.Err()
	}
	return nil
}
