package squashfs

import (
	"bytes"
	"fmt"
	"os"
	"os/exec"
	"path/filepath"
	"strconv"
	"strings"

	"github.com/Justype/condatainer/internal/image/tool"
)

// CatFile reads a file out of a SquashFS archive, distinguishing a missing file
// from a missing tool, an unreadable file and a corrupt archive, so a caller can
// tell "predates the manifest format" from "unsquashfs is missing".
//
// offset skips that many bytes, which is how a SquashFS partition inside a SIF
// is read in place; pass 0 for a plain .sqf.
func CatFile(sqfPath, filePath string, offset int64) ([]byte, error) {
	if !unsquashfsAvailable() {
		return nil, fmt.Errorf("%w: unsquashfs", tool.ErrToolMissing)
	}
	if _, err := os.Stat(sqfPath); err != nil {
		return nil, fmt.Errorf("%w: %s: %w", tool.ErrUnreadable, sqfPath, err)
	}
	inner := strings.TrimPrefix(filePath, "/")

	args := offsetArgs(offset)
	cmd := exec.Command("unsquashfs", append(append(args, "-cat", sqfPath), inner)...)
	cmd.Env = append(os.Environ(), "LC_ALL=C")
	var stderr bytes.Buffer
	cmd.Stderr = &stderr
	out, err := cmd.Output()
	if err == nil {
		return out, nil
	}
	if cause := classifyStderr(stderr.String(), sqfPath); cause != nil {
		return nil, cause
	}

	// Either the file is absent, or this unsquashfs predates -cat (4.5). Only
	// extraction tells the two apart, and it is the fallback either way.
	return catExtractFile(sqfPath, inner, offset)
}

// offsetArgs returns the unsquashfs flags that skip a leading byte offset.
func offsetArgs(offset int64) []string {
	if offset <= 0 {
		return nil
	}
	return []string{"-offset", strconv.FormatInt(offset, 10)}
}

// classifyStderr maps unsquashfs diagnostics to a sentinel, or returns nil when
// the failure is not conclusive and the extraction fallback should decide.
func classifyStderr(stderr, sqfPath string) error {
	low := strings.ToLower(stderr)
	switch {
	case strings.Contains(low, "squashfs superblock"),
		strings.Contains(low, "can't find a squashfs"),
		strings.Contains(low, "not a squashfs"):
		return fmt.Errorf("%w: %s: %s", tool.ErrCorrupt, sqfPath, strings.TrimSpace(stderr))
	case strings.Contains(low, "permission denied"):
		return fmt.Errorf("%w: %s: %s", tool.ErrUnreadable, sqfPath, strings.TrimSpace(stderr))
	case strings.Contains(low, "-offset") && strings.Contains(low, "invalid"):
		return fmt.Errorf("%w: unsquashfs has no -offset support (needs squashfs-tools 4.4+)", tool.ErrToolMissing)
	}
	return nil
}

// catExtractFile is catExtract with the cause preserved. In-archive symlinks are
// resolved manually, since single-file extraction yields a dangling link while
// -cat follows links itself.
func catExtractFile(sqfPath, filePath string, offset int64) ([]byte, error) {
	tmpDir, err := os.MkdirTemp("", "cnt-sqf-cat-")
	if err != nil {
		return nil, fmt.Errorf("%w: %w", tool.ErrUnreadable, err)
	}
	defer os.RemoveAll(tmpDir) //nolint:errcheck

	base := offsetArgs(offset)
	for hop := 0; hop < 4; hop++ {
		dest := filepath.Join(tmpDir, strconv.Itoa(hop))
		args := append(append([]string{}, base...), "-q", "-n", "-d", dest, sqfPath, filePath)
		cmd := exec.Command("unsquashfs", args...)
		cmd.Env = append(os.Environ(), "LC_ALL=C")
		var stderr bytes.Buffer
		cmd.Stderr = &stderr
		if err := cmd.Run(); err != nil {
			if cause := classifyStderr(stderr.String(), sqfPath); cause != nil {
				return nil, cause
			}
			return nil, fmt.Errorf("%w: %s in %s", tool.ErrFileNotFound, filePath, sqfPath)
		}
		extracted := filepath.Join(dest, filePath)
		fi, err := os.Lstat(extracted)
		if err != nil {
			// unsquashfs exits 0 when nothing matched, so an absent output path
			// is how a missing entry actually shows up.
			return nil, fmt.Errorf("%w: %s in %s", tool.ErrFileNotFound, filePath, sqfPath)
		}
		if fi.Mode()&os.ModeSymlink == 0 {
			data, err := os.ReadFile(extracted)
			if err != nil {
				return nil, fmt.Errorf("%w: %w", tool.ErrUnreadable, err)
			}
			return data, nil
		}
		target, err := os.Readlink(extracted)
		if err != nil {
			return nil, fmt.Errorf("%w: %w", tool.ErrUnreadable, err)
		}
		if filepath.IsAbs(target) {
			filePath = strings.TrimPrefix(filepath.Clean(target), "/")
		} else {
			filePath = filepath.Clean(filepath.Join(filepath.Dir(filePath), target))
		}
	}
	return nil, fmt.Errorf("%w: %s in %s: too many symlink hops", tool.ErrFileNotFound, filePath, sqfPath)
}
