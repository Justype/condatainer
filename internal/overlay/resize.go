package overlay

import (
	"context"
	"fmt"
	"os"
	"path/filepath"

	"github.com/Justype/condatainer/internal/logging"
	"github.com/Justype/condatainer/internal/utils"
)

// Resize adjusts the size of an existing ext3 overlay image to newSizeMB.
//
// Flow: fsck → grow the file (if growing) → resize2fs to the exact target →
// shrink the file (if shrinking) → fsck. Sizes are compared against the ext3
// filesystem, not the container file's byte size, since the two can diverge.
func Resize(ctx context.Context, imagePath string, newSizeMB int) error {
	if err := checkDependencies([]string{"resize2fs"}); err != nil {
		return err
	}

	absPath, err := filepath.Abs(imagePath)
	if err != nil {
		return fmt.Errorf("failed to resolve path %s: %w", imagePath, err)
	}

	info, err := os.Stat(absPath)
	if os.IsNotExist(err) {
		return fmt.Errorf("overlay image not found: %s", absPath)
	}
	if err != nil {
		return fmt.Errorf("failed to stat %s: %w", absPath, err)
	}

	currentFileBytes := info.Size()
	newSizeBytes := int64(newSizeMB) * 1024 * 1024

	// Current filesystem size from the superblock (block count × block size).
	stats, err := GetStats(absPath)
	if err != nil {
		return err
	}
	fsBytes := stats.TotalBlocks * stats.BlockSize

	log := logging.FromContext(ctx)

	name := filepath.Base(absPath)

	// No-op only when both filesystem and file already match the target.
	if fsBytes == newSizeBytes && currentFileBytes == newSizeBytes {
		log.Info(fmt.Sprintf("size unchanged (%s) for %s",
			utils.StyleNumber(fmt.Sprintf("%d MiB", newSizeMB)), utils.StylePath(name)))
		return nil
	}

	// Reject a shrink below current usage up front — cheap, works on a dirty fs
	// (no fsck). resize2fs still enforces the exact minimum later.
	if newSizeBytes < fsBytes {
		if usedBytes, _ := stats.Usage(); newSizeBytes < usedBytes {
			return fmt.Errorf("cannot shrink %s to %s: filesystem already uses %s",
				utils.StylePath(name),
				utils.StyleNumber(fmt.Sprintf("%d MiB", newSizeMB)),
				utils.StyleNumber(fmt.Sprintf("%d MiB", (usedBytes+1024*1024-1)/(1024*1024))))
		}
	}

	log.Info(fmt.Sprintf("resizing %s to %s (filesystem currently %s)",
		utils.StylePath(name),
		utils.StyleNumber(fmt.Sprintf("%d MiB", newSizeMB)),
		utils.StyleNumber(fmt.Sprintf("%d MiB", fsBytes/(1024*1024)))))

	// resize2fs refuses to run unless the filesystem was force-checked first.
	if err := CheckIntegrity(ctx, absPath, true); err != nil {
		return err
	}

	// Grow the file first so resize2fs has room to expand into.
	if newSizeBytes > currentFileBytes {
		if err := os.Truncate(absPath, newSizeBytes); err != nil {
			return fmt.Errorf("failed to expand file: %w", err)
		}
	}

	// Resize the filesystem to the exact target (handles grow and shrink).
	action := "expand filesystem"
	if newSizeBytes < fsBytes {
		action = "shrink filesystem"
	}
	log.Info(fmt.Sprintf("%s to %s", action, utils.StyleNumber(fmt.Sprintf("%d MiB", newSizeMB))))
	sizeArg := fmt.Sprintf("%dM", newSizeMB)
	if err := runCommand(ctx, action, absPath, "resize2fs", "-p", absPath, sizeArg); err != nil {
		return err
	}

	// Shrink the file down to reclaim the freed space.
	if newSizeBytes < currentFileBytes {
		if err := os.Truncate(absPath, newSizeBytes); err != nil {
			return fmt.Errorf("failed to truncate file: %w", err)
		}
	}

	if err := CheckIntegrity(ctx, absPath, true); err != nil {
		return err
	}

	log.Info(fmt.Sprintf("overlay image resized to %s: %s",
		utils.StyleNumber(fmt.Sprintf("%d MiB", newSizeMB)), utils.StylePath(name)), "kind", "success")
	return nil
}
