package overlay

import (
	"context"
	"fmt"
	"os"
	"path/filepath"

	"github.com/Justype/condatainer/internal/logging"
)

// Resize adjusts the size of an existing ext3 overlay image.
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

	currentSizeBytes := info.Size()
	newSizeBytes := int64(newSizeMB) * 1024 * 1024

	log := logging.FromContext(ctx)

	if currentSizeBytes == newSizeBytes {
		log.Info(fmt.Sprintf("size unchanged (%d MiB) for %s", newSizeMB, absPath))
		return nil
	}

	log.Info(fmt.Sprintf("resizing %s to %d MiB", absPath, newSizeMB))

	if err := CheckIntegrity(ctx, absPath, true); err != nil {
		return err
	}

	if newSizeBytes < currentSizeBytes {
		log.Info(fmt.Sprintf("shrinking overlay image to %d MiB", newSizeMB))
		sizeArg := fmt.Sprintf("%dM", newSizeMB)
		if err := runCommand(ctx, "shrink filesystem", absPath, "resize2fs", "-p", absPath, sizeArg); err != nil {
			return err
		}
		if err := os.Truncate(absPath, newSizeBytes); err != nil {
			return fmt.Errorf("failed to truncate file: %w", err)
		}
	} else {
		log.Info(fmt.Sprintf("expanding overlay image to %d MiB", newSizeMB))
		if err := os.Truncate(absPath, newSizeBytes); err != nil {
			return fmt.Errorf("failed to expand file: %w", err)
		}
		if err := runCommand(ctx, "expand filesystem", absPath, "resize2fs", "-p", absPath); err != nil {
			return err
		}
	}

	if err := CheckIntegrity(ctx, absPath, true); err != nil {
		return err
	}

	log.Info(fmt.Sprintf("overlay image resized to %d MiB: %s", newSizeMB, absPath), "kind", "success")
	return nil
}
