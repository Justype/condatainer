package overlay

import (
	"fmt"
	"os"
	"path/filepath"

	"github.com/Justype/condatainer/internal/utils"
)

// Resize adjusts the size of an existing ext3 overlay image.
func Resize(imagePath string, newSizeMB int) error {
	// 0. Check Dependencies (resize2fs needs to be checked here, e2fsck is checked in CheckIntegrity)
	if err := checkDependencies([]string{"resize2fs"}); err != nil {
		return err
	}

	absPath, err := filepath.Abs(imagePath)
	if err != nil {
		return fmt.Errorf("failed to resolve path %s: %w", imagePath, err)
	}

	// 1. Validate Target
	info, err := os.Stat(absPath)
	if os.IsNotExist(err) {
		return fmt.Errorf("overlay image not found: %s", absPath)
	}

	currentSizeBytes := info.Size()
	newSizeBytes := int64(newSizeMB) * 1024 * 1024

	if currentSizeBytes == newSizeBytes {
		utils.PrintDebug("Size unchanged (%d MiB) for %s", newSizeMB, utils.StylePath(absPath))
		return nil
	}

	utils.PrintDebug("Resizing %s to %d MiB",
		utils.StylePath(absPath), newSizeMB)

	// 2. Pre-check Integrity (Force = true is recommended before resizing)
	if err := CheckIntegrity(absPath, true); err != nil {
		return err
	}

	// 3. Perform Resize
	if newSizeBytes < currentSizeBytes {
		// --- SHRINK ---
		utils.PrintDebug("Shrinking overlay image to %d MiB", newSizeMB)

		sizeArg := fmt.Sprintf("%dM", newSizeMB)
		if err := runCommand("shrink filesystem", absPath, "resize2fs", "-p", absPath, sizeArg); err != nil {
			return err
		}

		if err := os.Truncate(absPath, newSizeBytes); err != nil {
			return fmt.Errorf("failed to truncate file: %w", err)
		}

	} else {
		// --- EXPAND ---
		utils.PrintDebug("Expanding overlay image to %d MiB", newSizeMB)

		if err := os.Truncate(absPath, newSizeBytes); err != nil {
			return fmt.Errorf("failed to expand file: %w", err)
		}

		if err := runCommand("expand filesystem", absPath, "resize2fs", "-p", absPath); err != nil {
			return err
		}
	}

	// 4. Final Verification
	// We run it again to ensure the resize didn't corrupt anything.
	if err := CheckIntegrity(absPath, true); err != nil {
		return err
	}

	utils.PrintDebug("Overlay image resized successfully to %d MiB: %s", newSizeMB, utils.StylePath(absPath))
	return nil
}
