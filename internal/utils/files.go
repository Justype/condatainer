package utils

import (
	"fmt"
	"io/fs"
	"os"
	"path/filepath"
	"strings"
)

// Standard default permissions
// File: u=rw, g=rw, o=r
const PermFile os.FileMode = 0664

// Dir:  u=rwx, g=rwx, o=rx (Requires +x to traverse)
const PermDir os.FileMode = 0775

// FixPermissionsDefault is a helper that applies standard permissions (0664/0775).
// It is a shorthand for FixPermissions(path, PermFile, PermDir).
func FixPermissionsDefault(path string) error {
	return FixPermissions(path, PermFile, PermDir)
}

// FixPermissions automatically adjusts permissions for a path using the provided modes.
// If 'path' is a file, it sets it to 'filePerm'.
// If 'path' is a directory, it recursively sets files to 'filePerm' and subdirs to 'dirPerm'.
func FixPermissions(path string, filePerm, dirPerm os.FileMode) error {
	info, err := os.Stat(path)
	if err != nil {
		return fmt.Errorf("could not stat path %s: %w", path, err)
	}

	// 1. Handle Single File
	if !info.IsDir() {
		PrintDebug("Fixing permissions for file: %s [%v]", StylePath(path), filePerm)
		if err := os.Chmod(path, filePerm); err != nil {
			return fmt.Errorf("failed to chmod file %s: %w", path, err)
		}
		return nil
	}

	// 2. Handle Directory (Recursive Walk)
	PrintDebug("Recursively fixing permissions for directory: %s [Dirs:%v, Files:%v]",
		StylePath(path), dirPerm, filePerm)

	// First, fix the root directory itself
	if err := os.Chmod(path, dirPerm); err != nil {
		return fmt.Errorf("failed to chmod root dir %s: %w", path, err)
	}

	// Then walk the contents
	return filepath.WalkDir(path, func(p string, d fs.DirEntry, err error) error {
		if err != nil {
			PrintWarning("Skipping inaccessible path during permission fix: %s", StylePath(p))
			return nil // Skip this file/dir but continue walking
		}

		// Calculate target mode based on the input parameters
		var targetMode os.FileMode
		if d.IsDir() {
			targetMode = dirPerm
		} else {
			targetMode = filePerm
		}

		// Optimization: Only run Chmod if permissions actually differ
		info, err := d.Info()
		if err == nil && info.Mode().Perm() != targetMode {
			if err := os.Chmod(p, targetMode); err != nil {
				PrintWarning("Could not chmod %s: %v", StylePath(p), err)
			}
		}

		return nil
	})
}

// ShareToUGORecursive recursively sets permissions to share files with user, group, and others.
// Files: ug+rw,o+r (0664)
// Directories: ug+rwx,o+rx (0775)
// This matches the Python implementation's share_to_ugo_recursive function.
func ShareToUGORecursive(path string) error {
	info, err := os.Stat(path)
	if err != nil {
		return fmt.Errorf("could not stat path %s: %w", path, err)
	}

	// 1. Handle Single File
	if !info.IsDir() {
		// Files: ug+rw,o+r
		currentMode := info.Mode()
		newMode := currentMode | 0644 // u+rw, o+r
		newMode = newMode | 0060      // g+rw
		if err := os.Chmod(path, newMode); err != nil {
			PrintWarning("Failed to set permissions for %s: %v", StylePath(path), err)
			return err
		}
		return nil
	}

	// 2. Handle Directory (Recursive Walk)
	return filepath.WalkDir(path, func(p string, d fs.DirEntry, err error) error {
		if err != nil {
			PrintWarning("Skipping inaccessible path: %s", StylePath(p))
			return nil // Continue walking
		}

		fileInfo, err := d.Info()
		if err != nil {
			PrintWarning("Could not get info for %s: %v", StylePath(p), err)
			return nil
		}

		currentMode := fileInfo.Mode()
		var newMode os.FileMode

		if d.IsDir() {
			// Directories: ug+rwx,o+rx
			newMode = currentMode | 0755 // u+rwx, o+rx
			newMode = newMode | 0070     // g+rwx
		} else {
			// Files: ug+rw,o+r
			newMode = currentMode | 0644 // u+rw, o+r
			newMode = newMode | 0060     // g+rw
		}

		if err := os.Chmod(p, newMode); err != nil {
			PrintWarning("Could not chmod %s: %v", StylePath(p), err)
		}

		return nil
	})
}

// --- Extension Checks (String-based) ---

// IsImg checks if the path has an ext3 overlay extension (.img).
// Note: In Apptainer context, .img usually implies a writable ext3 overlay.
func IsImg(path string) bool {
	ext := strings.ToLower(filepath.Ext(path))
	return ext == ".img"
}

// IsSqf checks if the path has a SquashFS extension (.sqf, .sqsh, .squashfs).
// These are read-only compressed images.
func IsSqf(path string) bool {
	ext := strings.ToLower(filepath.Ext(path))
	return ext == ".sqf" || ext == ".sqsh" || ext == ".squashfs"
}

// IsSif checks if the path has a Singularity Image Format extension (.sif).
// This is the native format for Apptainer/Singularity.
func IsSif(path string) bool {
	ext := strings.ToLower(filepath.Ext(path))
	return ext == ".sif"
}

// IsApptainerImage checks if the file matches any supported Apptainer image format.
// Returns true for .sif, .img, .sqf, .sqsh, .squashfs.
func IsApptainerImage(path string) bool {
	return IsSif(path) || IsImg(path) || IsSqf(path)
}

// IsOverlay checks if the path is an overlay file (.img, .sqf, .sqsh, .squashfs).
// This is used for CondaTainer overlay detection.
func IsOverlay(path string) bool {
	return IsImg(path) || IsSqf(path)
}

// IsYaml checks if the path has a YAML extension (.yaml, .yml).
// Useful for Conda environment definition files.
func IsYaml(path string) bool {
	ext := strings.ToLower(filepath.Ext(path))
	return ext == ".yaml" || ext == ".yml"
}

// --- Filesystem Checks (OS-based) ---

// FileExists checks if a file exists and is not a directory.
func FileExists(path string) bool {
	info, err := os.Stat(path)
	if os.IsNotExist(err) {
		return false
	}
	return !info.IsDir()
}

// DirExists checks if a path exists and is a directory.
func DirExists(path string) bool {
	info, err := os.Stat(path)
	if os.IsNotExist(err) {
		return false
	}
	return info.IsDir()
}

// EnsureDir checks if a directory exists, and creates it if it doesn't.
func EnsureDir(path string) error {
	if DirExists(path) {
		return nil
	}
	return os.MkdirAll(path, 0775)
}
