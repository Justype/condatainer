package utils

import (
	"compress/gzip"
	"encoding/json"
	"os"
	"path/filepath"
	"strings"

	"golang.org/x/sys/unix"
)

// Standard default permissions
// File: u=rw, g=rw, o=r
const PermFile os.FileMode = 0664

// Dir:  u=rwx, g=rwx, o=rx (Requires +x to traverse)
const PermDir os.FileMode = 0775

// Exec: u=rwx, g=rwx, o=rx (Executable files with group write access)
const PermExec os.FileMode = 0775

// BuildScriptName is the file a script/ref build embeds under its payload dir
// (/cnt/<name>/<version>/) to record the recipe that produced the overlay.
const BuildScriptName = ".cnt-build-script"

// BuildScriptDefName is the file a def or scheme:// (docker://) build embeds at
// the overlay root to record the definition that produced the overlay.
const BuildScriptDefName = ".cnt-build-script.def"

// --- Extension Checks (String-based) ---

// IsImg checks if the path has an ext3 overlay extension (.img, .ext3).
// Note: In Apptainer context, these imply a writable ext3 overlay.
func IsImg(path string) bool {
	ext := strings.ToLower(filepath.Ext(path))
	return ext == ".img" || ext == ".ext3"
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

// IsOverlay checks if the path is an overlay file (.img, .sqf, .sqsh, .squashfs).
// This is used for CondaTainer overlay detection.
func IsOverlay(path string) bool {
	return IsImg(path) || IsSqf(path)
}

// IsCondaFile reports whether path is a Conda environment input file:
// a YAML environment file (.yaml/.yml) or an explicit/spec text file (.txt,
// e.g. the @EXPLICIT output of `export -e`). All are accepted by
// `micromamba create -f`.
func IsCondaFile(path string) bool {
	ext := strings.ToLower(filepath.Ext(path))
	return ext == ".yaml" || ext == ".yml" || ext == ".txt"
}

// --- Filesystem Checks (OS-based) ---

// FindEnvOverlay returns the first env.img found under wd using the standard
// search order. At each location the user-specific variant (env-$USER.img) is
// checked before the generic one:
//
//	{wd}/env-$USER.img → {wd}/env.img
//	{wd}/overlay/env-$USER.img → {wd}/overlay/env.img
//	{wd}/src/overlay/env-$USER.img → {wd}/src/overlay/env.img
//
// envImg may be an explicit path (returned as-is) or "" / "env.img" (triggers search).
// Returns "" if nothing is found.
func FindEnvOverlay(envImg, wd string) string {
	if envImg != "" && envImg != "env.img" {
		return envImg
	}
	base := wd
	if base == "" {
		base, _ = os.Getwd()
	}
	userSuffix := os.Getenv("USER")
	dirs := []string{base, filepath.Join(base, "overlay"), filepath.Join(base, "src", "overlay")}
	for _, dir := range dirs {
		if userSuffix != "" {
			if p := filepath.Join(dir, "env-"+userSuffix+".img"); FileExists(p) {
				return p
			}
		}
		if p := filepath.Join(dir, "env.img"); FileExists(p) {
			return p
		}
	}
	return ""
}

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

// RemoveDirIfEmpty removes dir if it exists and contains no files or subdirectories.
// Silently does nothing if dir is not empty or does not exist.
func RemoveDirIfEmpty(dir string) {
	entries, err := os.ReadDir(dir)
	if err != nil || len(entries) > 0 {
		return
	}
	os.Remove(dir)
}

// CreateFileWritable creates or truncates a file using standard writable file permissions.
func CreateFileWritable(path string) (*os.File, error) {
	f, err := os.OpenFile(path, os.O_CREATE|os.O_WRONLY|os.O_TRUNC, PermFile)
	if err == nil {
		ShareWithParentGroup(path)
	}
	return f, err
}

// ShareWithParentGroup grants the group g+rw (plus g+x on dirs and executables) when
// path's parent is group-writable, and does nothing otherwise — so files land
// group-writable in a shared "2775" install and untouched in a personal one.
// Best-effort: errors are ignored (path may be on another owner's dir).
func ShareWithParentGroup(path string) {
	pi, err := os.Stat(filepath.Dir(path))
	if err != nil || pi.Mode()&0020 == 0 { // parent not group-writable
		return
	}
	fi, err := os.Stat(path)
	if err != nil {
		return
	}
	mode := fi.Mode().Perm()
	if fi.IsDir() {
		mode |= 0070 // g+rwx
	} else {
		mode |= 0060 // g+rw
		if mode&0100 != 0 {
			mode |= 0010 // g+x when u+x (executables)
		}
	}
	_ = os.Chmod(path, mode)
}

// MakeExecutable adds an execute bit wherever the matching read bit is set (x where r),
// then shares with the parent group. Mirroring the read bits keeps it umask-respecting,
// unlike os.Chmod(path, PermExec) which forces 0775 and leaks group/other-write into
// personal installs — so prefer this whenever a created file must be made executable.
func MakeExecutable(path string) error {
	fi, err := os.Stat(path)
	if err != nil {
		return err
	}
	mode := fi.Mode().Perm()
	if mode&0400 != 0 {
		mode |= 0100 // u+x when u+r
	}
	if mode&0040 != 0 {
		mode |= 0010 // g+x when g+r
	}
	if mode&0004 != 0 {
		mode |= 0001 // o+x when o+r
	}
	if err := os.Chmod(path, mode); err != nil {
		return err
	}
	ShareWithParentGroup(path)
	return nil
}

// MkdirAllShared is os.MkdirAll plus ShareWithParentGroup on every level it creates, so
// nested dirs under a shared "2775" tree all become group-writable. Pre-existing dirs
// are left untouched. Prefer this over os.MkdirAll for data/state directories.
func MkdirAllShared(dir string) error {
	// Collect the missing levels first, so we only share the ones we create.
	var created []string // deepest first
	for p := filepath.Clean(dir); ; {
		if DirExists(p) {
			break
		}
		created = append(created, p)
		parent := filepath.Dir(p)
		if parent == p { // reached filesystem root
			break
		}
		p = parent
	}
	if err := os.MkdirAll(dir, PermDir); err != nil {
		return err
	}
	// Shallowest-first, so each level's parent is already shared when we reach it.
	for i := len(created) - 1; i >= 0; i-- {
		ShareWithParentGroup(created[i])
	}
	return nil
}

// EnsureWritableDir creates dir if it does not exist, then checks it is writable.
// Permissions are set only when the directory is newly created.
// Returns true if the directory exists (or was created) and is writable.
func EnsureWritableDir(dir string) bool {
	if !DirExists(dir) {
		if err := os.MkdirAll(dir, PermDir); err != nil {
			return false
		}
		ShareWithParentGroup(dir)
	}
	return CanWriteToDir(dir)
}

// CanWriteToDir checks whether an existing directory is writable using a permission
// check syscall, without creating the directory or writing any files.
// Returns false if the directory does not exist.
// Use this for probe operations (display, bind decisions, search-path scanning).
// Use EnsureWritableDir when you intend to actually create the directory on first use.
func CanWriteToDir(dir string) bool {
	if _, err := os.Stat(dir); err != nil {
		return false
	}
	return unix.Access(dir, unix.W_OK|unix.X_OK) == nil
}

// CanWriteToFile reports whether path can be written to. An existing file must
// itself be writable — a writable parent directory is not enough, since a
// read-only file (e.g. a frozen shared config) cannot be opened for writing.
// A file that does not exist yet only needs a writable ancestor directory.
func CanWriteToFile(path string) bool {
	if _, err := os.Stat(path); err == nil {
		return unix.Access(path, unix.W_OK) == nil
	}
	return CanWriteToExistingAncestor(filepath.Dir(path))
}

// CanWriteToExistingAncestor walks up the path until it finds an existing
// directory and checks whether it is writable. Use this when dir may not exist
// yet but will be created by MkdirAll — a non-existent dir is not read-only,
// it just needs a writable parent.
func CanWriteToExistingAncestor(dir string) bool {
	for d := dir; d != filepath.Dir(d); d = filepath.Dir(d) {
		if _, err := os.Stat(d); err == nil {
			return unix.Access(d, unix.W_OK|unix.X_OK) == nil
		}
	}
	return false
}

// --- Gzip JSON Helpers ---

// ReadGzipJSONFile opens a .json.gz file and decodes JSON into out.
func ReadGzipJSONFile(path string, out any) error {
	f, err := os.Open(path)
	if err != nil {
		return err
	}
	defer f.Close()

	gzReader, err := gzip.NewReader(f)
	if err != nil {
		return err
	}
	defer gzReader.Close()

	return json.NewDecoder(gzReader).Decode(out)
}

// WriteGzipJSONFileAtomic encodes value as JSON, writes it as .json.gz to a temp file,
// then atomically renames it to path.
func WriteGzipJSONFileAtomic(path string, value any) error {
	tmp := path + ".tmp"
	f, err := os.OpenFile(tmp, os.O_CREATE|os.O_WRONLY|os.O_TRUNC, PermFile)
	if err != nil {
		return err
	}
	defer os.Remove(tmp) // no-op after successful rename; cleans up on any error path

	gzWriter, err := gzip.NewWriterLevel(f, gzip.BestSpeed)
	if err != nil {
		f.Close()
		return err
	}
	if err := json.NewEncoder(gzWriter).Encode(value); err != nil {
		gzWriter.Close()
		f.Close()
		return err
	}
	if err := gzWriter.Close(); err != nil {
		f.Close()
		return err
	}
	if err := f.Close(); err != nil {
		return err
	}

	if err := os.Rename(tmp, path); err != nil {
		return err
	}
	ShareWithParentGroup(path)
	return nil
}
