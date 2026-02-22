package overlay

import (
	"os"
	"os/exec"
	"path/filepath"
	"strconv"
	"strings"

	"github.com/Justype/condatainer/internal/utils"
)

// ============================================================================
// Low-level SquashFS operations
// ============================================================================

// PathExists returns true if the given entry exists at the top level of the
// SquashFS archive. unsquashfs always exits 0, but prints only "squashfs-root"
// when the entry is absent. A second line appears when found — so line count > 1 means found.
// The process is killed after the second newline for efficiency.
func PathExists(sqfPath, entry string) bool {
	cmd := exec.Command("unsquashfs", "-l", sqfPath, entry)
	cmd.Env = append(os.Environ(), "LC_ALL=C")
	pipe, err := cmd.StdoutPipe()
	if err != nil {
		return false
	}
	if err := cmd.Start(); err != nil {
		return false
	}
	defer cmd.Process.Kill() //nolint:errcheck

	buf := make([]byte, 4096)
	count := 0
	found := false
	for {
		n, err := pipe.Read(buf)
		for _, b := range buf[:n] {
			if b == '\n' {
				count++
				if count >= 2 {
					found = true
					goto done
				}
			}
		}
		if err != nil {
			break
		}
	}
done:
	return found
}

// Cat reads a file from a SquashFS archive using unsquashfs -cat.
// Returns nil if the file cannot be read.
func Cat(sqfPath, filePath string) []byte {
	cmd := exec.Command("unsquashfs", "-cat", sqfPath, filePath)
	cmd.Env = append(os.Environ(), "LC_ALL=C")
	out, err := cmd.Output()
	if err != nil {
		return nil
	}
	return out
}

// IsOSType returns true if the SquashFS archive is an Apptainer OS image,
// identified by the presence of .singularity.d at the archive root.
func IsOSType(sqfPath string) bool {
	return strings.HasSuffix(sqfPath, ".sqf") && PathExists(sqfPath, ".singularity.d")
}

// HasCntBin returns true if the SquashFS archive contains a bin directory at
// cnt/<name>/<version>/bin. nameVersion accepts both slash form ("samtools/1.22")
// and the normalized double-dash form ("samtools--1.22").
func HasCntBin(sqfPath, nameVersion string) bool {
	nv := utils.NormalizeNameVersion(nameVersion)
	return PathExists(sqfPath, "cnt/"+nv+"/bin")
}

// ============================================================================
// SquashFS stats
// ============================================================================

// SquashFSStats holds metadata parsed from `unsquashfs -stat`.
type SquashFSStats struct {
	CreatedTime         string
	FilesystemSizeBytes int64
	Compression         string
	CompressionLevel    int // 0 means not specified by the file
	BlockSize           int64
	NumFragments        int64
	NumInodes           int64
	NumIDs              int64
	DuplicatesRemoved   bool
	ExportableNFS       bool
}

// GetSquashFSStats runs `unsquashfs -stat` on the given path and parses its output.
func GetSquashFSStats(path string) (*SquashFSStats, error) {
	cmd := exec.Command("unsquashfs", "-stat", path)
	cmd.Env = append(os.Environ(), "LC_ALL=C", "LC_TIME=C")
	out, err := cmd.Output()
	if err != nil {
		return nil, &Error{Op: "stat", Path: path, Tool: "unsquashfs", BaseErr: err}
	}

	stats := &SquashFSStats{}
	lines := strings.Split(string(out), "\n")

	for i, raw := range lines {
		line := strings.TrimSpace(raw)
		switch {
		case strings.HasPrefix(line, "Creation or last append time "):
			stats.CreatedTime = strings.TrimPrefix(line, "Creation or last append time ")

		case strings.HasPrefix(line, "Filesystem size "):
			// "Filesystem size 24511225 bytes (23936.74 Kbytes / 23.38 Mbytes)"
			parts := strings.Fields(line)
			if len(parts) >= 3 {
				stats.FilesystemSizeBytes, _ = strconv.ParseInt(parts[2], 10, 64)
			}

		case strings.HasPrefix(line, "Compression "):
			stats.Compression = strings.TrimPrefix(line, "Compression ")
			// Next line may be indented "compression-level N"
			if i+1 < len(lines) {
				next := strings.TrimSpace(lines[i+1])
				if strings.HasPrefix(next, "compression-level ") {
					lvl, _ := strconv.Atoi(strings.TrimPrefix(next, "compression-level "))
					stats.CompressionLevel = lvl
				}
			}

		case strings.HasPrefix(line, "Block size "):
			stats.BlockSize, _ = strconv.ParseInt(strings.TrimPrefix(line, "Block size "), 10, 64)

		case strings.HasPrefix(line, "Number of fragments "):
			stats.NumFragments, _ = strconv.ParseInt(strings.TrimPrefix(line, "Number of fragments "), 10, 64)

		case strings.HasPrefix(line, "Number of inodes "):
			stats.NumInodes, _ = strconv.ParseInt(strings.TrimPrefix(line, "Number of inodes "), 10, 64)

		case strings.HasPrefix(line, "Number of ids "):
			stats.NumIDs, _ = strconv.ParseInt(strings.TrimPrefix(line, "Number of ids "), 10, 64)

		case strings.Contains(line, "Duplicates are removed"):
			stats.DuplicatesRemoved = true

		case strings.HasPrefix(line, "Filesystem is exportable via NFS"):
			stats.ExportableNFS = true
		}
	}

	return stats, nil
}

// GetOverlayType classifies a SquashFS overlay based on its contents and filename:
//   - "OS Overlay"     : contains .singularity.d (Apptainer-built)
//   - "Module Overlay" : contains cnt/ and filename uses -- name/version separator
//   - "Bundle Overlay" : contains cnt/ but filename has no --
//   - ""               : unrecognized
func GetOverlayType(sqfPath string) string {
	if PathExists(sqfPath, ".singularity.d") {
		return "OS Overlay"
	}
	if PathExists(sqfPath, "cnt") {
		base := strings.TrimSuffix(filepath.Base(sqfPath), filepath.Ext(sqfPath))
		if strings.Contains(base, "--") {
			return "Module Overlay"
		}
		return "Bundle Overlay"
	}
	return ""
}

// ============================================================================
// OS release info
// ============================================================================

// OSInfo holds distribution identity parsed from /etc/os-release.
type OSInfo struct {
	ID        string   // e.g. "ubuntu", "debian", "rhel", "rocky"
	IDLike    []string // e.g. ["debian"] for ubuntu, ["rhel", "fedora"] for rocky/alma
	VersionID string   // e.g. "24.04", "12", "9.3"
	Codename  string   // e.g. "Noble", "Bookworm", "Blue Onyx" (title-cased)
	Name      string   // e.g. "Ubuntu", "Debian", "Rocky Linux"
}

// GetOSInfo reads /etc/os-release from a SquashFS archive and returns an OSInfo.
// Returns nil if the file cannot be read or parsed.
func GetOSInfo(sqfPath string) *OSInfo {
	data := Cat(sqfPath, "etc/os-release")
	if data == nil {
		return nil
	}
	m := parseOSRelease(string(data))
	if len(m) == 0 {
		return nil
	}
	return newOSInfo(m)
}

// HostOSInfo reads /etc/os-release from the running host.
// Returns nil if the file cannot be read.
func HostOSInfo() *OSInfo {
	data, err := os.ReadFile("/etc/os-release")
	if err != nil {
		return nil
	}
	m := parseOSRelease(string(data))
	if len(m) == 0 {
		return nil
	}
	return newOSInfo(m)
}

// newOSInfo constructs an OSInfo from a parsed os-release key-value map.
func newOSInfo(m map[string]string) *OSInfo {
	o := &OSInfo{
		ID:        m["ID"],
		VersionID: m["VERSION_ID"],
		Name:      strings.ReplaceAll(m["NAME"], " GNU/Linux", ""), // "Debian GNU/Linux" → "Debian"
	}

	if idLike := m["ID_LIKE"]; idLike != "" {
		o.IDLike = strings.Fields(idLike)
	}

	// Prefer VERSION_CODENAME (Debian family); fall back to extracting from VERSION
	codename := m["VERSION_CODENAME"]
	if codename == "" {
		if v := m["VERSION"]; v != "" {
			if i := strings.Index(v, "("); i >= 0 {
				if j := strings.Index(v, ")"); j > i {
					codename = v[i+1 : j]
				}
			}
		}
	}

	// Title-case each word of the codename
	if codename != "" {
		words := strings.Fields(codename)
		for i, w := range words {
			if w != "" {
				words[i] = strings.ToUpper(w[:1]) + w[1:]
			}
		}
		o.Codename = strings.Join(words, " ")
	}

	return o
}

// parseOSRelease parses the KEY=VALUE content of an os-release file.
// Strips surrounding quotes from values; ignores comment and blank lines.
func parseOSRelease(content string) map[string]string {
	result := map[string]string{}
	for line := range strings.SplitSeq(content, "\n") {
		line = strings.TrimSpace(line)
		if line == "" || strings.HasPrefix(line, "#") {
			continue
		}
		k, v, ok := strings.Cut(line, "=")
		if !ok {
			continue
		}
		result[k] = strings.Trim(v, `"'`)
	}
	return result
}

// String returns a human-readable distro string suitable for CLI display.
//
//	"Ubuntu 24.04 (Noble)"
//	"Debian 12 (Bookworm)"
//	"Rocky Linux 9.3 (Blue Onyx)"
//	"Red Hat Enterprise Linux 9.3 (Plow)"
func (o *OSInfo) String() string {
	if o == nil {
		return ""
	}
	result := o.Name
	if o.VersionID != "" {
		result += " " + o.VersionID
	}
	if o.Codename != "" {
		result += " (" + o.Codename + ")"
	}
	return result
}

// MajorVersion returns the major component of VersionID (e.g., "24" from "24.04", "9" from "9.3").
func (o *OSInfo) MajorVersion() string {
	if o == nil || o.VersionID == "" {
		return ""
	}
	if i := strings.Index(o.VersionID, "."); i >= 0 {
		return o.VersionID[:i]
	}
	return o.VersionID
}

// SameFamily returns true if both OSInfo values share an ID or ID_LIKE entry,
// indicating compatibility within the same OS family.
//
//	Ubuntu + Debian  → true  (both ID_LIKE=debian)
//	Rocky  + RHEL    → true  (Rocky has ID_LIKE=rhel)
//	Ubuntu + Rocky   → false
func (o *OSInfo) SameFamily(other *OSInfo) bool {
	if o == nil || other == nil {
		return false
	}
	// Build a set of all IDs that o belongs to
	ids := map[string]bool{}
	if o.ID != "" {
		ids[o.ID] = true
	}
	for _, l := range o.IDLike {
		ids[l] = true
	}
	// Check if other's ID or any IDLike is in the set
	if other.ID != "" && ids[other.ID] {
		return true
	}
	for _, l := range other.IDLike {
		if ids[l] {
			return true
		}
	}
	return false
}

// SameMajorVersion returns true if both OSInfo values share the same major version number.
func (o *OSInfo) SameMajorVersion(other *OSInfo) bool {
	if o == nil || other == nil {
		return false
	}
	mv := o.MajorVersion()
	return mv != "" && mv == other.MajorVersion()
}
