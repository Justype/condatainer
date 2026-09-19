// Package helperhistory reads and writes cnt-lock/.helper-history/, the
// shared, cross-user record of which overlay combinations a project's
// helpers have been run with. It is a leaf package — internal/helper (the
// fresh-start reuse lookup) and internal/project (the project-wide "used but
// not pinned" index) both depend on it, so neither has to depend on the
// other just to read this data.
package helperhistory

import (
	"encoding/json"
	"fmt"
	"os"
	"path/filepath"
	"slices"
	"sort"
	"strings"
	"time"

	"github.com/Justype/condatainer/internal/project/lock"
	"github.com/Justype/condatainer/internal/utils"
	"golang.org/x/sys/unix"
)

// SubdirName is the gitignored dot directory inside a project's cnt-lock/
// holding one file per distinct overlay combination a helper has
// been run with, shared across every user of the project.
const SubdirName = ".helper-history"

// record is one shared history entry, marshaled as-is: no hash, no
// last-used field. Its filename exists only to be unique on disk — content
// is what a reader compares, and recency comes from the file's own mtime.
type record struct {
	Helper   string   `json:"helper"`
	Location string   `json:"location"`
	Required []string `json:"required,omitempty"`
	Overlays []string `json:"overlays"`
}

// UsedCombination is one shared record as ListUsed/ListAll return it.
// Helper/Location are already the map keys (ListAll) or the arguments
// (ListUsed) by the time a caller sees it, so only what's left is exposed.
type UsedCombination struct {
	Required []string  `json:"required,omitempty"`
	Overlays []string  `json:"overlays"`
	LastUsed time.Time `json:"last_used"`
}

// Dir is root's .helper-history/ directory, inside cnt-lock/.
func Dir(root string) string {
	return filepath.Join(root, lock.DirName, SubdirName)
}

// convertLocationForFilename makes location safe to embed in a filename.
// Purely for uniqueness — a reader never parses it back out.
func convertLocationForFilename(location string) string {
	if location == "" || location == "." {
		return "root"
	}
	return strings.ReplaceAll(location, "/", "-")
}

// RecordUsed records that name was run at location (project-root-relative,
// slash-separated) with the helper's required overlays and the overlays the
// user added, in root's shared history. A no-op when root is "" — no standing
// project means no shared audience to record for.
//
// Reads every file already in .helper-history/ and either bumps the matching
// combination's mtime (a repeat use) or writes a new entry (a fresh one).
func RecordUsed(root, name, location string, required, overlays []string) error {
	if root == "" {
		return nil
	}
	sortedRequired := sortedCopy(required)
	sorted := sortedCopy(overlays)

	dir := Dir(root)
	if err := utils.MkdirAllShared(dir); err != nil {
		return err
	}
	if err := lock.EnsureIgnore(root); err != nil {
		return err
	}
	entries, err := os.ReadDir(dir)
	if err != nil {
		return err
	}
	for _, e := range entries {
		if e.IsDir() || !strings.HasSuffix(e.Name(), ".json") {
			continue
		}
		path := filepath.Join(dir, e.Name())
		data, err := os.ReadFile(path)
		if err != nil {
			continue
		}
		var rec record
		if err := json.Unmarshal(data, &rec); err != nil {
			continue
		}
		if rec.Helper == name && rec.Location == location &&
			slices.Equal(rec.Required, sortedRequired) && slices.Equal(rec.Overlays, sorted) {
			return bumpMtime(path)
		}
	}

	rec := record{Helper: name, Location: location, Required: sortedRequired, Overlays: sorted}
	data, err := json.Marshal(rec)
	if err != nil {
		return err
	}
	filename := fmt.Sprintf("%s-%s-%d.json", name, convertLocationForFilename(location), time.Now().UnixNano())
	f, err := utils.CreateFileWritable(filepath.Join(dir, filename))
	if err != nil {
		return err
	}
	defer f.Close()
	_, err = f.Write(data)
	return err
}

func sortedCopy(in []string) []string {
	out := append([]string(nil), in...)
	sort.Strings(out)
	return out
}

// bumpMtime touches path's mtime to now, without reading or rewriting its
// content. Uses unix.UtimesNanoAt with UTIME_OMIT (atime) and the UTIME_NOW
// sentinel (mtime) rather than os.Chtimes: os.Chtimes always passes an
// explicit timestamp, which utimensat(2) only lets the file's owner set,
// group-write notwithstanding. UTIME_NOW waives that requirement and needs
// only write access — exactly what ShareWithParentGroup already guarantees
// every group member on a shared cnt-lock/.
func bumpMtime(path string) error {
	ts := []unix.Timespec{
		{Sec: 0, Nsec: unix.UTIME_OMIT},
		{Sec: 0, Nsec: unix.UTIME_NOW},
	}
	return unix.UtimesNanoAt(unix.AT_FDCWD, path, ts, 0)
}

// recordFile pairs a parsed record with the mtime of the file it came from —
// the sole source of recency; nothing about "last used" is stored in the
// record itself.
type recordFile struct {
	rec      record
	lastUsed time.Time
}

// readAll reads every entry in root's .helper-history/ once. Both ListUsed
// and ListAll build on this rather than scanning separately — finding one
// location's combinations, or aggregating across all of them, are two
// different questions asked of the same read.
func readAll(root string) ([]recordFile, error) {
	dir := Dir(root)
	entries, err := os.ReadDir(dir)
	if err != nil {
		if os.IsNotExist(err) {
			return nil, nil
		}
		return nil, err
	}
	out := make([]recordFile, 0, len(entries))
	for _, e := range entries {
		if e.IsDir() || !strings.HasSuffix(e.Name(), ".json") {
			continue
		}
		data, err := os.ReadFile(filepath.Join(dir, e.Name()))
		if err != nil {
			continue
		}
		var rec record
		if err := json.Unmarshal(data, &rec); err != nil {
			continue
		}
		info, err := e.Info()
		if err != nil {
			continue
		}
		out = append(out, recordFile{rec: rec, lastUsed: info.ModTime()})
	}
	return out, nil
}

// ListUsed returns every combination recorded for name at location, newest
// first — the fresh-start reuse lookup's one-location view. Consumed the
// moment someone actually starts that helper there; never shown as a list
// to browse.
func ListUsed(root, name, location string) ([]*UsedCombination, error) {
	all, err := readAll(root)
	if err != nil {
		return nil, err
	}
	var out []*UsedCombination
	for _, rf := range all {
		if rf.rec.Helper != name || rf.rec.Location != location {
			continue
		}
		out = append(out, &UsedCombination{Required: rf.rec.Required, Overlays: rf.rec.Overlays, LastUsed: rf.lastUsed})
	}
	sort.Slice(out, func(i, j int) bool { return out[i].LastUsed.After(out[j].LastUsed) })
	return out, nil
}

// ListAll groups every combination in root's shared history by helper, then
// location. Its only consumer is internal/project's usageIndex aggregate
// computation — never a caller that displays the combinations themselves.
func ListAll(root string) (map[string]map[string][]*UsedCombination, error) {
	all, err := readAll(root)
	if err != nil {
		return nil, err
	}
	out := make(map[string]map[string][]*UsedCombination, len(all))
	for _, rf := range all {
		byLocation := out[rf.rec.Helper]
		if byLocation == nil {
			byLocation = make(map[string][]*UsedCombination)
			out[rf.rec.Helper] = byLocation
		}
		byLocation[rf.rec.Location] = append(byLocation[rf.rec.Location],
			&UsedCombination{Required: rf.rec.Required, Overlays: rf.rec.Overlays, LastUsed: rf.lastUsed})
	}
	return out, nil
}
