package server

import (
	"encoding/json"
	"fmt"
	"io"
	"net/http"
	"os"
	"path/filepath"
	"strings"

	"github.com/Justype/condatainer/internal/utils"
)

// resolveUploadTarget validates that relPath is a safe relative path and
// returns the absolute target path joined under destAbs. Rejects empty or
// absolute paths, and any path that escapes destAbs after cleaning (e.g. via
// "..").
func resolveUploadTarget(destAbs, relPath string) (string, error) {
	if relPath == "" {
		return "", fmt.Errorf("empty file name")
	}
	if filepath.IsAbs(relPath) {
		return "", fmt.Errorf("absolute paths are not allowed: %s", relPath)
	}
	cleaned := filepath.Clean(relPath)
	target := filepath.Join(destAbs, cleaned)
	rel, err := filepath.Rel(destAbs, target)
	if err != nil || rel == ".." || strings.HasPrefix(rel, ".."+string(filepath.Separator)) {
		return "", fmt.Errorf("invalid path: %s", relPath)
	}
	return target, nil
}

// resolveUploadDest validates the dest query param: it must resolve to an
// existing directory. Returns its absolute path.
func resolveUploadDest(r *http.Request) (string, error) {
	dest := r.URL.Query().Get("dest")
	if dest == "" {
		return "", fmt.Errorf("dest is required")
	}
	abs, err := filepath.Abs(dest)
	if err != nil {
		return "", fmt.Errorf("invalid dest: %w", err)
	}
	info, err := os.Stat(abs)
	if err != nil {
		return "", fmt.Errorf("dest does not exist: %w", err)
	}
	if !info.IsDir() {
		return "", fmt.Errorf("dest is not a directory")
	}
	return abs, nil
}

// handleFSUploadCheck serves POST /api/fs/upload/check?dest=... — given a
// JSON body {"paths": [...]} of relative paths a client intends to upload,
// reports which of them already exist under dest. Also validates every path
// (rejecting traversal attempts) so a subsequent handleFSUpload call only
// needs to re-check, not re-derive, path safety.
func (s *srv) handleFSUploadCheck(w http.ResponseWriter, r *http.Request) {
	if r.Method != http.MethodPost {
		http.Error(w, "method not allowed", http.StatusMethodNotAllowed)
		return
	}
	destAbs, err := resolveUploadDest(r)
	if err != nil {
		http.Error(w, err.Error(), http.StatusBadRequest)
		return
	}

	var req struct {
		Paths []string `json:"paths"`
	}
	if err := json.NewDecoder(r.Body).Decode(&req); err != nil {
		http.Error(w, "invalid request body", http.StatusBadRequest)
		return
	}

	var conflicts []string
	for _, p := range req.Paths {
		target, err := resolveUploadTarget(destAbs, p)
		if err != nil {
			http.Error(w, err.Error(), http.StatusBadRequest)
			return
		}
		if _, err := os.Stat(target); err == nil {
			conflicts = append(conflicts, p)
		}
	}

	writeJSON(w, map[string]interface{}{"conflicts": conflicts})
}

// handleFSUpload serves POST /api/fs/upload?dest=...&overwrite=true|false —
// accepts a multipart/form-data body whose parts each carry a relative path
// in their form field name (may contain "/" to recreate subdirectories under
// dest). The relative path is carried in the field name rather than the
// filename because Go's mime/multipart deliberately sanitizes Part.FileName()
// down to its basename (a hardening measure against path traversal via a
// client-supplied filename) — Part.FormName() is not sanitized, so it's the
// only place a "/"-containing path survives parsing intact. Streams each
// part directly to disk via io.Copy (no ParseMultipartForm buffering) so
// large files aren't spooled in memory.
//
// A part whose copy fails mid-file (e.g. the client aborted the upload) is
// removed rather than left truncated on disk; fully-written earlier parts
// are kept.
//
// When overwrite=false, parts whose target already exists are skipped
// (recorded in the "skipped" response field) rather than overwritten; the
// dashboard's own upload UI always resolves conflicts client-side first
// (see the sequential Yes/Yes to All/No/No to All prompt in files.js) and
// calls with overwrite=true, only ever including files it decided to write.
func (s *srv) handleFSUpload(w http.ResponseWriter, r *http.Request) {
	if r.Method != http.MethodPost {
		http.Error(w, "method not allowed", http.StatusMethodNotAllowed)
		return
	}
	destAbs, err := resolveUploadDest(r)
	if err != nil {
		http.Error(w, err.Error(), http.StatusBadRequest)
		return
	}
	overwrite := r.URL.Query().Get("overwrite") == "true"

	mr, err := r.MultipartReader()
	if err != nil {
		http.Error(w, "expected multipart/form-data body", http.StatusBadRequest)
		return
	}

	var written, skipped []string
	for {
		part, err := mr.NextPart()
		if err == io.EOF {
			break
		}
		if err != nil {
			http.Error(w, err.Error(), http.StatusBadRequest)
			return
		}

		relPath := part.FormName()
		target, err := resolveUploadTarget(destAbs, relPath)
		if err != nil {
			part.Close()
			http.Error(w, fmt.Sprintf("%s (partial upload: %d files already written)", err.Error(), len(written)), http.StatusBadRequest)
			return
		}

		if !overwrite {
			if _, statErr := os.Stat(target); statErr == nil {
				skipped = append(skipped, relPath)
				io.Copy(io.Discard, part) //nolint:errcheck
				part.Close()
				continue
			}
		}

		if err := utils.MkdirAllShared(filepath.Dir(target)); err != nil {
			part.Close()
			http.Error(w, err.Error(), http.StatusInternalServerError)
			return
		}
		out, err := utils.CreateFileWritable(target)
		if err != nil {
			part.Close()
			http.Error(w, err.Error(), http.StatusInternalServerError)
			return
		}
		_, copyErr := io.Copy(out, part)
		out.Close()
		part.Close()
		if copyErr != nil {
			os.Remove(target) // drop the partial file rather than leaving it truncated
			http.Error(w, copyErr.Error(), http.StatusInternalServerError)
			return
		}
		written = append(written, relPath)
	}

	writeJSON(w, map[string]interface{}{"written": written, "skipped": skipped})
}
