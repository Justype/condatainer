package server

import (
	"encoding/json"
	"net/http"
	"os"
	"path/filepath"
	"sort"
	"strings"

	"github.com/Justype/condatainer/catalog"
	"github.com/Justype/condatainer/internal/artifact/capsule"
	"github.com/Justype/condatainer/internal/config"
	"github.com/Justype/condatainer/internal/image"
	"github.com/Justype/condatainer/internal/project"
	"github.com/Justype/condatainer/internal/project/lock"
	"github.com/Justype/condatainer/internal/project/orchestrate"
	"github.com/Justype/condatainer/internal/utils"
)

// requiredResolution is how one required-overlay name resolves for a launch.
// Resolved is empty when nothing answers it; Ident is set only when the
// project pins it.
type requiredResolution struct {
	Name     string `json:"name"`
	Resolved string `json:"resolved,omitempty"`
	Ident    string `json:"ident,omitempty"`
}

type projectInfo struct {
	Root     string               `json:"root"`
	CWD      string               `json:"cwd"`
	Distro   string               `json:"distro"`
	Required []requiredResolution `json:"required"`
}

// handleProject serves GET /api/project?cwd=&name=... (root, distro and how
// each required-overlay name resolves) and POST /api/project (create cnt-lock/
// at cwd). Reading resolves only what is already on disk.
func (s *srv) handleProject(w http.ResponseWriter, r *http.Request) {
	switch r.Method {
	case http.MethodGet:
		s.handleProjectInfo(w, r)
	case http.MethodPost:
		s.handleProjectCreate(w, r)
	default:
		http.Error(w, "method not allowed", http.StatusMethodNotAllowed)
	}
}

func (s *srv) handleProjectInfo(w http.ResponseWriter, r *http.Request) {
	q := r.URL.Query()
	cwd := q.Get("cwd")
	if cwd == "" {
		writeJSON(w, projectInfo{Required: []requiredResolution{}})
		return
	}
	cwd, err := filepath.Abs(cwd)
	if err != nil {
		http.Error(w, err.Error(), http.StatusBadRequest)
		return
	}
	standing, err := project.StandingAt(cwd)
	if err != nil {
		http.Error(w, err.Error(), http.StatusInternalServerError)
		return
	}
	info := projectInfo{CWD: cwd, Distro: config.ResolvedDefaultDistro(), Required: []requiredResolution{}}
	if standing != nil {
		info.Root = standing.Root
		if d := standing.SelectedDistro(); d != "" {
			info.Distro = d
		}
	}
	for _, name := range q["name"] {
		info.Required = append(info.Required, resolveRequired(r, standing, info.Distro, name))
	}
	writeJSON(w, info)
}

// resolveRequired resolves one name the way a launch does: through the
// project's pins when standing in one, by installed name otherwise. A name
// nothing answers comes back with Resolved empty.
func resolveRequired(r *http.Request, standing *project.Standing, distro, name string) requiredResolution {
	out := requiredResolution{Name: name}
	if standing != nil {
		request, reason := lock.ParseDeclaration(name)
		if reason != "" {
			return out
		}
		resolution, err := project.Resolve(r.Context(), standing.Root, standing.Lock,
			[]lock.Request{request}, project.ResolveOptions{})
		if err != nil || !resolution.Complete() {
			return out
		}
		mount := resolution.Mounts[0]
		out.Resolved = mount.Name
		out.Ident = shortIdentity(mount.Identity)
		return out
	}
	installed, err := image.ScanOverlays(image.ScanOptions{})
	if err != nil {
		return out
	}
	if resolved, _, found, err := catalog.SolveInstalled(r.Context(), installed, distro, name); err == nil && found {
		out.Resolved = resolved
	}
	return out
}

// handleProjectCreate creates a project at cwd with no scan: cnt-lock/ plus
// the default distro's base pin. A base that cannot be pinned still leaves the
// project created, and comes back as "warning".
func (s *srv) handleProjectCreate(w http.ResponseWriter, r *http.Request) {
	var req struct {
		CWD string `json:"cwd"`
	}
	if err := json.NewDecoder(r.Body).Decode(&req); err != nil || req.CWD == "" {
		http.Error(w, "cwd is required", http.StatusBadRequest)
		return
	}
	root, err := filepath.Abs(req.CWD)
	if err != nil {
		http.Error(w, err.Error(), http.StatusBadRequest)
		return
	}
	out := map[string]string{"root": root}
	if _, err := orchestrate.Init(r.Context(), root); err != nil {
		if _, statErr := os.Stat(lock.Dir(root)); statErr != nil {
			http.Error(w, err.Error(), http.StatusInternalServerError)
			return
		}
		out["warning"] = err.Error()
	}
	writeJSON(w, out)
}

type projectPin struct {
	Name  string `json:"name"`
	Path  string `json:"path"`
	Ident string `json:"ident"`
}

// handleProjectPins serves GET /api/project/pins?cwd= — the .sqf overlays
// pinned by the project cwd stands in that are available here, base excluded.
func (s *srv) handleProjectPins(w http.ResponseWriter, r *http.Request) {
	out := []projectPin{}
	cwd := r.URL.Query().Get("cwd")
	if cwd == "" {
		writeJSON(w, out)
		return
	}
	standing, err := project.StandingAt(cwd)
	if err != nil {
		http.Error(w, err.Error(), http.StatusInternalServerError)
		return
	}
	if standing == nil {
		writeJSON(w, out)
		return
	}
	var requests []lock.Request
	for key := range standing.Lock.Pins {
		if key == lock.BaseKey {
			continue
		}
		request, reason := lock.ParseDeclaration(strings.TrimPrefix(key, lock.PathPrefix))
		if reason != "" || request.Key != key {
			continue
		}
		requests = append(requests, request)
	}
	resolution, err := project.Resolve(r.Context(), standing.Root, standing.Lock, requests, project.ResolveOptions{})
	if err != nil {
		http.Error(w, err.Error(), http.StatusInternalServerError)
		return
	}
	for _, mount := range resolution.Mounts {
		if !utils.IsSqf(mount.Path) {
			continue
		}
		out = append(out, projectPin{Name: mount.Request, Path: mount.Path, Ident: shortIdentity(mount.Identity)})
	}
	sort.Slice(out, func(i, j int) bool { return out[i].Name < out[j].Name })
	writeJSON(w, out)
}

func shortIdentity(digest string) string {
	ident := strings.TrimPrefix(digest, "sha256:")
	return ident[:min(len(ident), capsule.IdentityChars)]
}
