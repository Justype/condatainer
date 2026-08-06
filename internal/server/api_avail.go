package server

import (
	"net/http"
	"sort"
	"strings"

	"github.com/Justype/condatainer/catalog"
	"github.com/Justype/condatainer/internal/config"
	"github.com/Justype/condatainer/internal/container"
	"github.com/Justype/condatainer/internal/scheduler"
	"github.com/Justype/condatainer/internal/utils"
)

// handleAvail serves GET /api/avail — lists available build scripts (local and
// remote, local wins on duplicates, mirroring `condatainer avail`). Template
// scripts are returned collapsed as a single entry carrying their placeholder
// metadata; expanded variants are suppressed.
func (s *srv) handleAvail(w http.ResponseWriter, r *http.Request) {
	type availEntry struct {
		Name           string              `json:"name"`
		Alias          string              `json:"alias,omitempty"` // bare name for default-distro scripts (ubuntu24/build-essential → build-essential)
		Whatis         string              `json:"whatis,omitempty"`
		Source         string              `json:"source,omitempty"`
		Container      bool                `json:"container,omitempty"`
		IsTemplate     bool                `json:"is_template,omitempty"`
		TargetTemplate string              `json:"target_template,omitempty"`
		PH             map[string][]string `json:"ph,omitempty"`
		PHNames        []string            `json:"ph_names,omitempty"`
		Installed      bool                `json:"installed,omitempty"`
	}

	installed := map[string]string{}
	if m, err := container.InstalledOverlays(); err == nil {
		installed = m
	}

	// Scripts under the default distro can be addressed by their bare name
	// (e.g. ubuntu24/build-essential → build-essential), matching `condatainer avail`.
	distroPrefix := ""
	if d := config.ResolvedBase(); d != "" {
		distroPrefix = d + "/"
	}

	entries := []availEntry{}
	cat, err := config.OpenCatalog(r.Context())
	if err != nil {
		http.Error(w, "failed to open recipe sources: "+err.Error(), http.StatusInternalServerError)
		return
	}
	seen := make(map[string]bool)
	for _, src := range cat {
		found, err := src.Entries(r.Context())
		if err != nil {
			continue
		}
		for name, e := range found {
			if seen[name] {
				continue
			}
			seen[name] = true
			_, isInstalled := installed[name]
			alias := ""
			if distroPrefix != "" {
				if a, ok := strings.CutPrefix(name, distroPrefix); ok && a != "" {
					alias = a
				}
			}
			var phNames []string
			if e.IsTemplate {
				phNames = catalog.NewTemplate(e.TargetTemplate).Names()
			}
			entries = append(entries, availEntry{
				Name:           name,
				Alias:          alias,
				Whatis:         e.Whatis,
				Source:         src.Name,
				Container:      strings.HasSuffix(e.Path, ".def"),
				IsTemplate:     e.IsTemplate,
				TargetTemplate: e.TargetTemplate,
				PH:             e.PH,
				PHNames:        phNames,
				Installed:      isInstalled,
			})
		}
	}

	sort.Slice(entries, func(i, j int) bool { return entries[i].Name < entries[j].Name })
	writeJSON(w, entries)
}

// handleAvailInspect serves GET /api/avail/inspect?name=<name/version> —
// pre-install info for one recipe: interactive prompts, scheduler directives,
// and whether the name would fall back to a conda install.
func (s *srv) handleAvailInspect(w http.ResponseWriter, r *http.Request) {
	name := catalog.Normalize(r.URL.Query().Get("name"))
	if name == "" {
		http.Error(w, "name required", http.StatusBadRequest)
		return
	}

	type inspectResp struct {
		Name                string   `json:"name"`
		Found               bool     `json:"found"`
		CondaFallback       bool     `json:"conda_fallback"`
		Prompts             []string `json:"prompts,omitempty"`
		SchedulerDirectives bool     `json:"scheduler_directives,omitempty"`
		Scheduler           string   `json:"scheduler,omitempty"`
		WillSubmit          bool     `json:"will_submit,omitempty"`
	}
	resp := inspectResp{Name: name}

	cat, err := config.OpenCatalog(r.Context())
	if err != nil {
		http.Error(w, "failed to open recipe sources: "+err.Error(), http.StatusInternalServerError)
		return
	}
	_, found, err := cat.Lookup(r.Context(), name)
	if err != nil {
		http.Error(w, "lookup failed: "+err.Error(), http.StatusBadGateway)
		return
	}
	resp.Found = found
	resp.CondaFallback = !found
	if !found {
		writeJSON(w, resp)
		return
	}

	// The recipe arrives expanded, so prompts and directives read straight off it.
	recipe, err := cat.Open(r.Context(), name, nil)
	if err != nil {
		http.Error(w, "failed to read recipe: "+err.Error(), http.StatusBadGateway)
		return
	}
	resp.Prompts = append(resp.Prompts, recipe.Inputs...)
	resp.SchedulerDirectives = len(recipe.Directives) > 0
	if sched := scheduler.ActiveScheduler(); sched != nil {
		resp.Scheduler = string(sched.GetType())
	}
	resp.WillSubmit = resp.SchedulerDirectives && resp.Scheduler != "" && config.Global.SubmitJob

	writeJSON(w, resp)
}

// handleSearch serves GET /api/search?q=<term> — exact conda package lookup on
// anaconda.org using the configured channels (first channel that has the
// package, matching install behaviour). Only an exact name match is returned;
// for fuzzy discovery, users go to anaconda.org directly.
func (s *srv) handleSearch(w http.ResponseWriter, r *http.Request) {
	q := r.URL.Query().Get("q")
	if q == "" {
		http.Error(w, "q required", http.StatusBadRequest)
		return
	}

	results, _, err := utils.SearchCondaPackages(q, config.Global.Build.Channels, false, 0)
	if err != nil {
		http.Error(w, "conda search failed: "+err.Error(), http.StatusBadGateway)
		return
	}
	if results == nil {
		results = []utils.CondaSearchResult{}
	}
	writeJSON(w, map[string]interface{}{
		"results":  results,
		"platform": utils.CurrentCondaPlatform(),
	})
}
