package server

import (
	"net/http"
	"sort"
	"strings"

	"github.com/Justype/condatainer/internal/build"
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
		Alias          string              `json:"alias,omitempty"` // bare name for default-distro scripts (ubuntu24/igv → igv)
		Whatis         string              `json:"whatis,omitempty"`
		Remote         bool                `json:"remote,omitempty"`
		Container      bool                `json:"container,omitempty"`
		IsTemplate     bool                `json:"is_template,omitempty"`
		TargetTemplate string              `json:"target_template,omitempty"`
		PL             map[string][]string `json:"pl,omitempty"`
		PLOrder        []string            `json:"pl_order,omitempty"`
		Installed      bool                `json:"installed,omitempty"`
	}

	installed := map[string]string{}
	if m, err := container.InstalledOverlays(); err == nil {
		installed = m
	}

	// Scripts under the default distro can be addressed by their bare name
	// (e.g. ubuntu24/igv → igv), matching `condatainer avail`.
	distroPrefix := ""
	if d := config.Global.DefaultDistro; d != "" {
		distroPrefix = d + "/"
	}

	entries := []availEntry{}
	seen := make(map[string]bool)
	add := func(scripts map[string]build.ScriptInfo, remote bool) {
		for name, info := range scripts {
			if seen[name] {
				continue
			}
			seen[name] = true
			// Suppress expanded template variants: the collapsed template entry
			// carries the full placeholder metadata instead.
			if !info.IsTemplate && len(info.PLOrder) > 0 {
				continue
			}
			_, isInstalled := installed[name]
			alias := ""
			if distroPrefix != "" {
				if a, ok := strings.CutPrefix(name, distroPrefix); ok && a != "" {
					alias = a
				}
			}
			entries = append(entries, availEntry{
				Name:           name,
				Alias:          alias,
				Whatis:         info.Whatis,
				Remote:         remote,
				Container:      info.IsContainer,
				IsTemplate:     info.IsTemplate,
				TargetTemplate: info.TargetTemplate,
				PL:             info.PL,
				PLOrder:        info.PLOrder,
				Installed:      isInstalled,
			})
		}
	}

	if local, err := build.GetLocalBuildScripts(); err == nil {
		add(local, false)
	}
	if remote, err := build.GetRemoteBuildScripts(); err == nil {
		add(remote, true)
	}

	sort.Slice(entries, func(i, j int) bool { return entries[i].Name < entries[j].Name })
	writeJSON(w, entries)
}

// handleAvailInspect serves GET /api/avail/inspect?name=<name/version> —
// pre-install info for one build script: interactive prompts, scheduler
// directives, and whether the name would fall back to a conda install.
// Remote scripts are downloaded to the tmp dir so their content can be parsed.
func (s *srv) handleAvailInspect(w http.ResponseWriter, r *http.Request) {
	name := utils.NormalizeNameVersion(r.URL.Query().Get("name"))
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

	info, found := build.FindBuildScript(name)
	resp.Found = found
	resp.CondaFallback = !found
	if !found {
		writeJSON(w, resp)
		return
	}

	scriptPath := info.Path
	if info.IsRemote {
		p, err := build.DownloadRemoteScript(info, config.GetWritableTmpDir())
		if err != nil {
			http.Error(w, "failed to download remote script: "+err.Error(), http.StatusBadGateway)
			return
		}
		scriptPath = p
	}

	if prompts, err := utils.GetInteractivePromptsFromScript(scriptPath); err == nil {
		vars := info.CurrentVars()
		for _, p := range prompts {
			if len(vars) > 0 {
				p = utils.InterpolateVars(p, vars)
			}
			resp.Prompts = append(resp.Prompts, p)
		}
	}

	if specs, err := scheduler.ReadScriptSpecsFromPath(scriptPath); err == nil && specs != nil {
		resp.SchedulerDirectives = specs.HasDirectives
	}
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
