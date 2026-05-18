package server

import (
	"context"
	"encoding/json"
	"errors"
	"fmt"
	"net/http"
	"os"
	"path/filepath"
	"strconv"
	"strings"
	"time"

	"github.com/Justype/condatainer/internal/config"
	"github.com/Justype/condatainer/internal/container"
	cntexec "github.com/Justype/condatainer/internal/exec"
	"github.com/Justype/condatainer/internal/helper"
	"github.com/Justype/condatainer/internal/overlay"
	"github.com/Justype/condatainer/internal/scheduler"
	"github.com/Justype/condatainer/internal/utils"
)

// handleHelpersAvailable serves GET /api/helpers/available
func (s *srv) handleHelpersAvailable(w http.ResponseWriter, r *http.Request) {
	type entry struct {
		Name             string               `json:"name"`
		Whatis           string               `json:"whatis,omitempty"`
		ImgPackages      string               `json:"img_packages,omitempty"`
		PostInstallCmd   string               `json:"post_install_cmd,omitempty"`
		RequiredOverlays string               `json:"required_overlays,omitempty"`
		Singleton        bool                 `json:"singleton,omitempty"`
		Params           []helper.HelperParam `json:"params,omitempty"`
		ParamValues      map[string][]string  `json:"param_values,omitempty"`
	}

	seen := make(map[string]bool)
	var scripts []entry

	for _, dir := range config.GetHelperScriptSearchPaths() {
		dirEntries, err := os.ReadDir(dir)
		if err != nil {
			continue
		}
		for _, e := range dirEntries {
			name := e.Name()
			if !e.IsDir() && !strings.HasPrefix(name, ".") && !seen[name] {
				seen[name] = true
				scriptPath := filepath.Join(dir, name)
				whatis := ""
				if data, err := os.ReadFile(scriptPath); err == nil {
					for _, line := range strings.Split(string(data), "\n") {
						if strings.HasPrefix(line, "#WHATIS:") {
							whatis = strings.TrimSpace(line[8:])
							break
						}
					}
				}
				meta, _ := helper.ParseHelperScriptMeta(scriptPath)
				params, _ := helper.ParseHelperParams(scriptPath)
				ent := entry{
					Name:           name,
					Whatis:         whatis,
					ImgPackages:    meta.ImgPackages,
					PostInstallCmd: meta.PostInstallCmd,
					Singleton:      meta.Singleton,
					Params:         params,
					ParamValues:    meta.ParamValues,
				}
				if meta.RequiredOverlays != "" {
					ent.RequiredOverlays = meta.RequiredOverlays
				}
				scripts = append(scripts, ent)
			}
		}
	}

	if r.Header.Get("Accept") == "text/html" || r.URL.Query().Get("format") == "html" {
		w.Header().Set("Content-Type", "text/html")
		fmt.Fprint(w, `<div class="helper-grid">`)
		for _, sc := range scripts {
			needsOverlay := sc.ImgPackages != ""
			fmt.Fprintf(w, `<div class="helper-item" data-name="%s" data-img-packages="%s" onclick="selectHelper('%s')"><strong>%s</strong>`,
				sc.Name, sc.ImgPackages, sc.Name, sc.Name)
			if sc.Whatis != "" {
				fmt.Fprintf(w, `<br><small style="color:var(--muted)">%s</small>`, sc.Whatis)
			}
			if needsOverlay {
				fmt.Fprint(w, `<br><span class="badge badge-overlay">overlay required</span>`)
			}
			fmt.Fprint(w, `</div>`)
		}
		fmt.Fprint(w, `</div>`)
		return
	}
	writeJSON(w, scripts)
}

// handleHelpersSub routes /api/helpers/{id-or-name}/{action}.
// Some actions key off the helper run ID (logs, joblog, stop) and others off
// the helper script name (start, find-env, missing-params, running, resources).
func (s *srv) handleHelpersSub(w http.ResponseWriter, r *http.Request) {
	sub := strings.TrimPrefix(r.URL.Path, "/api/helpers/")
	parts := strings.SplitN(sub, "/", 2)
	if len(parts) < 2 {
		http.NotFound(w, r)
		return
	}
	id, action := parts[0], parts[1]

	switch action {
	case "logs":
		s.handleHelperLogs(w, r, id)
	case "joblog":
		s.handleHelperJobLog(w, r, id)
	case "stop":
		s.handleHelperStop(w, r, id)
	case "delete":
		s.handleHelperDelete(w, r, id)
	case "start":
		s.handleHelperStart(w, r, id)
	case "find-env":
		s.handleHelperOverlay(w, r, id)
	case "missing-params":
		s.handleHelperMissingParams(w, r, id)
	case "running":
		s.handleHelperRunning(w, r, id)
	case "resources":
		s.handleHelperResources(w, r, id)
	default:
		http.NotFound(w, r)
	}
}

// handleEnvCheck serves GET /api/env/check?path=...&name=...&cwd=...
func (s *srv) handleEnvCheck(w http.ResponseWriter, r *http.Request) {
	q := r.URL.Query()
	path := q.Get("path")
	if path == "" {
		path = helper.ResolveEnv(q.Get("cwd"))
	}
	st, err := helper.CheckEnv(r.Context(), path)
	if err != nil {
		http.Error(w, err.Error(), http.StatusInternalServerError)
		return
	}
	writeJSON(w, st)
}

// handleHelperMissingParams serves GET /api/helpers/{name}/missing-params
func (s *srv) handleHelperMissingParams(w http.ResponseWriter, r *http.Request, name string) {
	scriptPath, err := config.FindHelperScript(name)
	if err != nil {
		http.Error(w, "helper not found: "+name, http.StatusNotFound)
		return
	}
	supplied := make(map[string]string)
	for k, v := range r.URL.Query() {
		if len(v) > 0 {
			supplied[k] = v[0]
		}
	}
	missing, err := helper.MissingParams(scriptPath, supplied)
	if err != nil {
		http.Error(w, err.Error(), http.StatusInternalServerError)
		return
	}
	writeJSON(w, missing)
}

// handleHelperRunning serves GET /api/helpers/{name}/running
func (s *srv) handleHelperRunning(w http.ResponseWriter, r *http.Request, name string) {
	running, err := helper.Running(name)
	if err != nil {
		http.Error(w, err.Error(), http.StatusInternalServerError)
		return
	}
	type resp struct {
		Running          []*helper.HelperRun `json:"running"`
		SingletonBlocked bool                `json:"singleton_blocked"`
	}
	out := resp{Running: running}
	if scriptPath, err := config.FindHelperScript(name); err == nil {
		if meta, err := helper.ParseHelperScriptMeta(scriptPath); err == nil {
			out.SingletonBlocked = helper.SingletonBlocked(meta, running)
		}
	}
	writeJSON(w, out)
}

// handleHelperResources serves GET /api/helpers/{name}/resources?cpus=&mem=&time=
func (s *srv) handleHelperResources(w http.ResponseWriter, r *http.Request, name string) {
	scriptPath, err := config.FindHelperScript(name)
	if err != nil {
		http.Error(w, "helper not found: "+name, http.StatusNotFound)
		return
	}
	q := r.URL.Query()
	cpus := 0
	if v := q.Get("cpus"); v != "" {
		if n, err := strconv.Atoi(v); err == nil {
			cpus = n
		}
	}
	var overrides *scheduler.ResourceSpec
	memStr := q.Get("mem")
	timeStr := q.Get("time")
	gpuStr := q.Get("gpu")
	if cpus > 0 || memStr != "" || timeStr != "" || gpuStr != "" {
		overrides = &scheduler.ResourceSpec{}
		if cpus > 0 {
			overrides.CpusPerTask = cpus
		}
		if memStr != "" {
			if mb, err := utils.ParseMemoryMB(memStr); err == nil {
				overrides.MemPerNodeMB = mb
			}
		}
		if timeStr != "" {
			if d, err := utils.ParseWalltime(timeStr); err == nil {
				overrides.Time = d
			}
		}
		if gpuStr != "" {
			overrides.Gpu = helper.ParseGPUSpec(gpuStr)
		}
	}
	spec := helper.ResolveResources(r.Context(), scriptPath, overrides)
	type resp struct {
		CPUs     int    `json:"cpus"`
		MemMB    int64  `json:"mem_mb"`
		Mem      string `json:"mem"`
		Walltime string `json:"walltime"`
		GPU      string `json:"gpu,omitempty"`
	}
	out := resp{}
	if spec != nil {
		out.CPUs = spec.GetNtasks() * spec.CpusPerTask
		out.MemMB = spec.GetMemPerNodeMB()
		if out.MemMB > 0 {
			out.Mem = utils.FormatMemoryMB(out.MemMB)
		}
		if spec.Time > 0 {
			out.Walltime = utils.FormatDuration(spec.Time)
		}
		out.GPU = helper.FormatGpuSpec(spec)
	}
	writeJSON(w, out)
}

// handleHelperOverlay serves GET /api/helpers/{name}/find-env
// cwd is required; if empty, returns no path rather than falling back to the
// server's working directory, which is unrelated to the user's job directory.
func (s *srv) handleHelperOverlay(w http.ResponseWriter, r *http.Request, name string) {
	cwd := r.URL.Query().Get("cwd")
	var path string
	if cwd != "" {
		path = helper.ResolveEnvOverlayInDir("", cwd)
	}
	w.Header().Set("Content-Type", "application/json")
	json.NewEncoder(w).Encode(map[string]string{"path": path}) //nolint:errcheck
}

// handleHelperStart serves POST /api/helpers/{name}/start — launches a new helper job.
func (s *srv) handleHelperStart(w http.ResponseWriter, r *http.Request, name string) {
	if r.Method != http.MethodPost {
		http.Error(w, "POST required", http.StatusMethodNotAllowed)
		return
	}

	type startReq struct {
		CPUs     int               `json:"cpus"`
		Mem      string            `json:"mem"`
		Time     string            `json:"time"`
		GPU      string            `json:"gpu"`
		CWD      string            `json:"cwd"`
		Overlay  string            `json:"overlay"`
		Overlays []string          `json:"overlays"`
		Params   map[string]string `json:"params"`
	}
	var req startReq
	if err := json.NewDecoder(r.Body).Decode(&req); err != nil {
		http.Error(w, "bad request body", http.StatusBadRequest)
		return
	}

	scriptPath, err := config.FindHelperScript(name)
	if err != nil {
		http.Error(w, "helper not found: "+name, http.StatusNotFound)
		return
	}

	var startOverrides *scheduler.ResourceSpec
	if req.CPUs > 0 || req.Mem != "" || req.Time != "" || req.GPU != "" {
		startOverrides = &scheduler.ResourceSpec{}
		if req.CPUs > 0 {
			startOverrides.CpusPerTask = req.CPUs
		}
		if req.Mem != "" {
			if mb, err := utils.ParseMemoryMB(req.Mem); err == nil {
				startOverrides.MemPerNodeMB = mb
			}
		}
		if req.Time != "" {
			if d, err := utils.ParseWalltime(req.Time); err == nil {
				startOverrides.Time = d
			}
		}
		if req.GPU != "" {
			startOverrides.Gpu = helper.ParseGPUSpec(req.GPU)
		}
	}

	opts := helper.RunOptions{
		ScriptPath: scriptPath,
		ScriptName: name,
		Resources:  startOverrides,
		CWD:        req.CWD,
		Overlays:   req.Overlays,
		Params:     req.Params,
		ForceNew:   true,
	}
	if req.Overlay != "" {
		opts.EnvImg = req.Overlay
	}

	plan, err := helper.PlanRun(s.ctx, opts)
	if err != nil {
		writePlanError(w, err)
		return
	}

	w.Header().Set("Content-Type", "text/plain; charset=utf-8")
	flusher, _ := w.(http.Flusher)
	fmt.Fprintf(w, "Starting %s...\n", name)
	if flusher != nil {
		flusher.Flush()
	}

	id, err := helper.ExecutePlan(s.ctx, plan, helper.IO{Stdout: w, Stderr: w})
	if err != nil {
		fmt.Fprintf(w, "ERROR: %v\n", err)
		return
	}
	fmt.Fprintf(w, "Submitted: %s\n", id)
}

// handleOverlayCreate serves POST /api/overlay/create.
// It validates the request synchronously, then starts a background task and
// immediately returns {"id": "<taskID>"}.  The caller should open
// GET /api/tasks/<id>/stream to receive SSE progress events.
func (s *srv) handleOverlayCreate(w http.ResponseWriter, r *http.Request) {
	if r.Method != http.MethodPost {
		http.Error(w, "POST required", http.StatusMethodNotAllowed)
		return
	}

	type createReq struct {
		Name          string `json:"name"`
		Size          string `json:"size"`
		BasePackages  string `json:"base_packages"`
		ExtraPackages string `json:"extra_packages"`
		BaseImage     string `json:"base_image"`
		PostInstall   string `json:"post_install_cmd"`
		CWD           string `json:"cwd"`
	}
	var req createReq
	if err := json.NewDecoder(r.Body).Decode(&req); err != nil {
		http.Error(w, "bad request body", http.StatusBadRequest)
		return
	}
	if req.Name == "" {
		http.Error(w, "name required", http.StatusBadRequest)
		return
	}

	imgPath := req.Name
	if !filepath.IsAbs(imgPath) {
		cwd := req.CWD
		if cwd == "" {
			cwd, _ = os.Getwd()
		}
		imgPath = filepath.Join(cwd, imgPath)
	}

	sizeStr := req.Size
	if sizeStr == "" {
		sizeStr = "20G"
	}
	sizeMB, err := utils.ParseSizeToMB(sizeStr)
	if err != nil {
		http.Error(w, "invalid size: "+err.Error(), http.StatusBadRequest)
		return
	}

	allPkgs := strings.Fields(req.BasePackages)
	if req.ExtraPackages != "" {
		allPkgs = append(allPkgs, strings.Fields(req.ExtraPackages)...)
	}
	if len(allPkgs) == 0 {
		http.Error(w, "no packages specified", http.StatusBadRequest)
		return
	}

	// Create the background task and return its ID immediately.
	taskID := fmt.Sprintf("overlay-%d", time.Now().UnixNano())
	broker := newSSEBroker()
	ctx, cancel := context.WithCancel(s.ctx)
	s.tasks.Store(taskID, &taskEntry{broker: broker, cancel: cancel})
	writeJSON(w, map[string]string{"id": taskID})

	bw := &brokerWriter{broker}
	go func() {
		defer cancel()
		defer s.tasks.Delete(taskID)

		fmt.Fprintf(bw, "Creating overlay at %s (%s)...\n", imgPath, sizeStr)
		if err := os.MkdirAll(filepath.Dir(imgPath), utils.PermDir); err != nil {
			broadcastDone(broker, fmt.Errorf("mkdir: %w", err))
			return
		}
		opts := &overlay.CreateOptions{
			Path:    imgPath,
			SizeMB:  sizeMB,
			Profile: overlay.ProfileSmall,
		}
		io := cntexec.IO{Stdout: bw, Stderr: bw}
		fmt.Fprintf(bw, "Installing packages: %s\n", strings.Join(allPkgs, " "))
		if err := cntexec.CreateCondaOverlay(ctx, opts, allPkgs, req.PostInstall, false, io); err != nil {
			broadcastDone(broker, fmt.Errorf("create overlay: %w", err))
			return
		}

		container.InvalidateInstalledOverlaysCache()
		fmt.Fprintf(bw, "Done. Overlay ready at: %s\n", imgPath)
		// Signal completion with the result path so the client can populate cfg-overlay.
		result, _ := json.Marshal(map[string]interface{}{"t": "done", "ok": true, "path": imgPath})
		broker.publish(result)
	}()
}

// writePlanError translates helper.PlanRun errors into typed JSON HTTP responses.
func writePlanError(w http.ResponseWriter, err error) {
	w.Header().Set("Content-Type", "application/json; charset=utf-8")
	var missing *helper.ErrMissingParam
	if errors.As(err, &missing) {
		w.WriteHeader(http.StatusBadRequest)
		keys := make([]string, len(missing.Params))
		for i, p := range missing.Params {
			keys[i] = p.Key
		}
		_ = json.NewEncoder(w).Encode(map[string]interface{}{
			"error":          "missing_params",
			"missing_params": keys,
			"detail":         missing.Error(),
		})
		return
	}
	var singleton *helper.ErrSingletonRunning
	if errors.As(err, &singleton) {
		w.WriteHeader(http.StatusConflict)
		ids := make([]string, len(singleton.Existing))
		for i, r := range singleton.Existing {
			ids[i] = r.ID
		}
		_ = json.NewEncoder(w).Encode(map[string]interface{}{
			"error":   "singleton_running",
			"running": ids,
			"detail":  singleton.Error(),
		})
		return
	}
	var inUse *helper.ErrEnvInUse
	if errors.As(err, &inUse) {
		w.WriteHeader(http.StatusConflict)
		_ = json.NewEncoder(w).Encode(map[string]interface{}{
			"error":  "env_in_use",
			"path":   inUse.Path,
			"detail": inUse.Error(),
		})
		return
	}
	w.WriteHeader(http.StatusInternalServerError)
	_ = json.NewEncoder(w).Encode(map[string]string{"error": err.Error()})
}
