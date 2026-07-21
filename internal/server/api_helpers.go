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
	"sync"
	"time"

	"log/slog"

	"github.com/Justype/condatainer/internal/config"
	"github.com/Justype/condatainer/internal/container"
	cntexec "github.com/Justype/condatainer/internal/exec"
	"github.com/Justype/condatainer/internal/helper"
	"github.com/Justype/condatainer/internal/logging"
	"github.com/Justype/condatainer/internal/overlay"
	"github.com/Justype/condatainer/internal/scheduler"
	"github.com/Justype/condatainer/internal/utils"
	"github.com/Justype/condatainer/internal/weblog"
)

const helperAvailableCacheTTL = 30 * time.Second

// flushWriter wraps an http.ResponseWriter and flushes after every Write call.
type flushWriter struct {
	w     http.ResponseWriter
	flush func()
}

func (fw flushWriter) Write(p []byte) (int, error) {
	n, err := fw.w.Write(p)
	fw.flush()
	return n, err
}

type helperAvailableEntry struct {
	Name             string               `json:"name"`
	Whatis           string               `json:"whatis,omitempty"`
	ImgRequired      bool                 `json:"img_required,omitempty"`
	ImgPackages      string               `json:"img_packages,omitempty"`
	PostInstallCmd   string               `json:"post_install_cmd,omitempty"`
	RequiredOverlays string               `json:"required_overlays,omitempty"`
	Singleton        bool                 `json:"singleton,omitempty"`
	Params           []helper.HelperParam `json:"params,omitempty"`
	ParamValues      map[string][]string  `json:"param_values,omitempty"`
}

var (
	helperAvailableCacheMu      sync.Mutex
	helperAvailableCacheEntries []helperAvailableEntry
	helperAvailableCacheAt      time.Time
)

func cachedAvailableHelpers() []helperAvailableEntry {
	helperAvailableCacheMu.Lock()
	defer helperAvailableCacheMu.Unlock()
	if helperAvailableCacheEntries != nil && time.Since(helperAvailableCacheAt) < helperAvailableCacheTTL {
		return helperAvailableCacheEntries
	}

	seen := make(map[string]bool)
	var scripts []helperAvailableEntry

	for _, dir := range config.GetHelperScriptSearchPaths() {
		dirEntries, err := os.ReadDir(dir)
		if err != nil {
			continue
		}
		for _, e := range dirEntries {
			name := e.Name()
			if e.IsDir() || strings.HasPrefix(name, ".") || seen[name] {
				continue
			}
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
			ent := helperAvailableEntry{
				Name:           name,
				Whatis:         whatis,
				ImgRequired:    meta.ImgRequired,
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

	helperAvailableCacheEntries = scripts
	helperAvailableCacheAt = time.Now()
	return scripts
}

func invalidateHelperAvailableCache() {
	helperAvailableCacheMu.Lock()
	helperAvailableCacheEntries = nil
	helperAvailableCacheAt = time.Time{}
	helperAvailableCacheMu.Unlock()
}

// handleHelpersAvailable serves GET /api/helpers/available
func (s *srv) handleHelpersAvailable(w http.ResponseWriter, r *http.Request) {
	if r.URL.Query().Get("refresh") == "1" {
		invalidateHelperAvailableCache()
	}
	scripts := cachedAvailableHelpers()

	writeJSON(w, scripts)
}

// handleHelpersUpdate serves POST /api/helpers/update — syncs helper scripts from remote sources.
func (s *srv) handleHelpersUpdate(w http.ResponseWriter, r *http.Request) {
	if r.Method != http.MethodPost {
		http.Error(w, "POST required", http.StatusMethodNotAllowed)
		return
	}
	type updateReq struct {
		Name string `json:"name"`
	}
	var req updateReq
	if r.Body != nil {
		_ = json.NewDecoder(r.Body).Decode(&req)
	}

	taskID := fmt.Sprintf("helpers-update-%d", time.Now().UnixNano())
	broker := newSSEBroker()
	ctx, cancel := context.WithCancel(s.ctx)
	s.tasks.Store(taskID, &taskEntry{broker: broker, cancel: cancel})
	writeJSON(w, map[string]string{"id": taskID})

	bw := &brokerWriter{broker}
	go func() {
		defer cancel()
		defer s.scheduleTaskCleanup(taskID)

		if err := helper.UpdateRemoteScripts(ctx, req.Name, false, bw); err != nil {
			broadcastResult(broker, ctx, err)
			return
		}
		invalidateHelperAvailableCache()
		broadcastDone(broker, nil)
	}()
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
	case "config":
		s.handleHelperConfig(w, r, id)
	default:
		http.NotFound(w, r)
	}
}

// handleEnvInfo serves GET /api/env/info?path=... — returns size and explicitly-installed
// conda specs for an existing .img overlay (read from conda-meta/history).
func (s *srv) handleEnvInfo(w http.ResponseWriter, r *http.Request) {
	path := r.URL.Query().Get("path")
	if path == "" {
		http.Error(w, "path required", http.StatusBadRequest)
		return
	}
	type envInfoResp struct {
		SizeMB   int64    `json:"size_mb"`
		Channels []string `json:"channels"`
		Specs    []string `json:"specs"`
	}
	var resp envInfoResp
	if info, err := os.Stat(path); err == nil {
		resp.SizeMB = info.Size() / (1024 * 1024)
	}
	if ci := overlay.ReadCondaInfo(path, "/ext3/env"); ci != nil {
		resp.Channels = ci.Channels
		resp.Specs = ci.Specs
	}
	writeJSON(w, resp)
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

// handleHelperConfig serves GET/POST /api/helpers/{name}/config.
// GET returns the saved per-helper key=value config as JSON.
// POST accepts a JSON object of updates; empty string value deletes the key.
func (s *srv) handleHelperConfig(w http.ResponseWriter, r *http.Request, name string) {
	switch r.Method {
	case http.MethodGet:
		cfg, _ := helper.LoadHelperConfig(name)
		if cfg == nil {
			cfg = map[string]string{}
		}
		writeJSON(w, cfg)
	case http.MethodPost:
		var updates map[string]string
		if err := json.NewDecoder(r.Body).Decode(&updates); err != nil {
			http.Error(w, "bad request body", http.StatusBadRequest)
			return
		}
		cfg, _ := helper.LoadHelperConfig(name)
		if cfg == nil {
			cfg = map[string]string{}
		}
		for k, v := range updates {
			lk := strings.ToLower(k)
			if v == "" {
				delete(cfg, lk)
			} else {
				cfg[lk] = v
			}
		}
		if err := helper.SaveHelperConfig(name, cfg); err != nil {
			http.Error(w, err.Error(), http.StatusInternalServerError)
			return
		}
		writeJSON(w, cfg)
	default:
		http.Error(w, "method not allowed", http.StatusMethodNotAllowed)
	}
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

	params := req.Params
	if saved, err := helper.LoadHelperConfig(name); err == nil && len(saved) > 0 {
		if params == nil {
			params = make(map[string]string)
		}
		// Saved-config keys are lowercased; param keys keep their original case
		// and PlanRun matches case-sensitively. Map each saved value onto the
		// actual param key so set defaults aren't lost to the case mismatch.
		scriptParams, _ := helper.ParseHelperParams(scriptPath)
		for _, p := range scriptParams {
			if v, ok := saved[strings.ToLower(p.Key)]; ok && v != "" {
				if _, exists := params[p.Key]; !exists {
					params[p.Key] = v
				}
			}
		}
	}

	opts := helper.RunOptions{
		ScriptPath: scriptPath,
		ScriptName: name,
		Resources:  startOverrides,
		CWD:        req.CWD,
		Overlays:   req.Overlays,
		Params:     params,
		ForceNew:   true,
	}
	if req.Overlay != "" {
		opts.EnvImg = req.Overlay
	}

	w.Header().Set("Content-Type", "text/plain; charset=utf-8")
	flusher, _ := w.(http.Flusher)
	flush := func() {
		if flusher != nil {
			flusher.Flush()
		}
	}
	fmt.Fprintf(w, "Preparing %s...\n", name)
	flush()

	// flushW flushes after every write so log lines reach the browser immediately.
	flushW := flushWriter{w: w, flush: flush}
	planCtx := logging.WithLogger(s.ctx, slog.New(weblog.New(flushW)))
	planCtx = logging.WithWriter(planCtx, flushW)

	plan, err := helper.PlanRun(planCtx, opts)
	if err != nil {
		flush()
		fmt.Fprintf(w, "ERROR: %s\n", planErrMsg(err))
		flush()
		return
	}

	id, err := helper.ExecutePlan(s.ctx, plan)
	if err != nil {
		fmt.Fprintf(w, "ERROR: %v\n", err)
		flush()
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

	// Abort early if the target file already exists to prevent accidental overwrites.
	if _, err := os.Stat(imgPath); err == nil {
		w.Header().Set("Content-Type", "application/json; charset=utf-8")
		w.WriteHeader(http.StatusConflict)
		json.NewEncoder(w).Encode(map[string]string{ //nolint:errcheck
			"error":  "file_exists",
			"detail": imgPath + " already exists; delete or rename it first",
		})
		return
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

	// Create the background task and return its ID immediately.
	taskID := fmt.Sprintf("overlay-%d", time.Now().UnixNano())
	broker := newSSEBroker()
	ctx, cancel := context.WithCancel(s.ctx)
	s.tasks.Store(taskID, &taskEntry{broker: broker, cancel: cancel})
	writeJSON(w, map[string]string{"id": taskID})

	bw := &brokerWriter{broker}
	ctx = logging.WithLogger(ctx, slog.New(weblog.New(bw)))
	ctx = logging.WithWriter(ctx, bw)
	go func() {
		defer cancel()
		defer s.scheduleTaskCleanup(taskID)

		fmt.Fprintf(bw, "Creating overlay at %s (%s)...\n", imgPath, sizeStr)
		if err := utils.MkdirAllShared(filepath.Dir(imgPath)); err != nil {
			broadcastDone(broker, fmt.Errorf("mkdir: %w", err))
			return
		}
		opts := &overlay.CreateOptions{
			Path:    imgPath,
			SizeMB:  sizeMB,
			Profile: overlay.ProfileSmall,
		}
		io := cntexec.IO{Stdout: bw, Stderr: bw}
		fmt.Fprintf(bw, "Installing packages: %s\n", cntexec.DescribeInitialCondaPackages(allPkgs))
		if err := cntexec.CreateCondaOverlay(ctx, opts, allPkgs, req.PostInstall, false, io); err != nil {
			broadcastResult(broker, ctx, fmt.Errorf("create overlay: %w", err))
			return
		}

		container.InvalidateInstalledOverlaysCache()
		fmt.Fprintf(bw, "Done. Overlay ready at: %s\n", imgPath)
		// Signal completion with the result path so the client can populate cfg-overlay.
		result, _ := json.Marshal(map[string]interface{}{"t": "done", "ok": true, "path": imgPath})
		broker.publishFinal(result)
	}()
}

// planErrMsg converts a PlanRun error into a human-readable string for the
// streaming start terminal (used when headers are already committed).
func planErrMsg(err error) string {
	var missing *helper.ErrMissingParam
	if errors.As(err, &missing) {
		keys := make([]string, len(missing.Params))
		for i, p := range missing.Params {
			keys[i] = p.Key
		}
		return "missing required param(s): " + strings.Join(keys, ", ")
	}
	var singleton *helper.ErrSingletonRunning
	if errors.As(err, &singleton) {
		return fmt.Sprintf("already running (%d instance(s)) — stop it first or open the existing session", len(singleton.Existing))
	}
	var inUse *helper.ErrEnvInUse
	if errors.As(err, &inUse) {
		return "env overlay is in use by another running instance — stop it first or choose a different overlay"
	}
	return err.Error()
}
