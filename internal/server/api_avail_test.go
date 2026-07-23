package server

import (
	"encoding/json"
	"net/http"
	"net/http/httptest"
	"net/url"
	"os"
	"path/filepath"
	"testing"

	"github.com/Justype/condatainer/internal/build"
	"github.com/Justype/condatainer/internal/config"
)

// setTestBuildScripts points the build-script scan at dir (remote disabled)
// and disables the on-disk scripts cache for the test.
func setTestBuildScripts(t *testing.T, dir string) {
	t.Helper()
	oldDirs := config.GlobalDataPaths.BuildScriptsDirs
	oldLinks := config.Global.ScriptsLinks
	oldRefresh := build.ForceRefresh
	oldDistro := config.Global.DefaultDistro
	config.GlobalDataPaths.BuildScriptsDirs = []string{dir}
	config.Global.ScriptsLinks = nil
	build.ForceRefresh = true
	if config.Global.DefaultDistro == "" {
		config.Global.DefaultDistro = "ubuntu24"
	}
	t.Cleanup(func() {
		config.GlobalDataPaths.BuildScriptsDirs = oldDirs
		config.Global.ScriptsLinks = oldLinks
		build.ForceRefresh = oldRefresh
		config.Global.DefaultDistro = oldDistro
	})
}

func writeScript(t *testing.T, dir, relPath, content string) {
	t.Helper()
	p := filepath.Join(dir, relPath)
	if err := os.MkdirAll(filepath.Dir(p), 0o775); err != nil {
		t.Fatal(err)
	}
	if err := os.WriteFile(p, []byte(content), 0o775); err != nil {
		t.Fatal(err)
	}
}

// TestHandleAvail checks that local scripts are listed, installed overlays are
// marked, and template scripts appear collapsed (variants suppressed).
func TestHandleAvail(t *testing.T) {
	scriptsDir := t.TempDir()
	writeScript(t, scriptsDir, "foo/1.0", "#!/usr/bin/bash\n#WHATIS:Test tool\ninstall_app() { :; }\n")
	writeScript(t, scriptsDir, "bar/gen",
		"#!/usr/bin/bash\n#PL:ver:2.0,1.0\n#TARGET:bar/{ver}\ninstall_app() { :; }\n")
	writeScript(t, scriptsDir, "ubuntu24/build-essential",
		"#!/usr/bin/bash\ninstall_app() { :; }\n")
	setTestBuildScripts(t, scriptsDir) // sets DefaultDistro=ubuntu24

	imgDir := t.TempDir()
	if err := os.WriteFile(filepath.Join(imgDir, "foo--1.0.sqf"), []byte("x"), 0o664); err != nil {
		t.Fatal(err)
	}
	setTestImagesDir(t, imgDir)

	r := httptest.NewRequest(http.MethodGet, "/api/avail", nil)
	w := httptest.NewRecorder()
	(&srv{}).handleAvail(w, r)
	if w.Code != http.StatusOK {
		t.Fatalf("status = %d (body: %s)", w.Code, w.Body.String())
	}

	type entry struct {
		Name       string              `json:"name"`
		Alias      string              `json:"alias"`
		Whatis     string              `json:"whatis"`
		IsTemplate bool                `json:"is_template"`
		PL         map[string][]string `json:"pl"`
		PLOrder    []string            `json:"pl_order"`
		Installed  bool                `json:"installed"`
	}
	var entries []entry
	if err := json.Unmarshal(w.Body.Bytes(), &entries); err != nil {
		t.Fatalf("bad JSON: %v", err)
	}

	byName := map[string]entry{}
	for _, e := range entries {
		byName[e.Name] = e
	}
	foo, ok := byName["foo/1.0"]
	if !ok {
		t.Fatalf("foo/1.0 missing from %v", byName)
	}
	if !foo.Installed || foo.Whatis != "Test tool" {
		t.Errorf("foo/1.0 = %+v, want installed with whatis", foo)
	}
	tmpl, ok := byName["bar/gen"]
	if !ok {
		t.Fatalf("template bar/gen missing from %v", byName)
	}
	if !tmpl.IsTemplate || len(tmpl.PLOrder) != 1 || tmpl.PLOrder[0] != "ver" || len(tmpl.PL["ver"]) != 2 {
		t.Errorf("template entry = %+v, want is_template with pl ver:[2.0 1.0]", tmpl)
	}
	// Expanded variants must be suppressed (collapsed view).
	for _, v := range []string{"bar/2.0", "bar/1.0"} {
		if _, exists := byName[v]; exists {
			t.Errorf("expanded variant %s should be suppressed", v)
		}
	}
	// Default-distro scripts expose a bare-name alias (ubuntu24/build-essential → build-essential).
	osEntry, ok := byName["ubuntu24/build-essential"]
	if !ok {
		t.Fatalf("ubuntu24/build-essential missing from %v", byName)
	}
	if osEntry.Alias != "build-essential" {
		t.Errorf("alias = %q, want %q", osEntry.Alias, "build-essential")
	}
}

// TestHandleAvailInspect checks prompt extraction and the conda fallback flag.
func TestHandleAvailInspect(t *testing.T) {
	scriptsDir := t.TempDir()
	writeScript(t, scriptsDir, "baz/1.0",
		"#!/usr/bin/bash\n#INTERACTIVE:Enter download URL\ninstall_app() { :; }\n")
	setTestBuildScripts(t, scriptsDir)
	setTestImagesDir(t, t.TempDir())

	get := func(name string) map[string]interface{} {
		r := httptest.NewRequest(http.MethodGet,
			"/api/avail/inspect?name="+url.QueryEscape(name), nil)
		w := httptest.NewRecorder()
		(&srv{}).handleAvailInspect(w, r)
		if w.Code != http.StatusOK {
			t.Fatalf("inspect %s: status = %d (body: %s)", name, w.Code, w.Body.String())
		}
		var resp map[string]interface{}
		if err := json.Unmarshal(w.Body.Bytes(), &resp); err != nil {
			t.Fatalf("bad JSON: %v", err)
		}
		return resp
	}

	resp := get("baz/1.0")
	if resp["found"] != true || resp["conda_fallback"] != false {
		t.Errorf("baz/1.0 = %v, want found without conda fallback", resp)
	}
	prompts, _ := resp["prompts"].([]interface{})
	if len(prompts) != 1 || prompts[0] != "Enter download URL" {
		t.Errorf("prompts = %v, want [Enter download URL]", prompts)
	}

	resp = get("nonexistent/9.9")
	if resp["found"] != false || resp["conda_fallback"] != true {
		t.Errorf("nonexistent = %v, want conda fallback", resp)
	}
}
