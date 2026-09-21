package server

import (
	"encoding/json"
	"net/http"
	"net/http/httptest"
	"net/url"
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/Justype/condatainer/internal/config"
	"github.com/Justype/condatainer/internal/project"
	"github.com/Justype/condatainer/internal/project/lock"
)

func TestProjectCreateThenInfo(t *testing.T) {
	dir := t.TempDir()
	s := &srv{}

	get := func() projectInfo {
		req := httptest.NewRequest(http.MethodGet, "/api/project?"+url.Values{"cwd": {dir}}.Encode(), nil)
		rec := httptest.NewRecorder()
		s.handleProject(rec, req)
		var info projectInfo
		if err := json.Unmarshal(rec.Body.Bytes(), &info); err != nil {
			t.Fatal(err)
		}
		return info
	}

	if info := get(); info.Root != "" {
		t.Fatalf("root before create = %q, want none", info.Root)
	}

	body := strings.NewReader(`{"cwd":"` + dir + `"}`)
	rec := httptest.NewRecorder()
	s.handleProject(rec, httptest.NewRequest(http.MethodPost, "/api/project", body))
	if rec.Code != http.StatusOK {
		t.Fatalf("create status = %d: %s", rec.Code, rec.Body)
	}

	if info := get(); info.Root != dir || info.Root != info.CWD {
		t.Fatalf("root after create = %q, cwd = %q", info.Root, info.CWD)
	}
}

// Inside a project, a required overlay the project has not pinned resolves by
// installed name, as the launch does, rather than showing as unresolved.
func TestResolveRequiredInAProjectUsesInstalled(t *testing.T) {
	imagesDir := t.TempDir()
	installed := "ubuntu24--rstudio-server--2026.09.0-174.sqf"
	if err := os.WriteFile(filepath.Join(imagesDir, installed), nil, 0o644); err != nil {
		t.Fatal(err)
	}
	prevPaths, prevDistro := config.GlobalDataPaths, config.Global.DefaultDistro
	config.GlobalDataPaths.ImagesDirs = []string{imagesDir}
	config.Global.DefaultDistro = "ubuntu24"
	t.Cleanup(func() { config.GlobalDataPaths, config.Global.DefaultDistro = prevPaths, prevDistro })

	root := t.TempDir()
	if err := lock.Publish(root, lock.New()); err != nil {
		t.Fatal(err)
	}
	standing, err := project.StandingAt(root)
	if err != nil || standing == nil {
		t.Fatalf("StandingAt: standing = %v, err = %v", standing, err)
	}
	req := httptest.NewRequest(http.MethodGet, "/api/project", nil)

	got := resolveRequired(req, standing, "ubuntu24", "rstudio-server")
	if got.Resolved != "ubuntu24/rstudio-server/2026.09.0-174" || got.Ident != "" {
		t.Errorf("an unpinned installed overlay resolved to %+v", got)
	}
	if got := resolveRequired(req, standing, "ubuntu24", "not-installed"); got.Resolved != "" {
		t.Errorf("an overlay nothing answers resolved to %+v", got)
	}
}
