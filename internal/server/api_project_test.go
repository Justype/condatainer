package server

import (
	"encoding/json"
	"net/http"
	"net/http/httptest"
	"net/url"
	"strings"
	"testing"
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
