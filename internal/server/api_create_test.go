package server

import (
	"context"
	"encoding/json"
	"net/http"
	"net/http/httptest"
	"os"
	"path/filepath"
	"strings"
	"testing"
	"time"
)

func postCreate(t *testing.T, s *srv, body string) *httptest.ResponseRecorder {
	t.Helper()
	r := httptest.NewRequest(http.MethodPost, "/api/create", strings.NewReader(body))
	r.Header.Set("Content-Type", "application/json")
	w := httptest.NewRecorder()
	s.handleCreate(w, r)
	return w
}

// waitForTaskDone subscribes to the task's broker and returns the parsed
// terminal "done" message.
func waitForTaskDone(t *testing.T, s *srv, taskID string) map[string]interface{} {
	t.Helper()
	v, ok := s.tasks.Load(taskID)
	if !ok {
		t.Fatalf("task %s not found", taskID)
	}
	broker := v.(*taskEntry).broker
	ch, finalMsg, finished := broker.subscribe()
	defer broker.unsubscribe(ch)

	parse := func(data []byte) (map[string]interface{}, bool) {
		var msg map[string]interface{}
		if err := json.Unmarshal(data, &msg); err != nil {
			return nil, false
		}
		return msg, msg["t"] == "done"
	}
	if finished {
		msg, _ := parse(finalMsg)
		return msg
	}
	deadline := time.After(10 * time.Second)
	for {
		select {
		case data := <-ch:
			if msg, done := parse(data); done {
				return msg
			}
		case <-deadline:
			t.Fatal("timed out waiting for task to finish")
		}
	}
}

// TestCreateValidation checks request validation before any task is started.
func TestCreateValidation(t *testing.T) {
	s := &srv{ctx: context.Background()}

	r := httptest.NewRequest(http.MethodGet, "/api/create", nil)
	w := httptest.NewRecorder()
	s.handleCreate(w, r)
	if w.Code != http.StatusMethodNotAllowed {
		t.Errorf("GET: status = %d, want 405", w.Code)
	}

	for _, body := range []string{"{not json", `{"name":""}`, `{"name":"  "}`, `{"name":"--update"}`} {
		if w := postCreate(t, s, body); w.Code != http.StatusBadRequest {
			t.Errorf("body %q: status = %d, want 400 (body: %s)", body, w.Code, w.Body.String())
		}
	}
}

// TestCreateRunsCLI stubs the executable and checks the subprocess flow:
// exit 0 → done ok, exit 3 (jobs submitted) → done ok + submitted.
func TestCreateRunsCLI(t *testing.T) {
	dir := t.TempDir()
	writeStub := func(name, script string) string {
		p := filepath.Join(dir, name)
		if err := os.WriteFile(p, []byte(script), 0o755); err != nil {
			t.Fatal(err)
		}
		return p
	}

	oldExe := createExecutable
	t.Cleanup(func() { createExecutable = oldExe })

	cases := []struct {
		name      string
		script    string
		submitted bool
		ok        bool
	}{
		{"exit0", "#!/bin/sh\necho building \"$@\"\nexit 0\n", false, true},
		{"exit3", "#!/bin/sh\necho submitted \"$@\"\nexit 3\n", true, true},
		{"exit1", "#!/bin/sh\necho failed \"$@\" >&2\nexit 1\n", false, false},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			stub := writeStub(tc.name+".sh", tc.script)
			createExecutable = func() (string, error) { return stub, nil }

			s := &srv{ctx: context.Background()}
			w := postCreate(t, s, `{"name":"foo/1.0"}`)
			if w.Code != http.StatusOK {
				t.Fatalf("status = %d, want 200 (body: %s)", w.Code, w.Body.String())
			}
			var resp map[string]string
			if err := json.Unmarshal(w.Body.Bytes(), &resp); err != nil || resp["id"] == "" {
				t.Fatalf("bad response: %s", w.Body.String())
			}

			msg := waitForTaskDone(t, s, resp["id"])
			if got := msg["ok"] == true; got != tc.ok {
				t.Errorf("ok = %v, want %v (msg: %v)", got, tc.ok, msg)
			}
			if got := msg["submitted"] == true; got != tc.submitted {
				t.Errorf("submitted = %v, want %v (msg: %v)", got, tc.submitted, msg)
			}
		})
	}
}
