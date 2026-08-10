package build

import (
	"strings"
	"testing"

	"github.com/Justype/condatainer/internal/artifact/meta"
)

func TestParseBootstrap(t *testing.T) {
	tests := []struct {
		name        string
		def         string
		agent, from string
	}{
		{"plain", "Bootstrap: docker\nFrom: ubuntu:24.04\n", "docker", "ubuntu:24.04"},
		{"case and spacing", "bootstrap:  Docker \nFROM:\tubuntu:24.04\n", "docker", "ubuntu:24.04"},
		{"headers before", "#DESC:x\nBootstrap: docker\nFrom: alpine:3.20\n", "docker", "alpine:3.20"},
		{"no upstream", "Bootstrap: scratch\n", "scratch", ""},
		// A From: inside a section is shell text, not a directive.
		{"from in a section", "Bootstrap: docker\nFrom: ubuntu:24.04\n%post\nFrom: not-a-directive\n", "docker", "ubuntu:24.04"},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			got := parseBootstrap([]byte(tt.def))
			if got.Agent != tt.agent || got.From != tt.from {
				t.Errorf("got %+v, want {%s %s}", got, tt.agent, tt.from)
			}
		})
	}
}

func TestBootstrapClassification(t *testing.T) {
	tests := []struct {
		def                string
		resolvable, remote bool
	}{
		{"Bootstrap: docker\nFrom: ubuntu:24.04\n", true, true},
		{"Bootstrap: oras\nFrom: reg.example.com/img:1\n", true, true},
		// Remote and moving, but not an OCI registry: unresolved, not absent.
		{"Bootstrap: library\nFrom: ubuntu:24.04\n", false, true},
		{"Bootstrap: scratch\n", false, false},
		{"Bootstrap: localimage\nFrom: /tmp/x.sif\n", false, false},
		{"Bootstrap: debootstrap\nOSVersion: noble\n", false, false},
		// An agent that pulls, with nothing to pull.
		{"Bootstrap: docker\n", false, false},
	}
	for _, tt := range tests {
		b := parseBootstrap([]byte(tt.def))
		if b.resolvable() != tt.resolvable || b.remote() != tt.remote {
			t.Errorf("%q: resolvable=%v remote=%v, want %v/%v",
				tt.def, b.resolvable(), b.remote(), tt.resolvable, tt.remote)
		}
	}
}

func TestPinned(t *testing.T) {
	digest := "sha256:" + strings.Repeat("a", 64)
	tests := []struct{ from, want string }{
		{"ubuntu:24.04", "ubuntu@" + digest},
		{"docker.io/library/ubuntu:24.04", "docker.io/library/ubuntu@" + digest},
		{"ubuntu", "ubuntu@" + digest},
		// A port is not a tag.
		{"reg.example.com:5000/img", "reg.example.com:5000/img@" + digest},
		{"reg.example.com:5000/img:1.2", "reg.example.com:5000/img@" + digest},
		// Already pinned: nothing to do.
		{"ubuntu@" + digest, "ubuntu@" + digest},
	}
	for _, tt := range tests {
		if got := (bootstrap{Agent: "docker", From: tt.from}).pinned(digest); got != tt.want {
			t.Errorf("pinned(%q) = %q, want %q", tt.from, got, tt.want)
		}
	}
}

// The definition handed to Apptainer names the resolved digest, so an upstream
// retagged mid-build cannot leave the record describing bytes it does not hold.
func TestPinDefinition(t *testing.T) {
	digest := "sha256:" + strings.Repeat("a", 64)
	def := []byte("#DESC:x\nBootstrap: docker\nFrom: ubuntu:24.04\n\n%post\n    echo hi\n")

	got, changed := pinDefinition(def, digest)
	if !changed {
		t.Fatal("the definition was not pinned")
	}
	text := string(got)
	if !strings.Contains(text, "From: ubuntu@"+digest) {
		t.Errorf("From: was not pinned:\n%s", text)
	}
	if !strings.Contains(text, "%post\n    echo hi\n") {
		t.Errorf("the body was disturbed:\n%s", text)
	}

	t.Run("nothing to pin", func(t *testing.T) {
		for _, tt := range []struct {
			name, def, digest string
		}{
			{"no digest", "Bootstrap: docker\nFrom: ubuntu:24.04\n", ""},
			{"unresolved", "Bootstrap: docker\nFrom: ubuntu:24.04\n", "unrecorded"},
			{"not a registry", "Bootstrap: scratch\n", digest},
			{"already pinned", "Bootstrap: docker\nFrom: ubuntu@" + digest + "\n", digest},
		} {
			if out, changed := pinDefinition([]byte(tt.def), tt.digest); changed {
				t.Errorf("%s: rewrote the definition:\n%s", tt.name, out)
			}
		}
	})
}

// docker and library serve different images under one name, so the Bootstrap:
// directive is recorded. It stays out of the key, where the digest already
// distinguishes them.
func TestUpstreamRecordsTheBootstrap(t *testing.T) {
	tests := []struct {
		def  string
		want string
	}{
		{"Bootstrap: docker\nFrom: alpine:3.20\n", "docker://alpine:3.20"},
		{"Bootstrap: library\nFrom: alpine:3.20\n", "library://alpine:3.20"},
		{"Bootstrap: oras\nFrom: reg.io/img:1\n", "oras://reg.io/img:1"},
	}
	for _, tt := range tests {
		b := parseBootstrap([]byte(tt.def))
		from := meta.From{Bootstrap: b.Agent, Ref: b.From}
		if got := from.URI(); got != tt.want {
			t.Errorf("URI() = %q, want %q", got, tt.want)
		}
	}
}
