package registry

import "testing"

func TestNormalize(t *testing.T) {
	tests := []struct{ in, want string }{
		// A bare name is a Docker Hub official image.
		{"ubuntu:24.04", "docker.io/library/ubuntu:24.04"},
		{"ubuntu", "docker.io/library/ubuntu"},
		// One component with no dot or colon is a Hub namespace, not a host.
		{"myorg/tool:1.0", "docker.io/myorg/tool:1.0"},
		// A dot, a colon, or localhost makes it a registry host.
		{"registry.example.com/tool:1.0", "registry.example.com/tool:1.0"},
		{"localhost/tool:1.0", "localhost/tool:1.0"},
		{"localhost:5000/tool:1.0", "localhost:5000/tool:1.0"},
		{"quay.io/biocontainers/star:2.7.11b", "quay.io/biocontainers/star:2.7.11b"},
		{"  ubuntu:24.04  ", "docker.io/library/ubuntu:24.04"},
	}
	for _, tt := range tests {
		got, err := Normalize(tt.in)
		if err != nil {
			t.Errorf("Normalize(%q): %v", tt.in, err)
			continue
		}
		if got != tt.want {
			t.Errorf("Normalize(%q) = %q, want %q", tt.in, got, tt.want)
		}
	}
}

func TestNormalizeRejects(t *testing.T) {
	for _, in := range []string{"", "   ", "docker://ubuntu:24.04"} {
		if got, err := Normalize(in); err == nil {
			t.Errorf("Normalize(%q) = %q, want an error", in, got)
		}
	}
}
