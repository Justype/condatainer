package registry

import "testing"

// Normalize and TrimBaseScheme both "clean up a reference" and do opposite things:
// one supplies a registry host that was left out, the other removes a scheme that
// was typed. Confusing them silently corrupts a base, so pin a case each way.
func TestNormalizeAndTrimBaseSchemeAreNotInterchangeable(t *testing.T) {
	t.Run("Normalize supplies a host, TrimBaseScheme does not", func(t *testing.T) {
		const short = "ubuntu:24.04"
		got, err := Normalize(short)
		if err != nil {
			t.Fatalf("Normalize: %v", err)
		}
		if got != "docker.io/library/ubuntu:24.04" {
			t.Errorf("Normalize(%q) = %q", short, got)
		}
		if got := TrimBaseScheme(short); got != short {
			t.Errorf("TrimBaseScheme(%q) = %q, want it unchanged", short, got)
		}
	})

	t.Run("TrimBaseScheme removes a scheme, Normalize rejects one", func(t *testing.T) {
		const scheme = "oci://ghcr.io/lab/cnt"
		if got := TrimBaseScheme(scheme); got != "ghcr.io/lab/cnt" {
			t.Errorf("TrimBaseScheme(%q) = %q", scheme, got)
		}
		if _, err := Normalize(scheme); err == nil {
			t.Errorf("Normalize(%q) accepted a URI", scheme)
		}
	})
}

func TestTrimBaseScheme(t *testing.T) {
	tests := []struct{ in, want string }{
		{"ghcr.io/lab/cnt", "ghcr.io/lab/cnt"},
		{"oci://ghcr.io/lab/cnt", "ghcr.io/lab/cnt"},
		{"oras://ghcr.io/lab/cnt", "ghcr.io/lab/cnt"},
		{"  ghcr.io/lab/cnt  ", "ghcr.io/lab/cnt"},
		{"ghcr.io/lab/cnt/", "ghcr.io/lab/cnt"},
		{"oci://ghcr.io/lab/cnt///", "ghcr.io/lab/cnt"},
		{"localhost:5000/cnt", "localhost:5000/cnt"},
		{"", ""},
	}
	for _, tt := range tests {
		if got := TrimBaseScheme(tt.in); got != tt.want {
			t.Errorf("TrimBaseScheme(%q) = %q, want %q", tt.in, got, tt.want)
		}
	}
}

// A digest joins with @ and a tag with :, so one function serves both selectors.
func TestFullRef(t *testing.T) {
	const digest = "sha256:0123456789abcdef0123456789abcdef0123456789abcdef0123456789abcdef"
	tests := []struct{ base, repo, tag, want string }{
		{"ghcr.io/lab/cnt", "cellranger", "9.0.1", "ghcr.io/lab/cnt/cellranger:9.0.1"},
		{"oci://ghcr.io/lab/cnt/", "ubuntu24/base", "latest", "ghcr.io/lab/cnt/ubuntu24/base:latest"},
		{"ghcr.io/lab/cnt", "ubuntu24/base", "20260721", "ghcr.io/lab/cnt/ubuntu24/base:20260721"},
		{"ghcr.io/lab/cnt", "cellranger", digest, "ghcr.io/lab/cnt/cellranger@" + digest},
	}
	for _, tt := range tests {
		if got := FullRef(tt.base, tt.repo, tt.tag); got != tt.want {
			t.Errorf("FullRef(%q, %q, %q) = %q, want %q", tt.base, tt.repo, tt.tag, got, tt.want)
		}
	}
}

// Plain HTTP is a downgrade, so it must reach loopback and stop there.
func TestIsLoopback(t *testing.T) {
	for _, host := range []string{"localhost", "localhost:5000", "127.0.0.1", "127.0.0.1:5000", "::1", "[::1]:5000"} {
		if !isLoopback(host) {
			t.Errorf("isLoopback(%q) = false, want true", host)
		}
	}
	for _, host := range []string{
		"ghcr.io", "registry.example.test:5000", "localhost.example.test", "127.0.0.1.example.test",
	} {
		if isLoopback(host) {
			t.Errorf("isLoopback(%q) = true — TLS must not be downgraded for it", host)
		}
	}
}
