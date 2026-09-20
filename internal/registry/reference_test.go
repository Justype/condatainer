package registry

import (
	"slices"
	"strings"
	"testing"

	"github.com/Justype/condatainer/catalog"
	"github.com/Justype/condatainer/internal/artifact/meta"
)

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

func projectManifest(name, sha string) meta.Manifest {
	return meta.Manifest{
		Name: name,
		Type: catalog.TypeApp,
		Keys: meta.Keys{Identity: meta.KeyRef{Scheme: "script-identity-v1", SHA256: sha}},
	}
}

func TestProjectTags(t *testing.T) {
	sha := "a31f902c12ab" + strings.Repeat("0", 52)

	// A selection gets both: the qualified tag is the retention anchor every
	// artifact needs, the plain one is the human handle that says which identity
	// the project uses now.
	got, err := ProjectTags(projectManifest("star/2.7.11b", sha), true)
	if err != nil {
		t.Fatalf("ProjectTags: %v", err)
	}
	want := []string{"star--2.7.11b__a31f902c12ab", "star--2.7.11b"}
	if !slices.Equal(got, want) {
		t.Errorf("tags = %v, want %v", got, want)
	}

	// A closure-only artifact gets only the anchor: nothing points at it by
	// name, and a plain tag would claim it is what the project uses.
	got, err = ProjectTags(projectManifest("grch38/genome/gencode49", sha), false)
	if err != nil {
		t.Fatalf("ProjectTags: %v", err)
	}
	if !slices.Equal(got, []string{"grch38--genome--gencode49__a31f902c12ab"}) {
		t.Errorf("closure tags = %v", got)
	}
}

func TestProjectTagsRefuseWhatCannotBeATag(t *testing.T) {
	sha := strings.Repeat("a", 64)
	// Refused, never truncated: a truncated tag is a different artifact's
	// address, and the name is what a puller reads back out of it.
	long := strings.Repeat("segment/", 20) + "1.0"
	if _, err := ProjectTags(projectManifest(long, sha), true); err == nil {
		t.Error("a name too long for an OCI tag must be refused")
	}
	if _, err := ProjectTags(projectManifest("star/2.7.11b", ""), true); err == nil {
		t.Error("an artifact with no identity has no retention anchor and must be refused")
	}
	if _, err := ProjectTags(projectManifest("", sha), true); err == nil {
		t.Error("an artifact with no name must be refused")
	}
}

func TestParseProjectTag(t *testing.T) {
	for tag, want := range map[string]string{
		"star--2.7.11b__a31f902c12ab":             "star/2.7.11b",
		"star--2.7.11b":                           "star/2.7.11b",
		"grch38--genome--gencode49__a31f902c12ab": "grch38/genome/gencode49",
		// A single-component name has no `--` to be recognised by. The `__` is
		// the proof for the qualified tag; the plain one is accepted because a
		// frozen environment is the only artifact named without a version.
		"env__a31f902c12ab": "env",
		"env":               "env",
	} {
		got, _, ok := ParseProjectTag(tag)
		if !ok || got != want {
			t.Errorf("ParseProjectTag(%q) = %q, %v; want %q", tag, got, ok, want)
		}
	}
	// A catalog tag is a single version segment and never contains `--`, so the
	// two namespaces cannot be confused into inventing a name.
	for _, notAProjectTag := range []string{"gencode49", "latest", "20260721", "4.4.4", "", "__abc", "star--2.7.11b__zz"} {
		if _, _, ok := ParseProjectTag(notAProjectTag); ok {
			t.Errorf("ParseProjectTag(%q) claimed to be a project tag", notAProjectTag)
		}
	}
}

// Round trip: what a push writes is what a pull reads back as the name.
func TestProjectTagRoundTrip(t *testing.T) {
	sha := strings.Repeat("c", 64)
	for _, name := range []string{"star/2.7.11b", "grch38/genome/gencode49", "ubuntu24/r/4.4.4"} {
		tags, err := ProjectTags(projectManifest(name, sha), true)
		if err != nil {
			t.Fatalf("%s: %v", name, err)
		}
		for _, tag := range tags {
			got, _, ok := ParseProjectTag(tag)
			if !ok || got != name {
				t.Errorf("%s -> %q -> %q, %v", name, tag, got, ok)
			}
		}
	}
}

// A frozen environment publishes under a name with no version, which is the one
// shape the plain `--` rule cannot recognise. What a project push writes for it
// must still read back as `env`, or a pull would install it under a name the
// project never used.
func TestProjectTagRoundTripSnapshot(t *testing.T) {
	m := meta.Manifest{Name: meta.EnvName, Type: catalog.TypeEnv, BuildType: meta.BuildTypeSnapshot}
	m.Keys.Identity = meta.KeyRef{Scheme: "payload-tree-v1", SHA256: strings.Repeat("d", 64)}
	m.Keys.Equiv = m.Keys.Identity

	tags, err := ProjectTags(m, true)
	if err != nil {
		t.Fatalf("ProjectTags: %v", err)
	}
	if !slices.Equal(tags, []string{"env__dddddddddddd", "env"}) {
		t.Fatalf("tags = %v", tags)
	}
	for _, tag := range tags {
		got, _, ok := ParseProjectTag(tag)
		if !ok || got != meta.EnvName {
			t.Errorf("%q -> %q, %v; want %q", tag, got, ok, meta.EnvName)
		}
	}
}
