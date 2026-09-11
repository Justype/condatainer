package build

import (
	"errors"
	"os"
	"path/filepath"
	"testing"

	"github.com/Justype/condatainer/internal/artifact/meta"
	"github.com/Justype/condatainer/internal/image/tool"
)

func TestAllowedForeignBootstrap(t *testing.T) {
	for _, agent := range []string{"docker", "oras", "library"} {
		if err := allowedForeignBootstrap(agent); err != nil {
			t.Errorf("%s: got %v, want nil", agent, err)
		}
	}
	for _, agent := range []string{"shub", "yum", "zypper", "debootstrap", "localimage", "scratch", ""} {
		if err := allowedForeignBootstrap(agent); err == nil {
			t.Errorf("%q: want an error, got nil", agent)
		}
	}
}

func TestDigestFromLabels(t *testing.T) {
	const want = "sha256:6baf435"
	got := digestFromLabels([]byte(`{"org.opencontainers.image.base.digest":"` + want + `"}`))
	if got != want {
		t.Errorf("digest = %q, want %q", got, want)
	}

	// Older Apptainer/Singularity builds carry only org.label-schema.* labels —
	// no base-digest at all — which must fall back, not error.
	if got := digestFromLabels([]byte(`{"org.label-schema.build-arch":"amd64"}`)); got != meta.Unrecorded {
		t.Errorf("missing label: got %q, want %q", got, meta.Unrecorded)
	}
	if got := digestFromLabels([]byte("not json")); got != meta.Unrecorded {
		t.Errorf("unparsable labels: got %q, want %q", got, meta.Unrecorded)
	}
}

func writeSandboxRoot(t *testing.T, singularity, labels string) string {
	t.Helper()
	root := t.TempDir()
	dir := filepath.Join(root, ".singularity.d")
	if err := os.MkdirAll(dir, 0o755); err != nil {
		t.Fatal(err)
	}
	if singularity != "" {
		if err := os.WriteFile(filepath.Join(dir, "Singularity"), []byte(singularity), 0o644); err != nil {
			t.Fatal(err)
		}
	}
	if labels != "" {
		if err := os.WriteFile(filepath.Join(dir, "labels.json"), []byte(labels), 0o644); err != nil {
			t.Fatal(err)
		}
	}
	return root
}

func TestReadForeignBootstrapSandbox(t *testing.T) {
	root := writeSandboxRoot(t, "bootstrap: docker\nfrom: alpine:3.19\n",
		`{"org.opencontainers.image.base.digest":"sha256:abc"}`)

	boot, digest, err := readForeignBootstrap(root, true)
	if err != nil {
		t.Fatalf("readForeignBootstrap: %v", err)
	}
	if boot.Agent != "docker" || boot.From != "alpine:3.19" {
		t.Errorf("boot = %+v, want {docker alpine:3.19}", boot)
	}
	if digest != "sha256:abc" {
		t.Errorf("digest = %q, want sha256:abc", digest)
	}
}

// A sandbox from an older Apptainer/Singularity build carries no
// labels.json at all, which readForeignBootstrap must treat as unrecorded
// rather than a failure — the bootstrap record is still there.
func TestReadForeignBootstrapSandboxNoLabels(t *testing.T) {
	root := writeSandboxRoot(t, "bootstrap: library\nfrom: alpine:3.11.5\n", "")

	boot, digest, err := readForeignBootstrap(root, true)
	if err != nil {
		t.Fatalf("readForeignBootstrap: %v", err)
	}
	if boot.Agent != "library" || boot.From != "alpine:3.11.5" {
		t.Errorf("boot = %+v, want {library alpine:3.11.5}", boot)
	}
	if digest != meta.Unrecorded {
		t.Errorf("digest = %q, want %q", digest, meta.Unrecorded)
	}
}

func TestReadForeignBootstrapSandboxNoRecord(t *testing.T) {
	root := writeSandboxRoot(t, "", "")
	if _, _, err := readForeignBootstrap(root, true); err == nil {
		t.Error("no .singularity.d/Singularity: want an error")
	}
}

func TestReadSandboxFileMissing(t *testing.T) {
	_, err := readSandboxFile(t.TempDir(), "nope")
	if !errors.Is(err, tool.ErrFileNotFound) {
		t.Errorf("err = %v, want ErrFileNotFound", err)
	}
}
