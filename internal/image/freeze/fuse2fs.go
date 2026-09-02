package freeze

import (
	"bytes"
	"context"
	"errors"
	"fmt"
	"os/exec"
	"path/filepath"
	"strings"
)

// ErrNoFuse2fs reports that no fuse2fs could be found to read the overlay with.
var ErrNoFuse2fs = errors.New("no fuse2fs available to read the overlay")

// ErrNoSquashfuse reports that no squashfuse could be found to read an artifact.
var ErrNoSquashfuse = errors.New("no squashfuse available to read the artifact")

// findSquashfuse locates the squashfuse that mounts a frozen artifact, the same
// way findFuse2fs locates its counterpart: Apptainer bundles one and uses it to
// mount .sqf overlays, so it is present wherever one can be mounted at all.
func findSquashfuse(ctx context.Context, apptainerBin string) (string, error) {
	for _, name := range []string{"squashfuse_ll", "squashfuse"} {
		if p := findHelper(ctx, apptainerBin, name); p != "" {
			return p, nil
		}
	}
	return "", fmt.Errorf("%w: looked in apptainer's libexec and on PATH", ErrNoSquashfuse)
}

// findFuse2fs locates the fuse2fs that will read the image.
//
// Apptainer bundles one and uses it to mount ext3 overlays unprivileged, so it is
// present wherever a writable overlay works — but installs differ on where.
// `apptainer buildcfg` is authoritative; the rest are layouts seen in the wild.
//
// Only the path is resolved. Whether it can mount is left to the pack, which runs
// it through Apptainer: mounting it ourselves fails inside a nested container,
// where a setuid helper cannot gain privilege.
func findFuse2fs(ctx context.Context, apptainerBin string) (string, error) {
	if p := findHelper(ctx, apptainerBin, "fuse2fs"); p != "" {
		return p, nil
	}
	return "", fmt.Errorf("%w: looked in apptainer's libexec and on PATH", ErrNoFuse2fs)
}

// findHelper locates one of Apptainer's bundled FUSE helpers. buildcfg is
// authoritative; the rest are layouts seen in the wild, PATH last. Only the path
// is resolved — whether it can mount is left to the caller, which runs it through
// Apptainer.
func findHelper(ctx context.Context, apptainerBin, name string) string {
	var candidates []string

	if libexec := buildcfgLibexec(ctx, apptainerBin); libexec != "" {
		candidates = append(candidates,
			filepath.Join(libexec, "apptainer", "bin", name),
			filepath.Join(libexec, "singularity", "bin", name))
	}
	if apptainerBin != "" {
		candidates = append(candidates,
			filepath.Join(filepath.Dir(apptainerBin), "..", "libexec", "apptainer", "bin", name))
	}
	candidates = append(candidates,
		filepath.Join("/usr/libexec/apptainer/bin", name),
		filepath.Join("/usr/local/libexec/apptainer/bin", name),
		filepath.Join("/opt/apptainer/libexec/apptainer/bin", name),
		filepath.Join("/usr/libexec/singularity/bin", name),
	)
	if p, err := exec.LookPath(name); err == nil {
		candidates = append(candidates, p)
	}
	candidates = append(candidates,
		filepath.Join("/usr/bin", name),
		filepath.Join("/usr/sbin", name),
		filepath.Join("/sbin", name))

	seen := map[string]bool{}
	for _, c := range candidates {
		abs, err := filepath.Abs(c)
		if err != nil || seen[abs] {
			continue
		}
		seen[abs] = true
		if isExecutable(abs) {
			return abs
		}
	}
	return ""
}

// buildcfgLibexec asks apptainer where it installed its helpers.
func buildcfgLibexec(ctx context.Context, apptainerBin string) string {
	if apptainerBin == "" {
		apptainerBin = "apptainer"
	}
	cmd := exec.CommandContext(ctx, apptainerBin, "buildcfg")
	var out bytes.Buffer
	cmd.Stdout = &out
	if err := cmd.Run(); err != nil {
		return ""
	}
	for _, line := range strings.Split(out.String(), "\n") {
		if v, ok := strings.CutPrefix(strings.TrimSpace(line), "LIBEXECDIR="); ok {
			return v
		}
	}
	return ""
}
