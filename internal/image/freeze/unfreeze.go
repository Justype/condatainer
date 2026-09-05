package freeze

import (
	"bufio"
	"bytes"
	"context"
	"errors"
	"fmt"
	"os"
	"os/exec"
	"path/filepath"
	"strconv"
	"strings"
	"time"

	"github.com/Justype/condatainer/internal/artifact/meta"
	"github.com/Justype/condatainer/internal/image/ext3"
	"github.com/Justype/condatainer/internal/image/tool"
	"github.com/Justype/condatainer/internal/logging"
	execpkg "github.com/Justype/condatainer/internal/runtime/exec"
	"github.com/Justype/condatainer/internal/utils"
)

// Errors a caller distinguishes.
var (
	// ErrNotFrozen reports an artifact that is not a frozen environment.
	ErrNotFrozen = errors.New("not a frozen environment")
	// ErrTooSmall reports a requested size the payload does not fit in.
	ErrTooSmall = errors.New("requested size is smaller than the payload")
)

// headroomMB is added to the payload when no size is given. An environment being
// unfrozen is an environment about to be installed into.
const headroomMB = 5 * 1024

// UnfreezeOptions is one unfreeze.
type UnfreezeOptions struct {
	// Artifact is the frozen .sqf to open.
	Artifact string
	// Target is the writable .img to create.
	Target string
	// SizeMB is the image size; zero derives it from the payload plus headroom.
	SizeMB int
	// UID and GID own the result, fixed up after the build — see buildImage.
	UID, GID int
	// Sparse allocates on write rather than up front.
	Sparse bool
	// StageDir is where the skeleton is assembled. It holds work/ and an empty
	// mountpoint — a few directories, never the payload.
	StageDir string
}

// UnfreezeResult is what an unfreeze produced.
type UnfreezeResult struct {
	Path string
	// From is the artifact's name, recorded in the sidecar so that three months
	// later a user can tell what they are developing against.
	From string
	// PayloadMB is what the archive holds; SizeMB is the image built for it.
	PayloadMB int
	SizeMB    int
	// Entries is how many things the payload contains, which decides the inode
	// table.
	Entries int
}

// Unfreeze turns a frozen environment back into a writable overlay.
//
// The image is built straight from the mounted artifact — mke2fs -d takes a
// directory tree, and a mounted SquashFS is one — so the payload is never staged.
// That also keeps it correct: mke2fs records whiteout device nodes an unprivileged
// unsquashfs cannot create. The artifact is packed -all-root and the mount cannot
// be told to report a different owner (see buildImage), so ownership is fixed
// afterward with a plain chown pass instead.
func Unfreeze(ctx context.Context, opts UnfreezeOptions) (UnfreezeResult, error) {
	log := logging.FromContext(ctx)

	if err := tool.CheckDependencies([]string{"unsquashfs", "debugfs"}); err != nil {
		return UnfreezeResult{}, err
	}

	m, err := meta.ReadManifest(opts.Artifact)
	if err != nil {
		return UnfreezeResult{}, fmt.Errorf("%s: %w", opts.Artifact, err)
	}
	if m.BuildType != meta.BuildTypeSnapshot {
		return UnfreezeResult{}, fmt.Errorf("%w: %s is a %s build; unfreeze is the inverse of freeze and takes what freeze produced",
			ErrNotFrozen, opts.Artifact, m.BuildType)
	}

	artifact, err := filepath.Abs(opts.Artifact)
	if err != nil {
		return UnfreezeResult{}, err
	}
	target, err := filepath.Abs(opts.Target)
	if err != nil {
		return UnfreezeResult{}, err
	}

	// The freeze recorded what the payload needs, having walked the tree to
	// produce the archive. Reading the listing back would recover the same two
	// numbers at the cost of decompressing the metadata of every entry.
	payloadMB, entries := 0, 0
	if m.Snapshot != nil {
		payloadMB, entries = m.Snapshot.PayloadMB, m.Snapshot.Entries
	}
	if payloadMB == 0 {
		// An artifact frozen before those were recorded.
		payloadMB, entries, err = archiveSize(ctx, artifact)
		if err != nil {
			return UnfreezeResult{}, err
		}
	}
	size := opts.SizeMB
	if size == 0 {
		size = payloadMB + headroomMB
	}
	// mke2fs -d accepts an undersized target and only discovers the problem
	// while writing, so the refusal has to be ours and has to come first.
	if size <= payloadMB {
		return UnfreezeResult{}, fmt.Errorf("%w: %s holds %d MB and %d MB was asked for",
			ErrTooSmall, opts.Artifact, payloadMB, size)
	}

	stageRoot := opts.StageDir
	if stageRoot == "" {
		stageRoot = utils.GetTmpDir()
	}
	if err := os.MkdirAll(stageRoot, 0o755); err != nil {
		return UnfreezeResult{}, fmt.Errorf("stage skeleton: %w", err)
	}
	stage, err := os.MkdirTemp(stageRoot, "cnt-unfreeze-")
	if err != nil {
		return UnfreezeResult{}, fmt.Errorf("stage skeleton: %w", err)
	}
	defer os.RemoveAll(stage)

	// The skeleton is what an overlay image holds beside its payload: an empty
	// upper/ for the artifact to be mounted onto, and the work/ OverlayFS needs.
	for _, dir := range []string{"upper", "work/work"} {
		if err := os.MkdirAll(filepath.Join(stage, dir), 0o755); err != nil {
			return UnfreezeResult{}, err
		}
	}

	log.Info("rebuilding writable overlay", "artifact", opts.Artifact, "target", opts.Target,
		"payload_mb", payloadMB, "size_mb", size, "entries", entries)
	if err := buildImage(ctx, opts, artifact, target, stage, size, entries); err != nil {
		return UnfreezeResult{}, err
	}
	// Read before the metadata is dropped: an unfrozen overlay has nowhere else
	// to keep its environment, and the artifact's copy is about to be deleted.
	var env []meta.EnvVar
	if rt, rtErr := meta.ReadRuntime(artifact); rtErr == nil {
		env = rt.Env
	}
	if err := dropMetadata(ctx, target); err != nil {
		return UnfreezeResult{}, err
	}
	if err := writeSidecar(target, m, env); err != nil {
		return UnfreezeResult{}, err
	}

	return UnfreezeResult{
		Path: target, From: m.Name,
		PayloadMB: payloadMB, SizeMB: size, Entries: entries,
	}, nil
}

// archiveSize measures the payload as a filesystem will hold it, from the archive
// listing rather than the mount: the size is needed before the image exists.
// Sizes round up to a block, since a conda prefix is mostly small files and
// summing apparent bytes understates it by orders of magnitude.
func archiveSize(ctx context.Context, sqf string) (sizeMB, entries int, err error) {
	out, err := exec.CommandContext(ctx, "unsquashfs", "-ll", sqf).Output()
	if err != nil {
		return 0, 0, &tool.Error{Op: "list", Path: sqf, Tool: "unsquashfs", BaseErr: err}
	}

	const block = 4096
	var bytesUsed int64
	scanner := bufio.NewScanner(bytes.NewReader(out))
	scanner.Buffer(make([]byte, 0, 64*1024), 4*1024*1024)
	for scanner.Scan() {
		fields := strings.Fields(scanner.Text())
		if len(fields) < 6 || !strings.Contains(scanner.Text(), "squashfs-root") {
			continue
		}
		entries++
		if fields[0][0] != '-' { // only a regular file occupies data blocks
			bytesUsed += block
			continue
		}
		size, convErr := strconv.ParseInt(fields[2], 10, 64)
		if convErr != nil {
			continue
		}
		bytesUsed += (size + block - 1) / block * block
	}
	if err := scanner.Err(); err != nil {
		return 0, 0, err
	}
	return int(bytesUsed/(1024*1024)) + 1, entries, nil
}

// buildImage creates the ext3 image with the artifact mounted as its payload,
// then fixes its ownership.
//
// The inode count comes from the payload rather than mke2fs's one-per-16KB
// ratio, which a conda prefix exhausts while plenty of space remains. The
// artifact is packed -all-root, and mountedRun's namespace maps exactly one
// uid (0) correctly — any -o uid override besides that comes back as the
// overflow uid instead of the one asked for — so the payload lands owned by
// root, and ChownRecursively repairs it afterward on the plain host, where no
// such mapping limit applies.
func buildImage(ctx context.Context, opts UnfreezeOptions, artifact, target, stage string, sizeMB, entries int) error {
	squashfuse, err := findSquashfuse()
	if err != nil {
		return err
	}
	upper := filepath.Join(stage, "upper")

	// -m matches what overlay create makes: the mke2fs default reserves blocks
	// for root, which nothing in an unprivileged overlay can use.
	args := []string{"-q", "-t", "ext3", "-d", stage,
		"-m", strconv.Itoa(ext3.ProfileSmall.ReservedPerc), "-F", target}
	if entries > 0 {
		args = append(args, "-N", strconv.Itoa(max(entries*2, 1024)))
	}
	args = append(args, fmt.Sprintf("%dM", sizeMB))

	script := "set -e\n"
	if !opts.Sparse {
		script += fmt.Sprintf("dd if=/dev/zero of=%s bs=1M count=%d status=none\n", target, sizeMB)
	}
	script += "mke2fs " + strings.Join(args, " ") + "\n"

	if err := mountedRun(ctx, squashfuse, []string{artifact}, upper, script, execpkg.IOFromContext(ctx)); err != nil {
		os.Remove(target)
		return fmt.Errorf("build %s: %w", target, err)
	}
	if err := ext3.ChownRecursively(ctx, target, opts.UID, opts.GID, "/"); err != nil {
		os.Remove(target)
		return fmt.Errorf("fix ownership on %s: %w", target, err)
	}
	return nil
}

// dropMetadata removes the embedded .cnt from the rebuilt overlay.
//
// A writable overlay carries no embedded metadata by design: it is mutable
// working state, so what it needs travels beside it in the sidecar instead. The
// artifact was mounted read-only, so this is deleted from the image afterwards
// rather than skipped on the way in.
func dropMetadata(ctx context.Context, img string) error {
	dir := "/upper/" + meta.DirName
	names, err := metadataNames(ctx, img, dir)
	if err != nil || len(names) == 0 {
		return err
	}
	var script bytes.Buffer
	for _, name := range names {
		fmt.Fprintf(&script, "rm %s/%s\n", dir, name)
	}
	fmt.Fprintf(&script, "rmdir %s\nquit\n", dir)
	return runDebugfs(ctx, img, script.String(), "remove embedded metadata")
}

// metadataNames lists what the embedded metadata directory holds.
func metadataNames(ctx context.Context, img, dir string) ([]string, error) {
	out, err := exec.CommandContext(ctx, "debugfs", "-R", "ls -p "+dir, img).Output()
	if err != nil {
		return nil, nil // no metadata directory is not a problem
	}
	var names []string
	for _, line := range strings.Split(string(out), "\n") {
		fields := strings.Split(strings.ReplaceAll(strings.TrimSpace(line), " ", ""), "/")
		if len(fields) < 7 || fields[5] == "" || fields[5] == "." || fields[5] == ".." {
			continue
		}
		names = append(names, fields[5])
	}
	return names, nil
}

func runDebugfs(ctx context.Context, img, script, op string) error {
	cmd := exec.CommandContext(ctx, "debugfs", "-w", img)
	cmd.Stdin = strings.NewReader(script)
	var out bytes.Buffer
	cmd.Stdout = &out
	cmd.Stderr = &out
	if err := cmd.Run(); err != nil {
		return &tool.Error{Op: op, Path: img, Tool: "debugfs", Output: out.String(), BaseErr: err}
	}
	return nil
}

// writeSidecar records what this overlay was unfrozen from.
//
// A writable overlay carries no embedded metadata, so the sidecar is the only
// place this can live. It is written as a comment: the sidecar's parser skips
// those, so provenance cannot be mistaken for an environment variable, and a
// person reading the file three months later still sees the answer.
func writeSidecar(img string, m meta.Manifest, env []meta.EnvVar) error {
	path := img + ".env"
	if _, err := os.Stat(path); err == nil {
		return nil // never overwrite an environment the user maintains
	}
	var b strings.Builder
	fmt.Fprintf(&b, "# unfrozen from %s on %s\n", m.Name, time.Now().UTC().Format("2006-01-02"))
	b.WriteString("# it is a new development line: this overlay has no identity, is not lockable,\n")
	b.WriteString("# and is not publishable until it is frozen again.\n")

	// The artifact's environment comes back out. A frozen artifact keeps its
	// variables inside, where an immutable image's metadata belongs; a writable
	// overlay reads them from here and from nowhere else, so a round trip that
	// did not write them would quietly lose them. {prefix} stays a token, as it
	// is on both sides.
	for _, v := range env {
		if v.Note != "" {
			fmt.Fprintf(&b, "%s=%s  ## %s\n", v.Key, v.Value, v.Note)
			continue
		}
		fmt.Fprintf(&b, "%s=%s\n", v.Key, v.Value)
	}
	return os.WriteFile(path, []byte(b.String()), 0o644)
}
