package freeze

import (
	"bytes"
	"context"
	"errors"
	"fmt"
	"os"
	"os/exec"
	"path"
	"path/filepath"
	"sort"
	"strings"
	"time"

	"github.com/Justype/condatainer/catalog"
	"github.com/Justype/condatainer/internal/artifact/meta"
	"github.com/Justype/condatainer/internal/logging"
	"github.com/Justype/condatainer/internal/runtime/container"
	"github.com/Justype/condatainer/internal/toolpath"
	"github.com/Justype/condatainer/internal/utils"
)

// Errors a caller distinguishes.
var (
	// ErrEmptyOverlay reports an overlay nobody has written to.
	ErrEmptyOverlay = errors.New("overlay has no payload to freeze")
)

// Options is one freeze.
type Options struct {
	// Image is the writable overlay to freeze. It is read, never written.
	Image string
	// Target is the artifact to write.
	Target string
	// Description is optional prose for `info`.
	Description string
	// Base is the base .sqf an opaque translation resolves its directories
	// against. Empty skips that translation.
	Base string
	// CompressArgs, BlockSize and Processors tune the pack.
	CompressArgs string
	BlockSize    string
	Processors   int
	// UseTmp copies the payload to fast scratch and packs it from there instead
	// of reading the image through a FUSE mount. It is faster and needs room for
	// the payload a second time, so it is the caller's choice, never inferred.
	UseTmp bool
}

// Result is what a freeze produced.
type Result struct {
	// Path is the artifact.
	Path string
	// Identity is what the artifact was keyed as, hashed from the packed payload.
	// It is also embedded in Manifest, which is why the pack runs in two phases.
	Identity meta.KeyRef
	// Manifest is what was embedded.
	Manifest meta.Manifest
	// Translation is what the deletions became.
	Translation Translation
}

// Freeze packs a writable overlay into an immutable artifact.
//
// It does not rebuild, re-solve or relocate: the payload keeps the paths it has,
// because a conda prefix is full of absolute ones.
func Freeze(ctx context.Context, opts Options) (Result, error) {
	log := logging.FromContext(ctx)

	if _, err := os.Stat(opts.Image); err != nil {
		return Result{}, err
	}

	// The scan reads every inode in the image and resolves opaque directories
	// against the base, which on a real environment is seconds of silence before
	// the pack's own progress starts.
	log.Info("scanning overlay", "image", opts.Image)
	entries, err := Walk(ctx, opts.Image)
	if err != nil {
		return Result{}, err
	}
	if len(entries) == 0 {
		return Result{}, fmt.Errorf("%w: %s", ErrEmptyOverlay, opts.Image)
	}

	tr, err := Translate(ctx, opts.Image, entries, BaseDirLister(opts.Base))
	if err != nil {
		return Result{}, err
	}
	if tr.Deletions() > 0 {
		log.Info("translating deletions", "count", tr.Deletions(), "convention", string(tr.Convention))
	}

	roots := archiveSources(entries, tr)
	if len(roots) == 0 {
		return Result{}, fmt.Errorf("%w: %s", ErrEmptyOverlay, opts.Image)
	}

	// A writable overlay keeps its environment in a sidecar, because it is mutable
	// working state. A frozen one cannot: an installed image is immutable and its
	// metadata travels inside it, and the mount path ignores a sidecar beside a
	// .sqf for exactly that reason. So the variables move into the artifact here,
	// or they are lost at the freeze and lost for good.
	env, err := container.ReadEnvSidecar(opts.Image)
	if err != nil {
		return Result{}, fmt.Errorf("read %s.env: %w", opts.Image, err)
	}
	if len(env) > 0 {
		log.Info("carrying environment from the sidecar", "variables", len(env))
	}

	m, rt := describe(opts, tr, env, entries, buildTools(ctx))
	if err := meta.ValidateManifest(m); err != nil {
		return Result{}, err
	}
	if err := meta.ValidateRuntime(rt); err != nil {
		return Result{}, err
	}

	scratch := utils.GetTmpDir()
	if err := os.MkdirAll(scratch, 0o755); err != nil {
		return Result{}, fmt.Errorf("stage metadata: %w", err)
	}
	metaDir, err := os.MkdirTemp(scratch, "cnt-freeze-meta-")
	if err != nil {
		return Result{}, fmt.Errorf("stage metadata: %w", err)
	}
	defer os.RemoveAll(metaDir)

	packOpts := PackOptions{
		Image: opts.Image, Target: opts.Target,
		CompressArgs: opts.CompressArgs,
		BlockSize:    opts.BlockSize, Processors: opts.Processors,
	}
	packTr := tr
	if opts.UseTmp {
		payload, err := stagePayload(ctx, opts.Image, scratch, payloadSizeMB(entries))
		if err != nil {
			return Result{}, err
		}
		defer os.RemoveAll(payload)
		packOpts.StageDir = payload
		// rdump left every char 0:0 whiteout behind; they have to be injected.
		packTr = tr.ForCopy()
	}

	// The payload is packed first, on its own. Its identity is hashed from the
	// finished archive, and only then can the manifest carrying that identity be
	// written and appended — a manifest cannot describe the archive it is inside.
	if err := Pack(ctx, packOpts, entries, packTr); err != nil {
		return Result{}, err
	}
	log.Info("identifying the packed payload", "target", opts.Target)
	m.Keys.Identity, err = TreeIdentity(ctx, opts.Target)
	if err != nil {
		os.Remove(opts.Target)
		return Result{}, err
	}
	// Equivalence is the same value. A snapshot has no inputs to abstract away,
	// so "can this substitute" and "is this the same" are one question here.
	m.Keys.Equiv = m.Keys.Identity
	if err := meta.ValidateManifest(m); err != nil {
		os.Remove(opts.Target)
		return Result{}, err
	}

	staged := filepath.Join(metaDir, meta.DirName)
	if err := os.MkdirAll(staged, 0o755); err != nil {
		return Result{}, fmt.Errorf("stage metadata: %w", err)
	}
	if err := meta.StageManifest(staged, m); err != nil {
		return Result{}, err
	}
	if err := meta.StageRuntime(staged, rt); err != nil {
		return Result{}, err
	}
	packOpts.MetaDir = staged
	if err := AppendMeta(ctx, packOpts); err != nil {
		os.Remove(opts.Target)
		return Result{}, err
	}
	return Result{Path: opts.Target, Identity: m.Keys.Identity, Manifest: m, Translation: tr}, nil
}

// stagePayload copies the overlay's upper/ to fast scratch for the copy route,
// refusing before it starts if the space is not there.
func stagePayload(ctx context.Context, img, scratch string, payloadMB int) (string, error) {
	log := logging.FromContext(ctx)
	if err := checkStageSpace(scratch, payloadMB); err != nil {
		return "", err
	}
	dir, err := os.MkdirTemp(scratch, "cnt-freeze-payload-")
	if err != nil {
		return "", fmt.Errorf("stage payload: %w", err)
	}
	utils.WarnNetworkScratch(scratch, "freeze")
	log.Info("staging payload for the copy route", "dir", dir, "payload_mb", payloadMB)
	if err := dumpUpper(ctx, img, dir); err != nil {
		os.RemoveAll(dir)
		return "", err
	}
	return dir, nil
}

// describe builds the manifest and runtime a frozen environment embeds. It
// carries no keys: identity is the digest of the file this manifest is inside, so
// it is recorded by whoever refers to the artifact.
func describe(opts Options, tr Translation, env []meta.EnvVar, entries []Entry, tools meta.BuildTools) (meta.Manifest, meta.Runtime) {
	payloadMB := payloadSizeMB(entries)
	// Uname form, which is what recipes, Conda subdirs and users all say — and
	// what MountAllowed compares against, so GOARCH here would make every frozen
	// artifact look like it was built for another machine.
	platform := meta.NativePlatform()

	m := meta.Manifest{
		SchemaVersion: meta.SchemaVersion,
		Name:          meta.EnvName,
		Type:          catalog.TypeEnv,
		BuildType:     meta.BuildTypeSnapshot,
		Description:   opts.Description,
		Platform:      platform,
		Build: meta.Build{
			Tools:   tools,
			Created: time.Now().UTC(),
		},
		Snapshot: &meta.Snapshot{
			Convention: string(tr.Convention),
			Whiteouts:  tr.Deletions(),
			PayloadMB:  payloadMB,
			Entries:    len(entries),
		},
	}

	rt := meta.Runtime{
		SchemaVersion: meta.SchemaVersion,
		Name:          meta.EnvName,
		Type:          catalog.TypeEnv,
		Description:   opts.Description,
		Platform:      platform,
		// {prefix} is kept intact, as it is for every other artifact: the token
		// is substituted when the image is loaded.
		Env: env,
	}
	// The prefix is recorded only when the payload actually has one, since an
	// environment that is all apt packages and no conda prefix would otherwise
	// claim a path it does not provide.
	for _, r := range archiveSources(entries, tr) {
		if "/"+r == meta.EnvPrefix {
			rt.Prefix = meta.EnvPrefix
			break
		}
	}
	return m, rt
}

// payloadSizeMB is what the payload needs as a filesystem: every entry occupies
// at least a block, and a regular file occupies as many as its size rounds up to.
func payloadSizeMB(entries []Entry) int {
	const block = 4096
	var used int64
	for _, e := range entries {
		if e.Size > block {
			used += (e.Size + block - 1) / block * block
			continue
		}
		used += block
	}
	return int(used/(1024*1024)) + 1
}

// BaseDirLister lists directories inside the base image, all in one pass.
//
// "Hide everything the base has here" has no pseudo-file form, so an opaque
// directory becomes one whiteout per base entry — exact against this base and
// no other. base is a plain .sqf, read directly with unsquashfs — the same
// tool internal/image/squashfs already uses raw and unwrapped for every other
// .sqf read in this codebase — rather than mounted or run inside.
func BaseDirLister(base string) BaseLister {
	return func(ctx context.Context, dirs []string) (map[string][]string, error) {
		out := map[string][]string{}
		if base == "" || len(dirs) == 0 {
			return out, nil
		}

		bin, err := toolpath.Resolve("unsquashfs")
		if err != nil {
			return nil, err
		}
		args := append([]string{"-l", "-d", "", "-no-progress", base}, dirs...)
		cmd := exec.CommandContext(ctx, bin, args...)
		var stdout, stderr bytes.Buffer
		cmd.Stdout = &stdout
		cmd.Stderr = &stderr
		if err := cmd.Run(); err != nil {
			return nil, fmt.Errorf("list base directories: %w: %s", err, stderr.String())
		}

		// unsquashfs -l lists every entry under a requested directory
		// recursively; grouping each line by its immediate parent recovers "the
		// direct children of X" for every directory the output mentions at all,
		// including ones nobody asked for — harmless extra data. A directory the
		// base does not have produces nothing, the same answer as one that
		// exists and is empty: neither hides anything.
		for _, line := range strings.Split(stdout.String(), "\n") {
			if !strings.HasPrefix(line, "/") {
				continue
			}
			parent := path.Dir(line)
			out[parent] = append(out[parent], path.Base(line))
		}
		for _, dir := range dirs {
			if _, ok := out[dir]; !ok {
				out[dir] = nil
			}
		}
		for dir := range out {
			sort.Strings(out[dir])
		}
		return out, nil
	}
}

// shellQuote renders a path as a single-quoted shell word.
func shellQuote(s string) string {
	return "'" + strings.ReplaceAll(s, "'", `'\''`) + "'"
}
