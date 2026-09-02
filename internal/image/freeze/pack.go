package freeze

import (
	"bytes"
	"context"
	"fmt"
	"os"
	"path"
	"path/filepath"
	"sort"
	"strings"

	"github.com/Justype/condatainer/internal/artifact/meta"
	"github.com/Justype/condatainer/internal/config"
	"github.com/Justype/condatainer/internal/logging"
	"github.com/Justype/condatainer/internal/runtime/container"
	execpkg "github.com/Justype/condatainer/internal/runtime/exec"
	"github.com/Justype/condatainer/internal/utils"
)

// mountPoint is where the overlay image is mounted inside the pack container.
const mountPoint = "/cnt_freeze"

// metaMountPath is where the staged metadata is bound. mksquashfs names an
// archive root after its source's basename and cannot rename one, so the
// directory must already be called .cnt when mksquashfs sees it.
const metaMountPath = "/" + meta.DirName

// PackOptions is one freeze pack.
type PackOptions struct {
	// Image is the source overlay. It is read, never written.
	Image string
	// Target is the .sqf to write.
	Target string
	// Base is the container the pack runs in — mksquashfs and fuse2fs come from
	// there, so a freeze uses the same tool versions a build does.
	Base string
	// MetaDir is the staged .cnt directory. AppendMeta adds it to a finished
	// archive; Pack itself never carries it, because the manifest cannot be
	// written until the payload it identifies has been packed.
	MetaDir string
	// ApptainerBin is the apptainer to run; empty finds it on PATH.
	ApptainerBin string
	// CompressArgs and BlockSize are mksquashfs tuning, as the build spells them.
	CompressArgs string
	BlockSize    string
	Processors   int
	// StageDir holds a copy of the payload made with rdump, laid out as the image
	// is: stage/upper/... Set, the pack reads that copy and never mounts the
	// image. Empty, it mounts.
	StageDir string
}

// Pack writes the overlay's upper/ layer to a SquashFS artifact, translating
// whiteouts on the way (§2.4a).
//
// By default the payload is read through a fuse2fs mount Apptainer performs in its
// own user namespace, so nothing is copied. The image is made read-only meanwhile:
// fuse2fs decides whether it may write by the permission bits, not by -o ro.
//
// With StageDir it reads a copy made earlier — faster, and costs the payload twice
// on disk. Both routes produce the same archive, given Translation.ForCopy.
func Pack(ctx context.Context, opts PackOptions, entries []Entry, tr Translation) error {
	log := logging.FromContext(ctx)

	img, err := filepath.Abs(opts.Image)
	if err != nil {
		return err
	}
	target, err := filepath.Abs(opts.Target)
	if err != nil {
		return err
	}

	sources := archiveSources(entries, tr)
	if len(sources) == 0 {
		return fmt.Errorf("%s has an empty upper/; there is nothing to freeze", opts.Image)
	}

	// A pack writes four small files and nothing payload-sized: the exclude list,
	// the pseudo-file definitions, and the staged manifest and runtime. They go to
	// the scratch a build uses rather than beside the artifact, which would
	// scatter working directories into an images directory and, on NFS, leave
	// .nfs* placeholders behind.
	scratch := utils.GetTmpDir()
	if err := os.MkdirAll(scratch, 0o755); err != nil {
		return fmt.Errorf("stage translation: %w", err)
	}
	stage, err := os.MkdirTemp(scratch, "cnt-freeze-")
	if err != nil {
		return fmt.Errorf("stage translation: %w", err)
	}
	defer os.RemoveAll(stage)

	// Where the payload's upper/ is readable from inside the pack container, and
	// what it takes to get it there.
	upper := path.Join(mountPoint, UpperDir)
	binds := []string{filepath.Dir(target), stage}
	var apptainerFlags []string
	if opts.StageDir == "" {
		fuse2fs, err := findFuse2fs(ctx, opts.ApptainerBin)
		if err != nil {
			return err
		}
		restore, err := protect(img)
		if err != nil {
			return err
		}
		defer restore()
		apptainerFlags = []string{
			"--fusemount",
			fmt.Sprintf("container:%s -o ro %s %s", fuse2fs, img, mountPoint),
		}
	} else {
		upper = filepath.Join(opts.StageDir, UpperDir)
		binds = append(binds, opts.StageDir)
	}

	args, err := stageTranslation(stage, upper, tr)
	if err != nil {
		return err
	}

	packSources := make([]string, 0, len(sources)+1)
	for _, s := range sources {
		packSources = append(packSources, path.Join(upper, s))
	}

	script := packScript(packSources, target, args, opts)
	route := "mount"
	if opts.StageDir != "" {
		route = "staged copy"
	}
	log.Info("packing frozen environment",
		"source", opts.Image, "target", opts.Target, "route", route,
		"roots", strings.Join(sources, " "), "deletions", tr.Deletions())

	runOpts := execpkg.Options{
		BaseImage:      opts.Base,
		ApptainerBin:   opts.ApptainerBin,
		BindPaths:      container.DeduplicateBindPaths(binds),
		ApptainerFlags: apptainerFlags,
		Command:        []string{"/bin/bash", "-c", script},
		HidePrompt:     true,
	}
	if err := execpkg.Run(ctx, runOpts, execpkg.IOFromContext(ctx)); err != nil {
		os.Remove(target)
		return fmt.Errorf("pack %s: %w", opts.Target, err)
	}
	if _, err := os.Stat(target); err != nil {
		return fmt.Errorf("pack produced no artifact at %s: %w", opts.Target, err)
	}
	return nil
}

// AppendMeta adds the staged .cnt directory to a finished archive.
//
// It is a second mksquashfs run because the manifest records the identity of the
// payload, and that is only knowable once the payload is packed — the manifest
// cannot describe an archive it is already inside. Appending costs a fraction of
// the pack: it adds a handful of small files and rewrites the metadata tables.
func AppendMeta(ctx context.Context, opts PackOptions) error {
	if opts.MetaDir == "" {
		return fmt.Errorf("no metadata staged to append to %s", opts.Target)
	}
	target, err := filepath.Abs(opts.Target)
	if err != nil {
		return err
	}

	// -noappend is deliberately absent: that flag is what makes mksquashfs
	// overwrite, and here the existing archive is the thing being added to.
	//
	// -no-recovery because appending otherwise writes a recovery copy into the
	// process working directory and says so on stdout. It buys nothing here: a
	// failed append deletes the artifact and the freeze is rerun from the
	// overlay, which is still there.
	args := []string{
		metaMountPath, target,
		"-keep-as-directory", "-all-root", "-no-xattrs", "-quiet", "-no-progress", "-no-recovery",
		"-processors", fmt.Sprint(processors(opts.Processors)),
	}
	if opts.BlockSize != "" {
		args = append(args, "-b", opts.BlockSize)
	}
	if opts.CompressArgs != "" {
		args = append(args, opts.CompressArgs)
	}

	runOpts := execpkg.Options{
		BaseImage:    opts.Base,
		ApptainerBin: opts.ApptainerBin,
		BindPaths: container.DeduplicateBindPaths([]string{
			filepath.Dir(target), opts.MetaDir + ":" + metaMountPath,
		}),
		Command:    []string{"/bin/bash", "-c", fmt.Sprintf("trap 'exit 130' INT TERM\nmksquashfs %s\n", strings.Join(args, " "))},
		HidePrompt: true,
	}
	// Silent unless it fails. The append writes two small files after the pack
	// has already reported, so Apptainer's greeting and a 2-entry progress bar
	// are the only things it would say; a failure needs all of both.
	var output bytes.Buffer
	if err := execpkg.Run(ctx, runOpts, execpkg.IO{Stdout: &output, Stderr: &output}); err != nil {
		if said := strings.TrimSpace(output.String()); said != "" {
			return fmt.Errorf("append metadata to %s: %w: %s", opts.Target, err, said)
		}
		return fmt.Errorf("append metadata to %s: %w", opts.Target, err)
	}
	return nil
}

// archiveSources are the top-level entries of upper/, which become the archive's
// roots: upper/cnt_env becomes /cnt_env at mount. Excluded markers are dropped
// here as well as in the exclude file, since a source that is itself a marker
// would otherwise be packed under its own name.
func archiveSources(entries []Entry, tr Translation) []string {
	excluded := map[string]bool{}
	for _, p := range tr.Exclude {
		excluded[p] = true
	}
	var out []string
	for _, e := range entries {
		if e.Dir() != "" || excluded[e.Path] {
			continue
		}
		out = append(out, e.Path)
	}
	sort.Strings(out)
	return out
}

// stageTranslation writes the -ef and -pf files and returns the arguments that
// reference them.
func stageTranslation(stage, upper string, tr Translation) ([]string, error) {
	var args []string
	write := func(name, content, flag string) error {
		if content == "" {
			return nil
		}
		p := filepath.Join(stage, name)
		if err := os.WriteFile(p, []byte(content), 0o644); err != nil {
			return fmt.Errorf("stage %s: %w", name, err)
		}
		args = append(args, flag, p)
		return nil
	}
	if err := write("exclude", tr.ExcludeFile(upper), "-ef"); err != nil {
		return nil, err
	}
	if err := write("pseudo", tr.PseudoFile(), "-pf"); err != nil {
		return nil, err
	}
	return args, nil
}

// protect clears the source's write bit for the duration of the pack and returns
// a function restoring it. An image that is already read-only is left alone,
// including its permissions afterwards — clearing the write bit is how an
// artifact is pinned, and freeze must not quietly unpin one.
func protect(img string) (func(), error) {
	info, err := os.Stat(img)
	if err != nil {
		return nil, err
	}
	mode := info.Mode().Perm()
	if mode&0o222 == 0 {
		return func() {}, nil
	}
	if err := os.Chmod(img, mode&^0o222); err != nil {
		return nil, fmt.Errorf("protect %s for reading: %w", img, err)
	}
	return func() { _ = os.Chmod(img, mode) }, nil
}

// processors is how many compressor threads the pack may use. Never zero and
// never omitted: without -processors mksquashfs takes every core on the machine.
// A freeze budgets them as a build does (build.ncpus).
func processors(n int) int {
	if n > 0 {
		return n
	}
	return config.DefaultNcpus
}

// packScript renders the mksquashfs invocation. -no-xattrs for the reason
// build/squashfs.go gives; nothing here needs them, since an opaque directory's
// xattr is re-expressed as whiteouts precisely because it could not be carried.
func packScript(sources []string, target string, translation []string, opts PackOptions) string {
	args := []string{
		strings.Join(sources, " "), target,
		"-noappend", "-keep-as-directory", "-all-root", "-no-xattrs", "-quiet",
	}
	args = append(args, "-processors", fmt.Sprint(processors(opts.Processors)))
	if opts.BlockSize != "" {
		args = append(args, "-b", opts.BlockSize)
	}
	args = append(args, translation...)
	if opts.CompressArgs != "" {
		args = append(args, opts.CompressArgs)
	}
	return fmt.Sprintf("trap 'exit 130' INT TERM\nmksquashfs %s\n", strings.Join(args, " "))
}
