package server

import (
	"context"
	"fmt"
	"io"
	"io/fs"
	"os"
	"path/filepath"
	"strings"

	"github.com/Justype/condatainer/internal/utils"
)

// countTreeEntries returns the number of filesystem entries under root
// (root itself plus everything beneath it), used to size a progress bar
// before copyTree starts.
func countTreeEntries(root string) (int, error) {
	total := 0
	err := filepath.WalkDir(root, func(_ string, _ fs.DirEntry, err error) error {
		if err != nil {
			return err
		}
		total++
		return nil
	})
	return total, err
}

// copyTree recursively copies src to dst, preserving each entry's
// permission bits like cp -R; with preserveTimes it also preserves
// modification times, matching what mv does for a cross-filesystem move.
// Symlinks are recreated as symlinks (not followed); FIFOs, sockets, and
// device nodes are skipped. onEntry is called once per entry, for progress
// reporting. Aborts as soon as ctx is done, leaving a partial copy at dst
// for the caller to clean up.
func copyTree(ctx context.Context, src, dst string, preserveTimes bool, onEntry func()) error {
	// Directories are created with owner-rwx added so their contents can be
	// written (a 0555 source dir would otherwise block its own children),
	// then fixed up children-first after the walk — which for preserveTimes
	// also keeps dir mtimes from being clobbered by the writes inside them.
	type dirFix struct {
		target string
		info   os.FileInfo
	}
	var dirs []dirFix

	err := filepath.WalkDir(src, func(p string, d fs.DirEntry, err error) error {
		if err != nil {
			return err
		}
		if ctx.Err() != nil {
			return ctx.Err()
		}
		rel, err := filepath.Rel(src, p)
		if err != nil {
			return err
		}
		target := filepath.Join(dst, rel)

		switch {
		case d.Type()&os.ModeSymlink != 0:
			linkTarget, err := os.Readlink(p)
			if err != nil {
				return err
			}
			if err := os.Symlink(linkTarget, target); err != nil {
				return err
			}
		case d.IsDir():
			info, err := d.Info()
			if err != nil {
				return err
			}
			if err := os.MkdirAll(target, info.Mode().Perm()|0o700); err != nil {
				return err
			}
			dirs = append(dirs, dirFix{target, info})
		case !d.Type().IsRegular():
			// FIFO/socket/device: opening one can block forever; skip it
			// (still counted for progress)
		default:
			info, err := d.Info()
			if err != nil {
				return err
			}
			if err := copyFileEntry(ctx, p, target, info, preserveTimes); err != nil {
				return err
			}
		}
		onEntry()
		return nil
	})
	if err != nil {
		return err
	}

	for i := len(dirs) - 1; i >= 0; i-- {
		os.Chmod(dirs[i].target, dirs[i].info.Mode().Perm()) //nolint:errcheck
		if preserveTimes {
			os.Chtimes(dirs[i].target, dirs[i].info.ModTime(), dirs[i].info.ModTime()) //nolint:errcheck
		}
	}
	return nil
}

// copyFileEntry copies the regular file at src to dst with the source's
// permission bits (and, with preserveTimes, its mtime), creating dst's
// parent if missing. Copies in chunks and aborts mid-file when ctx is
// done, so cancelling doesn't have to wait for a huge file to finish.
func copyFileEntry(ctx context.Context, src, dst string, info os.FileInfo, preserveTimes bool) error {
	in, err := os.Open(src)
	if err != nil {
		return err
	}
	defer in.Close()

	if err := os.MkdirAll(filepath.Dir(dst), utils.PermDir); err != nil {
		return err
	}
	out, err := os.OpenFile(dst, os.O_CREATE|os.O_WRONLY|os.O_TRUNC, info.Mode().Perm())
	if err != nil {
		return err
	}
	defer out.Close()
	out.Chmod(info.Mode().Perm()) //nolint:errcheck // exact bits — the create perm is masked by umask

	buf := make([]byte, 1<<20)
	for {
		if err := ctx.Err(); err != nil {
			return err
		}
		n, rerr := in.Read(buf)
		if n > 0 {
			if _, werr := out.Write(buf[:n]); werr != nil {
				return werr
			}
		}
		if rerr == io.EOF {
			break
		}
		if rerr != nil {
			return rerr
		}
	}
	if preserveTimes {
		os.Chtimes(dst, info.ModTime(), info.ModTime()) //nolint:errcheck
	}
	return nil
}

// removeTree recursively removes path (symlinks are removed as entries,
// not followed), calling onEntry after each removed entry and aborting
// between entries when ctx is done. Entries that disappear mid-delete are
// treated as already removed, matching os.RemoveAll.
func removeTree(ctx context.Context, path string, onEntry func()) error {
	if err := ctx.Err(); err != nil {
		return err
	}
	info, err := os.Lstat(path)
	if os.IsNotExist(err) {
		return nil
	}
	if err != nil {
		return err
	}
	if info.IsDir() {
		entries, err := os.ReadDir(path)
		if err != nil && !os.IsNotExist(err) {
			return err
		}
		for _, e := range entries {
			if err := removeTree(ctx, filepath.Join(path, e.Name()), onEntry); err != nil {
				return err
			}
		}
	}
	if err := os.Remove(path); err != nil && !os.IsNotExist(err) {
		return err
	}
	onEntry()
	return nil
}

// uniquePath returns path if nothing exists there, otherwise the first
// free "name_copy" / "name_copyN" variant (inserted before the extension
// for files; appended to the whole name for directories and dotfiles).
func uniquePath(path string, isDir bool) string {
	if _, err := os.Lstat(path); err != nil {
		return path
	}
	dir := filepath.Dir(path)
	name := filepath.Base(path)
	ext := ""
	if dot := strings.LastIndex(name, "."); !isDir && dot > 0 {
		ext = name[dot:]
		name = name[:dot]
	}
	for n := 1; ; n++ {
		suffix := "_copy"
		if n > 1 {
			suffix = fmt.Sprintf("_copy%d", n)
		}
		p := filepath.Join(dir, name+suffix+ext)
		if _, err := os.Lstat(p); err != nil {
			return p
		}
	}
}

// destInsideSource reports whether destDir is srcDir itself or inside its
// subtree, resolving symlinks so a linked destination can't dodge the
// check.
func destInsideSource(srcDir, destDir string) bool {
	src, err := filepath.EvalSymlinks(srcDir)
	if err != nil {
		src = filepath.Clean(srcDir)
	}
	dst, err := filepath.EvalSymlinks(destDir)
	if err != nil {
		dst = filepath.Clean(destDir)
	}
	sep := string(filepath.Separator)
	return dst == src || strings.HasPrefix(dst+sep, src+sep)
}
