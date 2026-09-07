// Package freeze converts a writable ext3 overlay into an immutable SquashFS
// artifact, and back.
//
// Everything here reads the image as a filesystem, never as a mounted overlay.
// What a process sees inside a container is the union of the base and the
// overlay, so packing from there would pack the base; the delta exists in exactly
// one place, the image's own upper/ directory, and reading it needs no inference
// and no comparison.
package freeze

import (
	"bufio"
	"bytes"
	"context"
	"fmt"
	"os/exec"
	"path"
	"strconv"
	"strings"

	"github.com/Justype/condatainer/internal/image/tool"
	"github.com/Justype/condatainer/internal/toolpath"
)

// UpperDir is the OverlayFS upper layer inside an ext3 overlay image. Everything
// the user has ever written lives beneath it.
const UpperDir = "/upper"

// Entry is one thing found in the overlay's upper layer.
//
// Path is relative to upper/, so it is also the path the payload will have once
// the wrapper is stripped: upper/cnt_env/bin/tool becomes cnt_env/bin/tool in the
// archive and /cnt_env/bin/tool at mount.
type Entry struct {
	Path string
	Mode uint32 // as debugfs reports it: 0100644, 040755, 020000, 0120777
	Size int64
}

// IsDir reports a directory.
func (e Entry) IsDir() bool { return e.Mode&0o170000 == 0o040000 }

// IsCharDev reports a character device. A char device 0:0 in an overlay upper
// layer is a whiteout — the record that a file present in the base was deleted.
func (e Entry) IsCharDev() bool { return e.Mode&0o170000 == 0o020000 }

// Name is the entry's final component.
func (e Entry) Name() string { return path.Base(e.Path) }

// Dir is the entry's parent, relative to upper/, or "" at the top level.
func (e Entry) Dir() string {
	d := path.Dir(e.Path)
	if d == "." || d == "/" {
		return ""
	}
	return d
}

// Walk lists everything under upper/ in the image, without mounting it.
//
// debugfs has no recursive listing, so this walks level by level, batching one
// session per level. Reading never writes, which is what lets enumeration run
// against an image whose write bit has not been cleared.
func Walk(ctx context.Context, imgPath string) ([]Entry, error) {
	debugfsPath, err := toolpath.Resolve("debugfs")
	if err != nil {
		return nil, err
	}

	var out []Entry
	level := []string{""} // relative to upper/; "" is upper/ itself

	for len(level) > 0 {
		if err := ctx.Err(); err != nil {
			return nil, err
		}
		entries, err := listDirs(ctx, debugfsPath, imgPath, level)
		if err != nil {
			return nil, err
		}
		var next []string
		for _, e := range entries {
			out = append(out, e)
			if e.IsDir() {
				next = append(next, e.Path)
			}
		}
		level = next
	}
	return out, nil
}

// listDirs runs one debugfs session listing every directory in dirs.
func listDirs(ctx context.Context, debugfsPath, imgPath string, dirs []string) ([]Entry, error) {
	var script bytes.Buffer
	for _, d := range dirs {
		fmt.Fprintf(&script, "ls -p %s\n", path.Join(UpperDir, d))
	}
	script.WriteString("quit\n")

	cmd := exec.CommandContext(ctx, debugfsPath, imgPath)
	cmd.Stdin = &script
	var stdout, stderr bytes.Buffer
	cmd.Stdout = &stdout
	cmd.Stderr = &stderr
	if err := cmd.Run(); err != nil {
		return nil, &tool.Error{
			Op: "list overlay", Path: imgPath, Tool: "debugfs",
			Output: stderr.String(), BaseErr: err,
		}
	}
	return parseListing(stdout.String(), dirs)
}

// debugfsPrompt precedes every command debugfs echoes, and so separates one
// directory's entries from the next.
const debugfsPrompt = "debugfs:"

// parseListing reads batched `ls -p` output, attributing each block to the
// directory asked for in that position.
//
// The echoed command is counted, never parsed for its path: debugfs uses
// readline, which redraws a long command and leaves cursor sequences in the
// middle of the echo, or truncates its front under TERM=dumb. Order is reliable
// instead — one prompt per command, in order.
func parseListing(out string, dirs []string) ([]Entry, error) {
	var entries []Entry
	idx := -1
	scanner := bufio.NewScanner(strings.NewReader(out))
	scanner.Buffer(make([]byte, 0, 64*1024), 4*1024*1024)

	for scanner.Scan() {
		line := scanner.Text()
		if strings.HasPrefix(line, debugfsPrompt) {
			idx++
			continue
		}
		if idx < 0 || idx >= len(dirs) {
			continue
		}
		e, ok := parseEntry(strings.TrimSpace(line), dirs[idx])
		if !ok {
			continue
		}
		entries = append(entries, e)
	}
	if err := scanner.Err(); err != nil {
		return nil, err
	}
	// Fewer blocks than commands means the output cannot be attributed, and
	// guessing would silently misplace every entry after the gap.
	if idx < len(dirs) {
		return nil, fmt.Errorf("debugfs listed %d of %d directories; its output cannot be attributed", idx+1, len(dirs))
	}
	return entries, nil
}

func parseEntry(line, dir string) (Entry, bool) {
	if !strings.HasPrefix(line, "/") {
		return Entry{}, false
	}
	// A leading empty field, then inode/mode/uid/gid/name/size and a trailing one.
	f := strings.Split(strings.ReplaceAll(line, " ", ""), "/")
	if len(f) < 7 {
		return Entry{}, false
	}
	name := f[5]
	if name == "" || name == "." || name == ".." {
		return Entry{}, false
	}
	mode, err := strconv.ParseUint(f[2], 8, 32)
	if err != nil {
		return Entry{}, false
	}
	size, _ := strconv.ParseInt(f[6], 10, 64)
	return Entry{Path: path.Join(dir, name), Mode: uint32(mode), Size: size}, true
}
