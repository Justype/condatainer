package freeze

import (
	"bufio"
	"bytes"
	"context"
	"fmt"
	"os/exec"
	"path"
	"sort"
	"strings"

	"github.com/Justype/condatainer/internal/image/tool"
	"github.com/Justype/condatainer/internal/logging"
)

// Whiteout conventions. Which one an overlay carries is decided by the driver
// Apptainer used, not by anything CondaTainer chose:
//
//	kernel overlay   a deleted file is a char device 0:0
//	                 a replaced directory carries trusted.overlay.opaque
//	fuse-overlayfs   a deleted file is a .wh.<name> plain file
//	                 a replaced directory carries .wh..wh..opq, and often
//	                 user.fuseoverlayfs.opaque as well
const (
	whPrefix    = ".wh."
	opqMarker   = ".wh..wh..opq"
	opqXattrKey = "overlay.opaque" // matches trusted.* and user.fuseoverlayfs.*
)

// Convention is the shape of the whiteouts an overlay was found to carry. It is
// recorded in the artifact as build diagnostics: nothing compares it, but the
// first person to debug a resurrected file needs to know which they are reading.
type Convention string

const (
	ConventionNone      Convention = "none"    // no deletions
	ConventionCharDev   Convention = "chardev" // kernel overlayfs
	ConventionWhiteFile Convention = "whfile"  // fuse-overlayfs
	ConventionMixed     Convention = "mixed"   // both, which one overlay can hold
)

// Translation is what the pack must do differently because of deletions.
//
// A marker file means "deleted" only in a writable upper layer; a frozen artifact
// is always a read-only lower one, where a char device 0:0 is what means it. So
// freeze rewrites deletions rather than carrying them across.
type Translation struct {
	// Exclude are paths relative to upper/ that must not be packed: the marker
	// files themselves, which have no meaning in the artifact.
	Exclude []string
	// Pseudo are mksquashfs pseudo-file definitions, one per whiteout, in the
	// form "<path> c <mode> <uid> <gid> <major> <minor>". They put a device node
	// in the archive without creating one on disk, which is what makes the whole
	// translation possible unprivileged.
	Pseudo []string
	// Devices are whiteouts the image already holds as char 0:0 nodes, by path
	// relative to upper/. They need no definition when the pack reads the image
	// directly — mksquashfs copies the node as it stands — so they are listed
	// rather than translated. The copy route cannot: see ForCopy.
	Devices []string
	// Convention is what the source overlay used.
	Convention Convention
	// Opaque are the directories whose contents replaced a base directory's.
	Opaque []string
}

// ForCopy returns the translation a pack needs when it reads a staged copy rather
// than the image itself. rdump cannot create a device node, so whiteouts the
// image held as char 0:0 are absent from the copy and must be injected.
func (t Translation) ForCopy() Translation {
	if len(t.Devices) == 0 {
		return t
	}
	out := t
	out.Pseudo = make([]string, 0, len(t.Pseudo)+len(t.Devices))
	out.Pseudo = append(out.Pseudo, t.Pseudo...)
	for _, p := range t.Devices {
		out.Pseudo = append(out.Pseudo, pseudoWhiteout(p))
	}
	sort.Strings(out.Pseudo)
	out.Devices = nil
	return out
}

// BaseLister answers what each of several directories holds in the base image.
// It takes them all at once because answering means entering the base container,
// and fuse-overlayfs marks every directory an overlay creates as opaque.
type BaseLister func(ctx context.Context, dirs []string) (map[string][]string, error)

// Translate reads an overlay's whiteouts and returns what the pack must do.
//
// The scan cannot be by filename alone: an opaque directory is a marker file on
// fuse-overlayfs and only an xattr on kernel overlayfs. Both are read from the raw
// filesystem with debugfs, so neither needs a mount or a privilege.
func Translate(ctx context.Context, imgPath string, entries []Entry, base BaseLister) (Translation, error) {
	t := Translation{Convention: ConventionNone}
	seenChar, seenFile := false, false
	opaque := map[string]bool{}

	for _, e := range entries {
		switch {
		case e.IsCharDev():
			// Already the form the artifact needs; pack it as it stands.
			seenChar = true
			t.Devices = append(t.Devices, e.Path)
		case e.Name() == opqMarker:
			seenFile = true
			t.Exclude = append(t.Exclude, e.Path)
			opaque[e.Dir()] = true
		case strings.HasPrefix(e.Name(), whPrefix):
			seenFile = true
			t.Exclude = append(t.Exclude, e.Path)
			t.Pseudo = append(t.Pseudo, pseudoWhiteout(path.Join(e.Dir(), strings.TrimPrefix(e.Name(), whPrefix))))
		}
	}

	xattrOpaque, err := opaqueByXattr(ctx, imgPath, entries)
	if err != nil {
		return Translation{}, err
	}
	for _, dir := range xattrOpaque {
		if !opaque[dir] {
			seenFile = true
			opaque[dir] = true
		}
	}

	// One whiteout per base entry the upper does not itself provide. This is
	// exact against the base it was frozen against and only that one: a
	// trusted.overlay.opaque xattr would be base-independent, but mksquashfs has
	// no pseudo syntax for an xattr and the staged route needs a filesystem that
	// stores them, which a shared NFS does not.
	present := map[string]bool{}
	for _, e := range entries {
		present[e.Path] = true
	}
	for dir := range opaque {
		t.Opaque = append(t.Opaque, dir)
	}
	sort.Strings(t.Opaque)

	if len(t.Opaque) > 0 {
		logging.FromContext(ctx).Info("resolving opaque directories against the base",
			"directories", len(t.Opaque))
		absolute := make([]string, 0, len(t.Opaque))
		for _, dir := range t.Opaque {
			absolute = append(absolute, "/"+dir)
		}
		contents, err := base(ctx, absolute)
		if err != nil {
			return Translation{}, fmt.Errorf("list base directories: %w", err)
		}
		for _, dir := range t.Opaque {
			for _, name := range contents["/"+dir] {
				p := path.Join(dir, name)
				if present[p] {
					continue // the overlay provides its own; nothing to hide
				}
				t.Pseudo = append(t.Pseudo, pseudoWhiteout(p))
			}
		}
	}

	switch {
	case seenChar && seenFile:
		t.Convention = ConventionMixed
	case seenChar:
		t.Convention = ConventionCharDev
	case seenFile:
		t.Convention = ConventionWhiteFile
	}

	sort.Strings(t.Exclude)
	sort.Strings(t.Pseudo)
	sort.Strings(t.Devices)
	return t, nil
}

// pseudoWhiteout renders one mksquashfs pseudo-file definition for a whiteout:
// a character device 0:0 with mode 0, which is what OverlayFS reads as "deleted".
func pseudoWhiteout(p string) string { return p + " c 0 0 0 0 0" }

// opaqueByXattr finds directories marked opaque by extended attribute rather
// than by a marker file.
//
// The kernel hides trusted.* from an unprivileged reader through a mount, which
// is why this reads the raw filesystem instead: debugfs is not subject to that
// restriction, so freeze needs no privilege on either driver.
func opaqueByXattr(ctx context.Context, imgPath string, entries []Entry) ([]string, error) {
	var dirs []string
	for _, e := range entries {
		if e.IsDir() {
			dirs = append(dirs, e.Path)
		}
	}
	if len(dirs) == 0 {
		return nil, nil
	}

	var script bytes.Buffer
	for _, d := range dirs {
		fmt.Fprintf(&script, "ea_list %s\n", path.Join(UpperDir, d))
	}
	script.WriteString("quit\n")

	cmd := exec.CommandContext(ctx, "debugfs", imgPath)
	cmd.Stdin = &script
	var stdout, stderr bytes.Buffer
	cmd.Stdout = &stdout
	cmd.Stderr = &stderr
	if err := cmd.Run(); err != nil {
		return nil, &tool.Error{
			Op: "read overlay xattrs", Path: imgPath, Tool: "debugfs",
			Output: stderr.String(), BaseErr: err,
		}
	}

	var out []string
	dir := ""
	scanner := bufio.NewScanner(&stdout)
	scanner.Buffer(make([]byte, 0, 64*1024), 4*1024*1024)
	for scanner.Scan() {
		line := scanner.Text()
		if i := strings.Index(line, "ea_list "); i >= 0 {
			dir = strings.TrimSpace(line[i+len("ea_list "):])
			dir = strings.TrimPrefix(strings.TrimPrefix(dir, UpperDir), "/")
			continue
		}
		if dir != "" && strings.Contains(line, opqXattrKey) {
			out = append(out, dir)
			dir = "" // one hit per directory is enough
		}
	}
	return out, scanner.Err()
}

// ExcludeFile renders the -ef list mksquashfs reads, one path per line, each
// resolved against the mount point the pack reads the payload from.
func (t Translation) ExcludeFile(mountUpper string) string {
	var b strings.Builder
	for _, p := range t.Exclude {
		b.WriteString(path.Join(mountUpper, p))
		b.WriteString("\n")
	}
	return b.String()
}

// PseudoFile renders the -pf definitions mksquashfs reads.
func (t Translation) PseudoFile() string {
	if len(t.Pseudo) == 0 {
		return ""
	}
	return strings.Join(t.Pseudo, "\n") + "\n"
}

// Deletions is how many whiteouts the artifact will carry, by either route: the
// ones being translated plus the ones already in the form the archive needs.
func (t Translation) Deletions() int { return len(t.Pseudo) + len(t.Devices) }
