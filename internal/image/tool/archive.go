package tool

import "errors"

// Sentinels for reading one file out of an image archive.
//
// A caller has to tell "this image does not contain that file" apart from "this
// image could not be read at all": the first is routine — an image built before
// the metadata format existed simply has no manifest — while the second means
// something is wrong with the host or the file. Collapsing them turns a missing
// unsquashfs into "no metadata" and hides the real fault. A missing tool itself
// is internal/toolpath.ErrToolMissing — not defined here, since finding a
// runnable tool at all is not an image-archive concern.
var (
	ErrFileNotFound = errors.New("file not found in image")
	ErrUnreadable   = errors.New("image could not be opened")
	ErrCorrupt      = errors.New("image is corrupt or not a supported format")
)
