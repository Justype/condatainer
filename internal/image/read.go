package image

import (
	"bytes"
	"fmt"
	"os/exec"
	"strings"

	"github.com/Justype/condatainer/internal/image/sif"
	"github.com/Justype/condatainer/internal/image/squashfs"
	"github.com/Justype/condatainer/internal/image/tool"
	"github.com/Justype/condatainer/internal/utils"
)

// ReadFile reads a file at innerPath from inside an image.
// For .sqf images it uses unsquashfs -cat; for .img images it uses debugfs.
// It returns nil if the file cannot be read.
func ReadFile(imagePath, innerPath string) []byte {
	if utils.IsSqf(imagePath) {
		return squashfs.Cat(imagePath, innerPath)
	}
	if utils.IsImg(imagePath) {
		return imgCat(imagePath, innerPath)
	}
	return nil
}

// imgCat reads a file from inside an ext3 overlay image using debugfs.
// Writable content lives under upper/, so /cnt_env/conda-meta/history is read
// as upper/cnt_env/conda-meta/history. It returns nil if the file is absent.
func imgCat(imgPath, innerPath string) []byte {
	dbg, err := exec.LookPath("debugfs")
	if err != nil {
		return nil
	}
	inner := strings.TrimPrefix(innerPath, "/")
	catArg := "cat upper/" + inner
	cmd := exec.Command(dbg, "-R", catArg, imgPath)
	var stderr bytes.Buffer
	cmd.Stderr = &stderr
	out, err := cmd.Output()
	if err != nil {
		return nil
	}
	if len(out) == 0 && strings.Contains(stderr.String(), "File not found") {
		return nil
	}
	return out
}

// ExtractDir extracts a directory out of an image into destDir, keeping the
// archive-relative path: extracting "/.cnt" yields destDir/.cnt/….
//
// It exists so a caller that wants several files from one directory pays for one
// unsquashfs process instead of one per file. A writable .img is not handled: it
// carries no embedded metadata.
func ExtractDir(imagePath, innerPath, destDir string) error {
	switch {
	case utils.IsSqf(imagePath):
		return squashfs.ExtractDir(imagePath, innerPath, destDir, 0)
	case utils.IsSif(imagePath):
		part, err := sif.PrimarySystemPartition(imagePath)
		if err != nil {
			return err
		}
		return squashfs.ExtractDir(imagePath, innerPath, destDir, part.Offset)
	default:
		return fmt.Errorf("%w: %s is not a readable image", tool.ErrCorrupt, imagePath)
	}
}
