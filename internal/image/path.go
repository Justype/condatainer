package image

import (
	"os/exec"
	"strings"

	"github.com/Justype/condatainer/internal/image/squashfs"
	"github.com/Justype/condatainer/internal/utils"
)

// imgPathExists checks whether entry exists inside an ext3 image.
// It uses debugfs to stat the upper/<entry> path in the image.
func imgPathExists(imgPath, entry string) bool {
	entry = strings.TrimPrefix(entry, "/")
	dbg, err := exec.LookPath("debugfs")
	if err != nil {
		return false
	}
	statArg := "stat upper/" + entry
	cmd := exec.Command(dbg, "-R", statArg, imgPath)
	out, err := cmd.CombinedOutput()
	if err != nil {
		outStr := string(out)
		if strings.Contains(outStr, "File not found") || strings.Contains(outStr, "No such file") {
			return false
		}
		return false
	}
	return true
}

// PathExists dispatches to the appropriate backend check depending on the image type.
func PathExists(imagePath, entry string) bool {
	if utils.IsSqf(imagePath) {
		return squashfs.PathExists(imagePath, entry)
	}
	if utils.IsImg(imagePath) {
		return imgPathExists(imagePath, entry)
	}
	return false
}
