package freeze

import (
	"context"
	"errors"
	"fmt"
	"os/exec"
	"syscall"

	"github.com/Justype/condatainer/internal/image/tool"
	"github.com/Justype/condatainer/internal/utils"
)

// ErrNoSpace reports scratch that cannot hold the staged payload.
var ErrNoSpace = errors.New("not enough scratch space to stage the payload")

// stageMarginMB is headroom over the payload: rdump writes into a filesystem
// whose block size and metadata cost are not the source's, and a stage that
// fills its filesystem is worse than one that was never started.
const stageMarginMB = 256

// dumpUpper copies the overlay's upper/ into stage, producing stage/upper. It
// reads the image as a filesystem, so it needs no mount and no privilege — but it
// cannot create a device node, which is why the copy route needs
// Translation.ForCopy.
func dumpUpper(ctx context.Context, img, stage string) error {
	cmd := exec.CommandContext(ctx, "debugfs", "-R", "rdump "+UpperDir+" "+stage, img)
	out, err := cmd.CombinedOutput()
	if err != nil {
		return &tool.Error{Op: "dump", Path: img, Tool: "debugfs", Output: string(out), BaseErr: err}
	}
	return nil
}

// checkStageSpace refuses a copy that will not fit, before anything is written.
//
// A staged freeze writes the whole payload a second time, and finding that out
// partway through leaves a half-copied tree in someone's scratch and no artifact.
func checkStageSpace(dir string, payloadMB int) error {
	var st syscall.Statfs_t
	if err := syscall.Statfs(dir, &st); err != nil {
		return fmt.Errorf("check space on %s: %w", dir, err)
	}
	availMB := int(st.Bavail * uint64(st.Bsize) / (1024 * 1024))
	needMB := payloadMB + stageMarginMB
	if availMB < needMB {
		return fmt.Errorf("%w: %s holds %s free, staging needs %s — freeze without --use-tmp to read the image directly",
			ErrNoSpace, dir, utils.FormatSize(int64(availMB)*1024*1024), utils.FormatSize(int64(needMB)*1024*1024))
	}
	return nil
}
