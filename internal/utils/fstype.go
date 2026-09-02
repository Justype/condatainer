package utils

import (
	"sync"
	"syscall"
)

// Filesystem magic numbers, from Linux statfs(2). Only the ones worth naming in
// a message are here; anything else is reported by its number.
const (
	fsAFS     = 0x5346414F
	fsBeeGFS  = 0x19830326
	fsBtrfs   = 0x9123683E
	fsCeph    = 0x00C36400
	fsCIFS    = 0xFF534D42
	fsExt     = 0xEF53
	fsGPFS    = 0x47504653
	fsLustre  = 0x0BD00BD0
	fsOverlay = 0x794C7630
	fsPanFS   = 0xAAD7AAEA
	fsTmpfs   = 0x01021994
	fsNFS     = 0x6969
	fsXFS     = 0x58465342
	fsZFS     = 0x2FC12FC1
)

// fsNames is what each magic number is called.
var fsNames = map[int64]string{
	fsAFS: "afs", fsBeeGFS: "beegfs", fsBtrfs: "btrfs", fsCeph: "ceph",
	fsCIFS: "cifs", fsExt: "ext", fsGPFS: "gpfs", fsLustre: "lustre",
	fsOverlay: "overlay", fsPanFS: "panfs", fsTmpfs: "tmpfs", fsNFS: "nfs",
	fsXFS: "xfs", fsZFS: "zfs",
}

// networkFS is which of them are reached over a network, where a workload of
// many small files pays a round trip per file.
var networkFS = map[int64]bool{
	fsAFS: true, fsBeeGFS: true, fsCeph: true, fsCIFS: true,
	fsGPFS: true, fsLustre: true, fsNFS: true, fsPanFS: true,
}

// FSKind names the filesystem a path is on and says whether it is reached over a
// network. An unreadable path reports nothing rather than failing: this exists to
// improve a message, never to decide anything.
func FSKind(path string) (name string, network bool) {
	var st syscall.Statfs_t
	if err := syscall.Statfs(path, &st); err != nil {
		return "", false
	}
	t := int64(st.Type)
	name = fsNames[t]
	if name == "" {
		name = "unknown"
	}
	return name, networkFS[t]
}

var scratchWarned sync.Once

// WarnNetworkScratch says once that the scratch a job is about to fill is on a
// network filesystem. It only reports: where temporary files belong is answered
// by CNT_TMPDIR, the scheduler, then TMPDIR, and an administrator who sets one
// has done so deliberately.
func WarnNetworkScratch(dir, what string) {
	name, network := FSKind(dir)
	if !network {
		return
	}
	scratchWarned.Do(func() {
		PrintWarning("%s is staging to %s, which is %s; many small files are slower over a network.", what, dir, name)
		PrintWarning("Set CNT_TMPDIR to a local path if this machine has one.")
	})
}
