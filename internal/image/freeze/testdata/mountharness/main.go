// Command mountharness is a test-only stand-in for the condatainer binary,
// used by internal/image/freeze's mount tests. MountedRun re-execs
// os.Executable() into the hidden _mount_sentinel command; the `go test`
// binary can't play that role since it has no cmd/cobra dispatch at all, so
// tests point MountedRun at this binary instead.
//
// Usage:
//
//	mountharness _mount_sentinel <fuseBin> <mnt> [fuseArgs...]   see freeze.RunSentinel; work comes from $CNT_MOUNT_WORK
//	mountharness drive <fuseBin> <sqf> <mnt> <work>              calls freeze.MountedRun directly and blocks until it returns
package main

import (
	"context"
	"fmt"
	"os"

	"github.com/Justype/condatainer/internal/image/freeze"
	execpkg "github.com/Justype/condatainer/internal/runtime/exec"
)

func main() {
	if len(os.Args) < 2 {
		fail("usage: mountharness _mount_sentinel <fuseBin> <mnt> [fuseArgs...] | drive <fuseBin> <sqf> <mnt> <work>")
	}

	switch os.Args[1] {
	case "_mount_sentinel":
		if len(os.Args) < 4 {
			fail("usage: mountharness _mount_sentinel <fuseBin> <mnt> [fuseArgs...]")
		}
		fuseBin, mnt := os.Args[2], os.Args[3]
		fuseArgs := os.Args[4:]
		if err := freeze.RunSentinel(fuseBin, mnt, os.Getenv(freeze.EnvSentinelWork), fuseArgs); err != nil {
			fail("%v", err)
		}
	case "drive":
		if len(os.Args) != 6 {
			fail("usage: mountharness drive <fuseBin> <sqf> <mnt> <work>")
		}
		io := execpkg.IO{Stdin: os.Stdin, Stdout: os.Stdout, Stderr: os.Stderr}
		if err := freeze.MountedRun(context.Background(), os.Args[2], []string{os.Args[3]}, os.Args[4], os.Args[5], io); err != nil {
			fail("%v", err)
		}
	default:
		fail("unknown mode: %s", os.Args[1])
	}
}

func fail(format string, a ...any) {
	fmt.Fprintf(os.Stderr, format+"\n", a...)
	os.Exit(1)
}
