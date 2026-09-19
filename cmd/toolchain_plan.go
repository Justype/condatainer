package cmd

import (
	"context"
	"fmt"

	"github.com/Justype/condatainer/internal/image/freeze"
	"github.com/Justype/condatainer/internal/libexec"
	"github.com/Justype/condatainer/internal/runtime/apptainer"
	"github.com/Justype/condatainer/internal/toolpath"
)

// systemProbe is what the host offers for condatainer's own tools.
type systemProbe struct {
	ApptainerErr error  // nil when a usable apptainer exists outside libexec
	Apptainer    string // its version, for the report
	Squashfs     bool   // mksquashfs and unsquashfs both resolve
	Squashfuse   bool   // squashfuse or squashfuse_ll resolves
	Fuse2fs      bool   // fuse2fs resolves

	// SystemApptainerErr is why no system/module apptainer or singularity can be
	// run, nil when one can. Building an os overlay needs it (libexec's
	// non-setuid apptainer cannot do that build).
	SystemApptainerErr error
}

// planToolchain decides which libexec packages the system lacks, and describes
// what it found. micromamba is always installed and is not listed. fuse2fs and
// a missing system apptainer cannot be provided by libexec, so they are only
// reported: the second as a warning.
func planToolchain(p systemProbe) (install, report, warnings []string) {
	if p.ApptainerErr != nil {
		install = append(install, "apptainer") // brings squashfs-tools and squashfuse with it
		report = append(report, fmt.Sprintf("apptainer: none usable (%v) -> install apptainer, with its squashfs tools", p.ApptainerErr))
	} else {
		report = append(report, "apptainer: "+p.Apptainer+" ok")
		if p.Squashfs {
			report = append(report, "mksquashfs/unsquashfs: ok")
		} else {
			install = append(install, "squashfs-tools")
			report = append(report, "mksquashfs/unsquashfs: not found -> install squashfs-tools")
		}
		if p.Squashfuse {
			report = append(report, "squashfuse: ok")
		} else {
			install = append(install, "squashfuse")
			report = append(report, "squashfuse: not found -> install squashfuse")
		}
	}
	if !p.Fuse2fs {
		report = append(report, "fuse2fs: not found; it cannot be installed here, and writable .img overlays need it")
	}
	if p.SystemApptainerErr != nil {
		warnings = append(warnings, fmt.Sprintf("no system apptainer found (%v): it is required to build os overlays", p.SystemApptainerErr))
	}
	return install, report, warnings
}

// probeSystem asks the same resolvers the tools themselves use.
func probeSystem(ctx context.Context) systemProbe {
	p := systemProbe{SystemApptainerErr: apptainer.EnsureApptainer()}
	if libexec.Installed("apptainer") {
		p.Apptainer = "installed in libexec"
	} else if p.ApptainerErr = apptainer.CheckSystemBin(); p.ApptainerErr == nil {
		p.Apptainer = "found"
		if _, v, err := apptainer.Current(); err == nil {
			p.Apptainer = v
		}
	}
	_, mkErr := toolpath.Resolve("mksquashfs")
	_, unErr := toolpath.Resolve("unsquashfs")
	p.Squashfs = mkErr == nil && unErr == nil
	_, sfErr := freeze.FindSquashfuse()
	p.Squashfuse = sfErr == nil
	_, fuseErr := toolpath.Resolve("fuse2fs")
	p.Fuse2fs = fuseErr == nil
	return p
}
