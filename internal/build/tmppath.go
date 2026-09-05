package build

import (
	"github.com/Justype/condatainer/catalog"
	"github.com/Justype/condatainer/internal/config"
	"github.com/Justype/condatainer/internal/utils"
)

// A build works under one of two roots — fast local scratch (utils.GetTmpDir,
// which CNT_TMPDIR and then TMPDIR select) or the stable writable tmp
// (config.GetWritableTmpDir). See the README's Workspace Strategy for which
// build gets which, and why.

// tmpRootForType picks the root for a catalog build: the stable writable tmp for
// data, fast local scratch for everything else. A definition build overrides
// this once its type is known — see tmpRootForDef.
func tmpRootForType(typ catalog.Type) string {
	if typ == catalog.TypeData {
		return config.GetWritableTmpDir()
	}
	scratch := utils.GetTmpDir()
	utils.WarnNetworkScratch(scratch, "build")
	return scratch
}

// tmpRootForDef is the root for a definition build. Apptainer writes its root as
// a sandbox of many small files and builds it under --fakeroot, which NFS,
// Lustre, GPFS and PanFS do not support, so this is fast local scratch and the
// warning says the build will fail rather than merely run slowly.
func tmpRootForDef() string {
	scratch := utils.GetTmpDir()
	utils.WarnUnfakerootableScratch(scratch, "building a container")
	return scratch
}

// tmpRootForExternal picks the root for an external build (-f). Data keeps its
// large intermediates beside the target, whose location the user chose; an app
// and a definition go to fast local scratch, a definition because it must —
// see tmpRootForDef.
func tmpRootForExternal(targetDir string, typ catalog.Type, isDef bool) string {
	if isDef {
		return tmpRootForDef()
	}
	if typ == catalog.TypeData {
		return targetDir
	}
	scratch := utils.GetTmpDir()
	utils.WarnNetworkScratch(scratch, "build")
	return scratch
}
