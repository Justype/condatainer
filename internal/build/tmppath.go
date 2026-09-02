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

// tmpRootForDef is the root for a definition build, which keeps its recipe and
// the multi-GB SIF apptainer writes in one place.
func tmpRootForDef() string {
	return config.GetWritableTmpDir()
}

// tmpRootForExternal picks the root for an external build (-f). An app goes to
// fast local scratch; data and definitions keep their large intermediates beside
// the target, whose location the user chose.
func tmpRootForExternal(targetDir string, typ catalog.Type, isDef bool) string {
	if isDef || typ == catalog.TypeData {
		return targetDir
	}
	scratch := utils.GetTmpDir()
	utils.WarnNetworkScratch(scratch, "build")
	return scratch
}
