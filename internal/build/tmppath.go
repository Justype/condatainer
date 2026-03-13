package build

import (
	"path/filepath"
	"strings"

	"github.com/Justype/condatainer/internal/config"
	"github.com/Justype/condatainer/internal/utils"
)

// resolveTmpDirForConda returns the fast local node scratch for conda/app builds.
// Prefers scheduler-assigned scratch (SLURM_TMPDIR, etc.), then TMPDIR, then /tmp/cnt-$USER.
func resolveTmpDirForConda() string {
	return utils.GetTmpDir()
}

// resolveTmpDirForRef returns the stable condatainer writable tmp for ref/data builds.
// Ref builds involve large datasets that benefit from a persistent, shared filesystem path.
func resolveTmpDirForRef() string {
	return config.GetWritableTmpDir()
}

// resolveTmpDirForDef returns the stable condatainer writable tmp for def/base-image builds.
// Def builds produce a SIF which is also stored here before extraction.
func resolveTmpDirForDef() string {
	return config.GetWritableTmpDir()
}

// resolveTmpDirForExternal returns the directory next to the target for external source builds.
// This keeps build artifacts alongside the user-specified output location.
func resolveTmpDirForExternal(targetDir string) string {
	return targetDir
}

// buildTmpPaths computes tmpOverlayPath and cntDirPath from tmpDir and nameVersion.
// ext is ".img" for conda/script builds, ".sif" for def builds, "" for dir-mode (tmpOverlay will be "").
func buildTmpPaths(nameVersion, tmpDir, ext string) (tmpOverlayPath, cntDirPath string) {
	cntDirPath = getCntDirPath(nameVersion, tmpDir)
	if ext != "" {
		filename := strings.ReplaceAll(nameVersion, "/", "--") + ext
		tmpOverlayPath = filepath.Join(tmpDir, filename)
	}
	return
}

// getCntDirPath returns the container directory path for a name/version.
// Format: <tmpDir>/build_<nameVersion>/cnt
func getCntDirPath(nameVersion, tmpDir string) string {
	buildDirName := "build_" + strings.ReplaceAll(nameVersion, "/", "_")
	return filepath.Join(tmpDir, buildDirName, "cnt")
}

// getBuildTmpDir returns the host tmp directory used for TMPDIR/micromamba root in dir-mode builds.
// Format: <tmpDir>/build_<nameVersion>/tmp
func getBuildTmpDir(b *BaseBuildObject) string {
	return filepath.Join(filepath.Dir(b.cntDirPath), "tmp")
}

// getTmpOverlayPath returns the temporary overlay path for ext3 builds (script/conda).
// Format: <tmpDir>/<nameVersion>.img (with / replaced by --)
// Kept for callers that need the .img path explicitly (e.g. NewBaseImageBuildObject).
func getTmpOverlayPath(nameVersion, tmpDir string) string {
	filename := strings.ReplaceAll(nameVersion, "/", "--") + ".img"
	return filepath.Join(tmpDir, filename)
}
