package build

import (
	"os"
	"path/filepath"
	"strings"

	"github.com/Justype/condatainer/internal/config"
	"github.com/Justype/condatainer/internal/utils"
)

// resolveTmpDirForConda returns the tmp directory for conda/app builds.
// If CNT_TMPDIR is set, it takes precedence over scheduler-assigned scratch and TMPDIR.
func resolveTmpDirForConda() string {
	return utils.GetTmpDir()
}

// resolveTmpDirForRef returns the writable tmp for ref/data builds.
// Ref builds involve large datasets that benefit from a persistent, shared filesystem path.
// If CNT_TMPDIR is set, it overrides the writable tmp path.
func resolveTmpDirForRef() string {
	return config.GetWritableTmpDir()
}

// resolveTmpDirForDef returns the writable tmp for def/base-image builds.
// Def builds produce a SIF which is also stored here before extraction.
// If CNT_TMPDIR is set, it overrides the writable tmp path.
func resolveTmpDirForDef() string {
	return config.GetWritableTmpDir()
}

// resolveTmpDirForExternal resolves tmp directory for external source builds by TYPE.
// CNT_TMPDIR has highest priority and overrides all external TYPE behaviors.
// Without CNT_TMPDIR: TYPE=app (default) uses dynamic scratch (utils.GetTmpDir), TYPE=data uses target-adjacent path.
func resolveTmpDirForExternal(targetDir, externalType string) string {
	if os.Getenv("CNT_TMPDIR") != "" {
		return utils.GetTmpDir()
	}
	if strings.ToLower(strings.TrimSpace(externalType)) == "data" {
		return targetDir
	}
	return utils.GetTmpDir()
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
