package build

import (
	"os"
	"path/filepath"
	"strings"

	"github.com/Justype/condatainer/catalog"
	"github.com/Justype/condatainer/internal/config"
	"github.com/Justype/condatainer/internal/utils"
)

// A build works in fast local scratch (utils.GetTmpDir) or in the stable
// writable tmp (config.GetWritableTmpDir). See the README's Workspace Strategy.

// tmpRootForType picks where a build does its work: the stable writable tmp for
// data, fast local scratch for everything else.
func tmpRootForType(typ catalog.Type) string {
	if typ == catalog.TypeData {
		return config.GetWritableTmpDir()
	}
	return utils.GetTmpDir()
}

// resolveTmpDirForDef returns the tmp directory for definition builds, which
// keep their recipe and the SIF apptainer writes in one place.
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

// getCntDirPath returns the container directory path for a name/version.
// Format: <tmpDir>/build_<nameVersion>/cnt
func getCntDirPath(nameVersion, tmpDir string) string {
	buildDirName := "build_" + strings.ReplaceAll(nameVersion, "/", "_")
	return filepath.Join(tmpDir, buildDirName, "cnt")
}
