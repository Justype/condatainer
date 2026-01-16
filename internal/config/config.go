package config

import (
	"os"
	"path/filepath"
)

const VERSION = "1.0.6"

// Config holds global application settings
type Config struct {
	Debug            bool
	SubmitJob        bool
	Version          string
	ProgramDir       string
	BaseDir          string
	ImagesDir        string
	RefImagesDir     string
	BuildScriptsDir  string
	HelperScriptsDir string
	TmpDir           string
	LogsDir          string
	CompressArgs     string
	TmpSizeMB        int

	ApptainerBin string
	BaseImage    string
}

// Global holds the singleton configuration instance
var Global Config

func LoadDefaults(executablePath string) {
	programDir := filepath.Dir(executablePath)
	baseDir := filepath.Dir(programDir)

	// Developer Mode Check
	if _, err := os.Stat(filepath.Join(baseDir, "images")); os.IsNotExist(err) {
		cwd, _ := os.Getwd()
		baseDir = cwd
		programDir = filepath.Join(cwd, "bin")
	}

	Global = Config{
		Debug:            false,
		SubmitJob:        true,
		Version:          VERSION,
		ProgramDir:       programDir,
		BaseDir:          baseDir,
		ImagesDir:        filepath.Join(baseDir, "images"),
		RefImagesDir:     filepath.Join(baseDir, "ref-images"),
		BuildScriptsDir:  filepath.Join(baseDir, "build-scripts"),
		HelperScriptsDir: filepath.Join(baseDir, "helper-scripts"),
		TmpDir:           filepath.Join(baseDir, "tmp"),
		LogsDir:          filepath.Join(os.Getenv("HOME"), "logs"),
		CompressArgs:     "-comp lz4", // zstd only compatible with apptainer version > 1.4
		TmpSizeMB:        20480,

		ApptainerBin: "apptainer",
		BaseImage:    filepath.Join(baseDir, "images", "base_image.sif"),
	}
}
