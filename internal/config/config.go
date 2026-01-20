package config

import (
	"os"
	"path/filepath"
	"runtime"
	"strconv"
)

const VERSION = "1.0.6"
const GitHubRepo = "Justype/condatainer"

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
	NCPUs        int
}

// Global holds the singleton configuration instance
var Global Config

func LoadDefaults(executablePath string) {
	programDir := filepath.Dir(executablePath)
	if absProgDir, err := filepath.Abs(programDir); err == nil {
		programDir = absProgDir
	}
	baseDir := filepath.Dir(programDir)
	baseDir = filepath.Clean(baseDir)

	// Developer Mode Check
	if _, err := os.Stat(filepath.Join(baseDir, "images")); os.IsNotExist(err) {
		cwd, err := os.Getwd()
		if err == nil {
			baseDir = cwd
			programDir = filepath.Join(cwd, "bin")
			if absProgDir, err := filepath.Abs(programDir); err == nil {
				programDir = absProgDir
			}
		}
	}

	if absBaseDir, err := filepath.Abs(baseDir); err == nil {
		baseDir = absBaseDir
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
		NCPUs:            computeNCPUs(),

		ApptainerBin: "apptainer",
		BaseImage:    filepath.Join(baseDir, "images", "base_image.sif"),
	}
}

func computeNCPUs() int {
	ncpus := 4
	if env := os.Getenv("SLURM_CPUS_PER_TASK"); env != "" {
		if parsed, err := strconv.Atoi(env); err == nil && parsed > 0 {
			ncpus = parsed
		}
	} else if runtime.NumCPU() > 0 {
		ncpus = runtime.NumCPU()
	}
	return ncpus
}
