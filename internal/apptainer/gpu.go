package apptainer

import (
	"os"
	"os/exec"
	"strings"

	"github.com/Justype/condatainer/internal/utils"
)

// DetectGPUFlags inspects the host for GPU support and returns the necessary apptainer flags.
func DetectGPUFlags() []string {
	flags := []string{}
	if hasNvidiaGPU() {
		flags = append(flags, "--nv")
	}
	if hasRocmGPU() {
		flags = append(flags, "--rocm")
	}
	if len(flags) > 0 {
		utils.PrintDebug("Detected GPU flags: %s", strings.Join(flags, " "))
	}
	return flags
}

func hasNvidiaGPU() bool {
	if _, err := exec.LookPath("nvidia-smi"); err == nil {
		return true
	}
	entries, err := os.ReadDir("/dev")
	if err != nil {
		return false
	}
	for _, entry := range entries {
		if strings.Contains(strings.ToLower(entry.Name()), "nvidia") {
			return true
		}
	}
	return false
}

func hasRocmGPU() bool {
	if _, err := exec.LookPath("rocm-smi"); err == nil {
		return true
	}
	entries, err := os.ReadDir("/dev")
	if err != nil {
		return false
	}
	for _, entry := range entries {
		if strings.Contains(strings.ToLower(entry.Name()), "rocm") {
			return true
		}
	}
	return false
}
