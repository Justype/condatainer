package container

import (
	"log/slog"
	"os"
	"strings"

	"github.com/Justype/condatainer/internal/config"
)

// DetectGPUFlags inspects the host for GPU support and returns the necessary apptainer flags.
//
// Detection is a stat of the device node, so it reports that a driver is
// installed, not that it works. On a node where the driver is loaded but
// unusable, the node is still there and the container then fails to start —
// which is what `autoload_gpu: false` is for.
func DetectGPUFlags() []string {
	if !config.Global.AutoloadGPU {
		slog.Default().Debug("GPU autoload disabled, passing no GPU flags")
		return nil
	}

	flags := []string{}
	if hasNvidiaGPU() {
		flags = append(flags, "--nv")
	}
	if hasRocmGPU() {
		flags = append(flags, "--rocm")
	}
	if len(flags) > 0 {
		slog.Default().Debug("detected GPU flags", "flags", strings.Join(flags, " "))
	}
	return flags
}

// hasNvidiaGPU checks for the presence of NVIDIA GPUs on the host by /dev/nvidiactl.
// nvidiactl is always present when the NVIDIA driver is loaded, regardless of GPU
// numbering or MIG configuration (unlike /dev/nvidia0 which may be absent on MIG nodes
// or when the allocated GPU index is not 0).
func hasNvidiaGPU() bool {
	_, err := os.Stat("/dev/nvidiactl")
	return err == nil
}

// hasRocmGPU checks for the presence of ROCm GPUs on the host by /dev/kfd.
func hasRocmGPU() bool {
	_, err := os.Stat("/dev/kfd")
	return err == nil
}
