package container

import (
	"testing"

	"github.com/Justype/condatainer/internal/config"
)

// autoload_gpu: false is the escape hatch for a node whose driver is installed
// but unusable — detection fires on the device node and the container then
// fails to start, with nothing else able to stop it.
func TestDetectGPUFlagsRespectsAutoload(t *testing.T) {
	prev := config.Global.AutoloadGPU
	t.Cleanup(func() { config.Global.AutoloadGPU = prev })

	config.Global.AutoloadGPU = false
	if flags := DetectGPUFlags(); len(flags) != 0 {
		t.Errorf("DetectGPUFlags() = %v with autoload off, want none", flags)
	}

	// With it on, the result is whatever this host has — the only invariant is
	// that it stops suppressing.
	config.Global.AutoloadGPU = true
	onFlags := DetectGPUFlags()
	if hasNvidiaGPU() && !contains(onFlags, "--nv") {
		t.Errorf("DetectGPUFlags() = %v, want --nv on a host with /dev/nvidiactl", onFlags)
	}
	if !hasNvidiaGPU() && !hasRocmGPU() && len(onFlags) != 0 {
		t.Errorf("DetectGPUFlags() = %v on a host with no GPU device nodes", onFlags)
	}
}

func contains(list []string, want string) bool {
	for _, v := range list {
		if v == want {
			return true
		}
	}
	return false
}

// The default has to stay on: turning it off by accident silently drops GPU
// access from every container.
func TestAutoloadGPUDefaultsOn(t *testing.T) {
	config.LoadDefaults("/usr/local/bin/condatainer")
	if !config.Global.AutoloadGPU {
		t.Error("autoload_gpu defaulted to false")
	}
}
