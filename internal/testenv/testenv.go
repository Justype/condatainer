// Package testenv isolates a test binary from the CondaTainer installation and
// configuration of the machine it runs on.
package testenv

import (
	"os"
	"strings"
)

// Isolate points HOME and every data root at a fresh temporary directory,
// clears every CNT_* variable, and fixes default_distro, which every command
// requires. It returns a func that removes the directory.
//
// Call it from TestMain before any test runs: config caches the installation
// root on first use, so a later t.Setenv has no effect.
func Isolate() (cleanup func()) {
	dir, err := os.MkdirTemp("", "condatainer-test-")
	if err != nil {
		panic(err)
	}
	for _, kv := range os.Environ() {
		if name, _, _ := strings.Cut(kv, "="); strings.HasPrefix(name, "CNT_") {
			os.Unsetenv(name)
		}
	}
	for name, value := range map[string]string{
		"HOME":               dir,
		"CNT_ROOT":           dir,
		"CNT_DEFAULT_DISTRO": "ubuntu24",
		"SCRATCH":            "",
		"XDG_DATA_HOME":      "",
		"XDG_CONFIG_HOME":    "",
	} {
		os.Setenv(name, value)
	}
	return func() { os.RemoveAll(dir) }
}
