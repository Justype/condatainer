package build

import (
	"os/exec"
	"testing"
)

// The toolchain is appended, so a tool the base or an overlay provides always
// wins, and the path is quoted for a directory a user named with a space or an
// apostrophe. Run through a real bash: quoting is exactly what a string
// comparison would let through.
func TestToolchainPathLineAppendsAndQuotes(t *testing.T) {
	bash, err := exec.LookPath("bash")
	if err != nil {
		t.Skip("bash not available")
	}
	const bin = "/opt/it's a dir/libexec/bin"

	cmd := exec.Command(bash, "-c", toolchainPathLine(bin)+`; printf %s "$PATH"`)
	cmd.Env = []string{"PATH=/usr/bin:/bin"}
	out, err := cmd.Output()
	if err != nil {
		t.Fatal(err)
	}
	if want := "/usr/bin:/bin:" + bin; string(out) != want {
		t.Errorf("PATH = %q, want %q", out, want)
	}
}

// Without a toolchain the wrapper script gains nothing.
func TestToolchainPathLineIsEmptyWithoutOne(t *testing.T) {
	if got := toolchainPathLine(""); got != "" {
		t.Errorf("toolchainPathLine(\"\") = %q, want nothing", got)
	}
}
