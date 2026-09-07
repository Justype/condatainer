package exec

import (
	"os/exec"
	"testing"
)

func TestWrapWithActivationRunsActivationThenExec(t *testing.T) {
	activation := "export CNT_TEST_ACTIVATED=1\n"
	command := []string{"bash", "-c", "echo activated=$CNT_TEST_ACTIVATED"}

	wrapped := wrapWithActivation(activation, command)
	out, err := exec.Command(wrapped[0], wrapped[1:]...).CombinedOutput()
	if err != nil {
		t.Fatalf("wrapped command: %v\n%s", err, out)
	}
	if got, want := string(out), "activated=1\n"; got != want {
		t.Errorf("output = %q, want %q", got, want)
	}
}
