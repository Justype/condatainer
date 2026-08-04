package cmd

import (
	"context"
	"testing"

	execpkg "github.com/Justype/condatainer/internal/exec"
)

func TestInitCondaInOverlaySkipsBlankOverlay(t *testing.T) {
	err := initCondaInOverlay(context.Background(), "/path/that/does/not/exist.img", "", nil, false, execpkg.IO{})
	if err != nil {
		t.Fatalf("blank overlay initialization returned %v", err)
	}
}
