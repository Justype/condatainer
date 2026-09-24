package cmd

import (
	"os"
	"testing"

	"github.com/Justype/condatainer/internal/testenv"
)

func TestMain(m *testing.M) {
	cleanup := testenv.Isolate()
	code := m.Run()
	cleanup()
	os.Exit(code)
}
