package apptainer

import (
	"os"
	"testing"
)

// TestMain keeps every test off the real per-user cache.
func TestMain(m *testing.M) {
	dir, err := os.MkdirTemp("", "cnt-apptainer-test")
	if err != nil {
		os.Exit(1)
	}
	os.Setenv("XDG_CACHE_HOME", dir)
	code := m.Run()
	os.RemoveAll(dir)
	os.Exit(code)
}
