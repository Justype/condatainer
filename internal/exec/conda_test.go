package exec

import (
	"context"
	"testing"
)

func TestInitCondaEnvSkipsEmptyPackageList(t *testing.T) {
	err := InitCondaEnv(context.Background(), "/path/that/does/not/exist.img", nil, false, IO{})
	if err != nil {
		t.Fatalf("empty initialization returned %v", err)
	}
}
