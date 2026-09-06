package exec

import (
	"context"
	"slices"
	"testing"
)

func TestInitCondaEnvSkipsEmptyPackageList(t *testing.T) {
	err := InitCondaEnv(context.Background(), "/path/that/does/not/exist.img", "", nil, false, IO{})
	if err != nil {
		t.Fatalf("empty initialization returned %v", err)
	}
}

func TestCondaScratchOverlaysAppendsSnapshotOnlyWhenGiven(t *testing.T) {
	if got := condaScratchOverlays("/tmp/x.img", ""); !slices.Equal(got, []string{"/tmp/x.img"}) {
		t.Fatalf("overlays = %v, want just the img", got)
	}
	want := []string{"/tmp/x.img", "/proj/env.sqf"}
	if got := condaScratchOverlays("/tmp/x.img", "/proj/env.sqf"); !slices.Equal(got, want) {
		t.Fatalf("overlays = %v, want %v", got, want)
	}
}
