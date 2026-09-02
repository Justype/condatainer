package utils

import "testing"

// The mapping is the whole point: a magic number nobody recognises must not be
// reported as a network filesystem, since the warning it drives should never
// fire on a guess.
func TestFSKindNamesAndClassifies(t *testing.T) {
	for _, tc := range []struct {
		magic   int64
		name    string
		network bool
	}{
		{fsNFS, "nfs", true},
		{fsLustre, "lustre", true},
		{fsGPFS, "gpfs", true},
		{fsXFS, "xfs", false},
		{fsTmpfs, "tmpfs", false},
		{fsExt, "ext", false},
	} {
		if got := fsNames[tc.magic]; got != tc.name {
			t.Errorf("fsNames[%#x] = %q, want %q", tc.magic, got, tc.name)
		}
		if got := networkFS[tc.magic]; got != tc.network {
			t.Errorf("networkFS[%#x] = %v, want %v", tc.magic, got, tc.network)
		}
	}
	if networkFS[0x12345678] {
		t.Error("an unrecognised filesystem was classified as network")
	}
}

// An unreadable path reports nothing rather than failing: FSKind only ever
// improves a message, so it must never become a reason a command errors.
func TestFSKindOnMissingPathIsQuiet(t *testing.T) {
	name, network := FSKind("/nonexistent/path/for/this/test")
	if name != "" || network {
		t.Errorf("FSKind on a missing path = (%q, %v), want (\"\", false)", name, network)
	}
}

// A real directory reports something, whatever this machine happens to run.
func TestFSKindOnARealDirectory(t *testing.T) {
	if name, _ := FSKind(t.TempDir()); name == "" {
		t.Error("FSKind reported nothing for a directory that exists")
	}
}
