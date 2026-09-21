package registry

import (
	"os"
	"path/filepath"
	"strings"
	"testing"
)

// credFile makes a credential file with the given file and directory modes.
func credFile(t *testing.T, fileMode, dirMode os.FileMode) string {
	t.Helper()
	dir := filepath.Join(t.TempDir(), "layer")
	if err := os.Mkdir(dir, 0o700); err != nil {
		t.Fatal(err)
	}
	path := filepath.Join(dir, credentialFileName)
	if err := os.WriteFile(path, []byte("{}"), 0o600); err != nil {
		t.Fatal(err)
	}
	if err := os.Chmod(path, fileMode); err != nil {
		t.Fatal(err)
	}
	if err := os.Chmod(dir, dirMode); err != nil {
		t.Fatal(err)
	}
	t.Cleanup(func() { os.Chmod(dir, 0o700) })
	return path
}

func TestAccessOfFollowsTheFileAndItsDirectory(t *testing.T) {
	tests := []struct {
		name      string
		file, dir os.FileMode
		readers   Reach
		writers   Reach
	}{
		{"private", 0o600, 0o700, ReachOwner, ReachOwner},
		{"group read", 0o640, 0o750, ReachGroup, ReachOwner},
		{"group read and write", 0o660, 0o770, ReachGroup, ReachGroup},
		{"world readable", 0o644, 0o755, ReachEveryone, ReachOwner},
		{"world writable file", 0o666, 0o755, ReachEveryone, ReachEveryone},
		// The directory decides who reaches the file at all.
		{"readable file in a closed directory", 0o644, 0o700, ReachOwner, ReachOwner},
		// Anyone who can write the directory can replace the file, whatever its mode.
		{"private file in a writable directory", 0o600, 0o777, ReachOwner, ReachEveryone},
		// A sticky directory refuses deleting another user's file.
		{"private file in a sticky directory", 0o600, 0o777 | os.ModeSticky, ReachOwner, ReachOwner},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			got, err := accessOf(credFile(t, tt.file, tt.dir))
			if err != nil {
				t.Fatal(err)
			}
			if got.Readers != tt.readers || got.Writers != tt.writers {
				t.Errorf("readers=%d writers=%d, want %d and %d", got.Readers, got.Writers, tt.readers, tt.writers)
			}
		})
	}
}

func TestFindingsWarnAboutAPersonalFileOthersCanRead(t *testing.T) {
	path := credFile(t, 0o644, 0o755)
	findings := Findings("user", path)
	if len(findings) != 1 || !findings[0].Warn || !strings.Contains(findings[0].Text, "chmod 600") {
		t.Fatalf("findings = %+v", findings)
	}
	if got := Findings("user", credFile(t, 0o600, 0o700)); len(got) != 0 {
		t.Errorf("a private file was reported: %+v", got)
	}
}

// A shared layer is meant to be read by others, so reading alone is not a finding.
func TestFindingsLeaveASharedFileReadableByTheGroupAlone(t *testing.T) {
	if got := Findings("extra-root", credFile(t, 0o640, 0o750)); len(got) != 0 {
		t.Errorf("findings = %+v", got)
	}
}

func TestFindingsReportWhoCanReplaceTheCredential(t *testing.T) {
	group := Findings("extra-root", credFile(t, 0o660, 0o770))
	if len(group) != 1 || group[0].Warn {
		t.Errorf("group replace in a shared layer should be a note: %+v", group)
	}
	if personal := Findings("user", credFile(t, 0o600, 0o770)); len(personal) != 1 || !personal[0].Warn {
		t.Errorf("group replace in the user layer should warn: %+v", personal)
	}
	for _, layer := range []string{"user", "extra-root", "app-root"} {
		anyone := Findings(layer, credFile(t, 0o666, 0o777))
		warned := false
		for _, finding := range anyone {
			warned = warned || finding.Warn && strings.Contains(finding.Text, "anyone can replace")
		}
		if !warned {
			t.Errorf("%s: no warning that anyone can replace it: %+v", layer, anyone)
		}
	}
}

func TestReadableByNamesTheAudience(t *testing.T) {
	if got := ReadableBy(credFile(t, 0o600, 0o700)); got != "you" {
		t.Errorf("private = %q", got)
	}
	if got := ReadableBy(credFile(t, 0o644, 0o755)); got != "everyone" {
		t.Errorf("world = %q", got)
	}
	if got := ReadableBy(credFile(t, 0o640, 0o750)); !strings.HasPrefix(got, "group ") {
		t.Errorf("group = %q", got)
	}
	if got := ReadableBy(filepath.Join(t.TempDir(), "missing")); got != "" {
		t.Errorf("missing = %q", got)
	}
}
