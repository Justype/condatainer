package cmd

import (
	"os"
	"path/filepath"
	"testing"
)

// A frozen environment is always called env, so position is the only question
// left: where the artifact lands. These are the forms a user can type.
func TestResolveFreezeTarget(t *testing.T) {
	dir := t.TempDir()
	source := filepath.Join(dir, "env.img")
	existingDir := filepath.Join(dir, "overlays")
	if err := os.MkdirAll(existingDir, 0o755); err != nil {
		t.Fatal(err)
	}

	cases := []struct {
		label    string
		dest     string
		wantFile string
	}{
		{"nothing given: beside the source, extension changed",
			"", filepath.Join(dir, "env.sqf")},
		{"a destination places it",
			filepath.Join(dir, "frozen.sqf"), filepath.Join(dir, "frozen.sqf")},
		{"a destination without the extension",
			filepath.Join(dir, "frozen"), filepath.Join(dir, "frozen.sqf")},
		{"an existing directory is a place, not the artifact",
			existingDir, filepath.Join(existingDir, "env.sqf")},
	}
	for _, tc := range cases {
		t.Run(tc.label, func(t *testing.T) {
			target, err := resolveFreezeTarget(source, tc.dest)
			if err != nil {
				t.Fatalf("unexpected error: %v", err)
			}
			if target != tc.wantFile {
				t.Errorf("target = %q, want %q", target, tc.wantFile)
			}
		})
	}
}

// A trailing separator says "put it in this directory", and it means that
// whether or not the directory is there yet. Without that reading, a missing one
// fell through to the extension branch and produced a file named ".sqf".
func TestResolveFreezeTargetTrailingSlashIsADirectory(t *testing.T) {
	dir := t.TempDir()
	source := filepath.Join(dir, "env.img")
	missing := filepath.Join(dir, "nodir") + string(filepath.Separator)

	if _, err := resolveFreezeTarget(source, missing); err == nil {
		t.Fatal("a destination in a directory that does not exist was accepted")
	}
	if err := os.MkdirAll(filepath.Join(dir, "nodir"), 0o755); err != nil {
		t.Fatal(err)
	}
	got, err := resolveFreezeTarget(source, missing)
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}
	if want := filepath.Join(dir, "nodir", "env.sqf"); got != want {
		t.Errorf("target = %q, want %q", got, want)
	}
}

// A destination whose directory is missing is refused here, with a sentence,
// rather than by Apptainer failing to bind it after the payload has been walked.
func TestResolveFreezeTargetRefusesAMissingDirectory(t *testing.T) {
	dir := t.TempDir()
	source := filepath.Join(dir, "env.img")
	if _, err := resolveFreezeTarget(source, filepath.Join(dir, "nope", "env.sqf")); err == nil {
		t.Fatal("a destination under a missing directory was accepted")
	}
}

// Unfreeze infers its destination the way freeze does: beside the artifact with
// the extension changed.
func TestResolveUnfreezeTarget(t *testing.T) {
	dir := t.TempDir()
	artifact := filepath.Join(dir, "env.sqf")
	if err := os.WriteFile(artifact, []byte("x"), 0o644); err != nil {
		t.Fatal(err)
	}

	got, err := resolveUnfreezeTarget(artifact, "")
	if err != nil {
		t.Fatalf("inferred destination: %v", err)
	}
	if want := filepath.Join(dir, "env.img"); got != want {
		t.Errorf("got %q, want %q", got, want)
	}

	// A directory is a place to put it, not the name of it.
	got, err = resolveUnfreezeTarget(artifact, dir)
	if err != nil {
		t.Fatalf("directory destination: %v", err)
	}
	if want := filepath.Join(dir, "env.img"); got != want {
		t.Errorf("directory destination = %q, want %q", got, want)
	}

	// An overlay already there is what someone is developing in.
	occupied := filepath.Join(dir, "env.img")
	if err := os.WriteFile(occupied, []byte("x"), 0o644); err != nil {
		t.Fatal(err)
	}
	if _, err := resolveUnfreezeTarget(artifact, ""); err == nil {
		t.Error("an existing overlay was not refused")
	}

	// Only an .img: the result is writable and nothing else can hold it.
	if _, err := resolveUnfreezeTarget(artifact, filepath.Join(dir, "dev.sqf")); err == nil {
		t.Error("a .sqf destination was accepted")
	}
}
