package utils

import (
	"os"
	"path/filepath"
	"reflect"
	"testing"
	"time"
)

func TestParseMemoryMB(t *testing.T) {
	tests := []struct {
		input  string
		wantMB int64
	}{
		{"8G", 8 * 1024},
		{"8GB", 8 * 1024},
		{"1024M", 1024},
		{"1024MB", 1024},
		{"4096K", 4},
		{"4096KB", 4},
		{"1T", 1024 * 1024},
		{"1TB", 1024 * 1024},
	}

	for _, tt := range tests {
		t.Run(tt.input, func(t *testing.T) {
			mb, err := ParseMemoryMB(tt.input)
			if err != nil {
				t.Errorf("ParseMemoryMB(%q) error: %v", tt.input, err)
				return
			}
			if mb != tt.wantMB {
				t.Errorf("ParseMemoryMB(%q) = %d MB; want %d MB", tt.input, mb, tt.wantMB)
			}
		})
	}
}

func TestParseWalltime(t *testing.T) {
	hour := time.Hour
	min := time.Minute
	sec := time.Second
	day := 24 * time.Hour

	tests := []struct {
		input   string
		want    time.Duration
		wantErr bool
	}{
		// Compound: Go-style with optional integer days
		{"4d12h", 4*day + 12*hour, false},
		{"2h30m", 2*hour + 30*min, false},
		{"3h", 3 * hour, false},
		{"3H", 3 * hour, false}, // case-insensitive
		{"90m", 90 * min, false},
		{"1.5h", 90 * min, false},
		{"1d2h30m45s", day + 2*hour + 30*min + 45*sec, false},
		{"4d", 4 * day, false},
		{"", 0, false},
		// Colon-separated
		{"01:30:00", hour + 30*min, false},
		{"1:30", hour + 30*min, false},
		{"02:30:00", 2*hour + 30*min, false},
		{"90", 90 * min, false}, // minutes only
		// D-HH:MM:SS
		{"1-12:00:00", day + 12*hour, false},
		{"2-06:00:00", 2*day + 6*hour, false},
		// Errors
		{"abc", 0, true},      // no valid unit letters
		{"1.5d", 0, true},     // fractional days not supported (integer only)
		{"bad:time", 0, true}, // letters in colon-separated path
	}
	for _, tt := range tests {
		t.Run(tt.input, func(t *testing.T) {
			dur, err := ParseWalltime(tt.input)
			if tt.wantErr {
				if err == nil {
					t.Errorf("ParseWalltime(%q): expected error, got %v", tt.input, dur)
				}
				return
			}
			if err != nil {
				t.Fatalf("ParseWalltime(%q) unexpected error: %v", tt.input, err)
			}
			if dur != tt.want {
				t.Errorf("ParseWalltime(%q) = %v; want %v", tt.input, dur, tt.want)
			}
		})
	}
}

func TestStripInlineComment(t *testing.T) {
	tests := []struct {
		name  string
		input string
		want  string
	}{
		{"No comment", "--cpus-per-task=8", "--cpus-per-task=8"},
		{"With comment", "--cpus-per-task=8  # This is a comment", "--cpus-per-task=8"},
		{"Comment only", "# Just a comment", ""},
		{"Multiple hashes", "--mem=16G # First # Second", "--mem=16G"},
		{"Hash in value needs escaping", "foo=bar#baz", "foo=bar"},
		{"Whitespace around comment", "--time=02:00:00   #   Time limit  ", "--time=02:00:00"},
		{"Empty after hash", "value #", "value"},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			got := StripInlineComment(tt.input)
			if got != tt.want {
				t.Errorf("StripInlineComment(%q) = %q, want %q", tt.input, got, tt.want)
			}
		})
	}
}

// Module lines name the site's own module tree, not CondaTainer artifacts, so
// they contribute nothing however closely they resemble a #DEP: name.
func TestGetDependenciesFromScriptIgnoresModuleLines(t *testing.T) {
	tmp, err := os.CreateTemp("", "script-deps-*.sh")
	if err != nil {
		t.Fatalf("failed to create temp file: %v", err)
	}
	defer os.Remove(tmp.Name())

	content := `#!/bin/bash
#DEP:foo/1.0
module load alpha/1.0 beta/2.0
ml load gamma/3.0 delta/4.0
ml purge
ml foo/5.0
`
	if _, err := tmp.WriteString(content); err != nil {
		t.Fatalf("failed to write to temp file: %v", err)
	}
	tmp.Close()

	deps, err := GetDependenciesFromScript(tmp.Name())
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}

	expected := []string{"foo/1.0"}
	if len(deps) != len(expected) {
		t.Fatalf("expected %d deps, got %d: %v", len(expected), len(deps), deps)
	}
	for i := range expected {
		if deps[i] != expected[i] {
			t.Fatalf("dep[%d] mismatch: expected %s, got %s", i, expected[i], deps[i])
		}
	}
}

func TestGetDependenciesFromScript_WithComments(t *testing.T) {
	tmp, err := os.CreateTemp("", "script-deps-comments-*.sh")
	if err != nil {
		t.Fatalf("failed to create temp file: %v", err)
	}
	defer os.Remove(tmp.Name())

	content := `#!/bin/bash
#DEP:foo/1.0  # This is a dependency comment
#DEP:bar/2.0 # Another comment
#DEP:baz/3.0#No space before hash
module load alpha/1.0 beta/2.0
`
	if _, err := tmp.WriteString(content); err != nil {
		t.Fatalf("failed to write to temp file: %v", err)
	}
	tmp.Close()

	deps, err := GetDependenciesFromScript(tmp.Name())
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}

	expected := []string{"foo/1.0", "bar/2.0", "baz/3.0"}
	if len(deps) != len(expected) {
		t.Fatalf("expected %d deps, got %d: %v", len(expected), len(deps), deps)
	}
	for i := range expected {
		if deps[i] != expected[i] {
			t.Fatalf("dep[%d] mismatch: expected %s, got %s", i, expected[i], deps[i])
		}
	}
}

func TestGetTypeFromScript(t *testing.T) {
	tests := []struct {
		name     string
		content  string
		wantType string
		wantErr  bool
	}{
		{
			name:     "NoTypeDefaultsToApp",
			content:  "#!/bin/bash\necho hello\n",
			wantType: "app",
		},
		{
			name:     "HashTypeData",
			content:  "#!/bin/bash\n#TYPE:data\n",
			wantType: "data",
		},
		{
			// A bare TYPE: has no '#', so it is not an annotation at all and the
			// script simply declares nothing.
			name:     "BareTypeIsNotAnAnnotation",
			content:  "#!/bin/bash\nTYPE:data\n",
			wantType: "app",
		},
		{
			name:    "InvalidType",
			content: "#!/bin/bash\n#TYPE:unknown\n",
			wantErr: true,
		},
		// Aliases are gone: #TYPE: has to mean the same thing here as it does to
		// catalog.DeriveType, which only ever accepted app and data.
		{
			name:    "RejectsFormerRefAlias",
			content: "#!/bin/bash\n#TYPE:ref\n",
			wantErr: true,
		},
		{
			name:    "RejectsFormerToolAlias",
			content: "#!/bin/bash\n#TYPE:tool\n",
			wantErr: true,
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			tmp, err := os.CreateTemp("", "script-type-*.sh")
			if err != nil {
				t.Fatalf("failed to create temp file: %v", err)
			}
			defer os.Remove(tmp.Name())

			if _, err := tmp.WriteString(tt.content); err != nil {
				t.Fatalf("failed to write temp file: %v", err)
			}
			tmp.Close()

			got, err := GetTypeFromScript(tmp.Name())
			if tt.wantErr {
				if err == nil {
					t.Fatalf("expected error, got nil")
				}
				return
			}
			if err != nil {
				t.Fatalf("unexpected error: %v", err)
			}
			if got != tt.wantType {
				t.Fatalf("GetTypeFromScript() = %q, want %q", got, tt.wantType)
			}
		})
	}
}

func TestSortVersionsDescending(t *testing.T) {
	tests := []struct {
		name  string
		input []string
		want  []string
	}{
		{"integers", []string{"22", "26", "23", "25", "24"}, []string{"26", "25", "24", "23", "22"}},
		{"semver", []string{"2.7.11a", "2.7.11b", "2.7.10"}, []string{"2.7.11b", "2.7.11a", "2.7.10"}},
		{"semver-mixed-suffix", []string{"2.7.9a", "2.7.11b", "2.7.10"}, []string{"2.7.11b", "2.7.10", "2.7.9a"}},
		{"mixed", []string{"101", "75", "151"}, []string{"151", "101", "75"}},
		{"single", []string{"1.0"}, []string{"1.0"}},
		{"empty", []string{}, []string{}},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			got := SortVersionsDescending(tt.input)
			if !reflect.DeepEqual(got, tt.want) {
				t.Errorf("SortVersionsDescending(%v) = %v; want %v", tt.input, got, tt.want)
			}
		})
	}
}

// Position carries no meaning: a declaration counts wherever it is written,
// including inside a heredoc that writes another script.
func TestGetDependenciesFromScriptReadsAnnotationsAnywhere(t *testing.T) {
	tmp, err := os.CreateTemp("", "script-deps-anywhere-*.sh")
	if err != nil {
		t.Fatalf("failed to create temp file: %v", err)
	}
	defer os.Remove(tmp.Name())

	content := "#!/bin/bash\n" +
		"#DEP:foo/1.0\n" +
		"echo running\n" +
		"#DEP:below/3.0\n" +
		"cat <<'END' > generated.sh\n" +
		"#DEP:heredoc/4.0\n" +
		"END\n"
	if _, err := tmp.WriteString(content); err != nil {
		t.Fatalf("failed to write to temp file: %v", err)
	}
	tmp.Close()

	deps, err := GetDependenciesFromScript(tmp.Name())
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}

	expected := []string{"foo/1.0", "below/3.0", "heredoc/4.0"}
	if len(deps) != len(expected) {
		t.Fatalf("expected %v, got %v", expected, deps)
	}
	for i := range expected {
		if deps[i] != expected[i] {
			t.Fatalf("dep[%d] = %s, want %s", i, deps[i], expected[i])
		}
	}
}

func TestGetTargetFromScript(t *testing.T) {
	tests := []struct {
		name    string
		content string
		want    string
		wantErr bool
	}{
		{
			name:    "NoTargetIsEmpty",
			content: "#!/bin/bash\necho hello\n",
			want:    "",
		},
		{
			name:    "DeclaredTarget",
			content: "#!/bin/bash\n#TARGET: star/2.7.11b/index\n",
			want:    "star/2.7.11b/index",
		},
		{
			// The same normalization every other name goes through, so a script
			// and a recipe cannot disagree about what they named.
			name:    "TargetIsNormalized",
			content: "#!/bin/bash\n#TARGET: /star/2.7.11b/index/\n",
			want:    "star/2.7.11b/index",
		},
		{
			// A bare TARGET: has no '#', so it is not an annotation at all.
			name:    "BareTargetIsNotAnAnnotation",
			content: "#!/bin/bash\nTARGET: star/2.7.11b\n",
			want:    "",
		},
		{
			// #TARGET: doubles as a template pattern for catalog recipes, where
			// #PH: supplies the values. An external build is addressed by path
			// and has neither, so a pattern could never be filled in.
			name:    "RejectsAPlaceholder",
			content: "#!/bin/bash\n#TARGET: star/{star_version}/index\n",
			wantErr: true,
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			path := filepath.Join(t.TempDir(), "build.sh")
			if err := os.WriteFile(path, []byte(tt.content), 0o644); err != nil {
				t.Fatalf("failed to write script: %v", err)
			}

			got, err := GetTargetFromScript(path)
			if tt.wantErr {
				if err == nil {
					t.Fatalf("expected error, got %q", got)
				}
				return
			}
			if err != nil {
				t.Fatalf("unexpected error: %v", err)
			}
			if got != tt.want {
				t.Fatalf("GetTargetFromScript() = %q, want %q", got, tt.want)
			}
		})
	}
}

// An empty component would make key.Role's component match silently fail,
// downgrading a real dependency to build history.
func TestGetTargetFromScriptRejectsAnEmptyComponent(t *testing.T) {
	path := filepath.Join(t.TempDir(), "build.sh")
	if err := os.WriteFile(path, []byte("#!/bin/bash\n#TARGET: star//index\n"), 0o644); err != nil {
		t.Fatal(err)
	}
	if got, err := GetTargetFromScript(path); err == nil {
		t.Fatalf("an empty component was accepted as %q", got)
	}
}

// A running script mounts what it names and records nothing, so an overlay path
// or an external .sqf is a legitimate declaration there. Only a build's #DEP:
// becomes an edge that has to mean the same thing on another machine, and that
// rule lives in catalog.ValidateDeps — never here.
func TestGetDependenciesFromScriptKeepsPathDeclarations(t *testing.T) {
	path := filepath.Join(t.TempDir(), "run.sh")
	content := "#!/bin/bash\n#DEP: samtools/1.21\n#DEP: env.img\n#DEP: env.ext3\n" +
		"#DEP: overlays/tool.sqf\n#DEP: /shared/vendor/tool.sqf\n"
	if err := os.WriteFile(path, []byte(content), 0o644); err != nil {
		t.Fatal(err)
	}

	deps, err := GetDependenciesFromScript(path)
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}
	want := []string{"samtools/1.21", "env.img", "env.ext3", "overlays/tool.sqf", "/shared/vendor/tool.sqf"}
	if !reflect.DeepEqual(deps, want) {
		t.Fatalf("deps = %v, want %v", deps, want)
	}
}

func TestFormatSpeed(t *testing.T) {
	tests := []struct {
		name    string
		bytes   int64
		elapsed time.Duration
		want    string
	}{
		{"steady", 100 << 20, time.Second, "100 MiB/s"},
		{"fraction of a second", 50 << 20, 500 * time.Millisecond, "100 MiB/s"},
		{"too short to quote", 100 << 20, 10 * time.Millisecond, ""},
		{"nothing moved", 0, time.Second, ""},
		{"a retried transfer withdrew bytes", -5, time.Second, ""},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			got, ok := FormatSpeed(tt.bytes, tt.elapsed)
			if got != tt.want || ok != (tt.want != "") {
				t.Errorf("FormatSpeed(%d, %v) = %q, %v; want %q", tt.bytes, tt.elapsed, got, ok, tt.want)
			}
		})
	}
}
