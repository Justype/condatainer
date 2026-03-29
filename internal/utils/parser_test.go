package utils

import (
	"os"
	"reflect"
	"sort"
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

func TestGetInteractivePromptsFromScript(t *testing.T) {
	tmp, err := os.CreateTemp("", "script-*.sh")
	if err != nil {
		t.Fatalf("failed to create temp file: %v", err)
	}
	defer os.Remove(tmp.Name())

	content := `#!/bin/bash
#DEP:foo/1.0
#INTERACTIVE:Please enter the license key
#INTERACTIVE:Second prompt text
echo hello
`
	if _, err := tmp.WriteString(content); err != nil {
		t.Fatalf("failed to write to temp file: %v", err)
	}
	tmp.Close()

	prompts, err := GetInteractivePromptsFromScript(tmp.Name())
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}
	if len(prompts) != 2 {
		t.Fatalf("expected 2 prompts, got %d", len(prompts))
	}
	if prompts[0] != "Please enter the license key" {
		t.Fatalf("unexpected first prompt: %s", prompts[0])
	}
	if prompts[1] != "Second prompt text" {
		t.Fatalf("unexpected second prompt: %s", prompts[1])
	}
}

func TestGetInteractivePromptsFromScript_Missing(t *testing.T) {
	_, err := GetInteractivePromptsFromScript("/nonexistent/path/to/script.sh")
	if err == nil {
		t.Fatalf("expected error for missing file")
	}
}

func TestGetDependenciesFromScript_MLAndModuleLoad(t *testing.T) {
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

	deps, err := GetDependenciesFromScript(tmp.Name(), true)
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}

	expected := []string{"foo/1.0", "alpha/1.0", "beta/2.0", "gamma/3.0", "delta/4.0", "foo/5.0"}
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

	deps, err := GetDependenciesFromScript(tmp.Name(), true)
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}

	expected := []string{"foo/1.0", "bar/2.0", "baz/3.0", "alpha/1.0", "beta/2.0"}
	if len(deps) != len(expected) {
		t.Fatalf("expected %d deps, got %d: %v", len(expected), len(deps), deps)
	}
	for i := range expected {
		if deps[i] != expected[i] {
			t.Fatalf("dep[%d] mismatch: expected %s, got %s", i, expected[i], deps[i])
		}
	}
}

func TestGetExternalBuildTypeFromScript(t *testing.T) {
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
			name:     "TypeRefAlias",
			content:  "#!/bin/bash\nTYPE:ref\n",
			wantType: "data",
		},
		{
			name:     "HashTypeToolAlias",
			content:  "#!/bin/bash\n#TYPE:tool\n",
			wantType: "app",
		},
		{
			name:    "InvalidType",
			content: "#!/bin/bash\n#TYPE:unknown\n",
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

			got, err := GetExternalBuildTypeFromScript(tmp.Name())
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
				t.Fatalf("GetExternalBuildTypeFromScript() = %q, want %q", got, tt.wantType)
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

func TestParsePLValues(t *testing.T) {
	tests := []struct {
		name      string
		raw       string
		wantVals  []string
		wantOpen  bool
		wantError bool
	}{
		{
			name:     "comma list",
			raw:      "2.7.11b,2.7.11a",
			wantVals: []string{"2.7.11b", "2.7.11a"},
			wantOpen: false,
		},
		{
			name:     "integer range",
			raw:      "22-25",
			wantVals: []string{"25", "24", "23", "22"},
			wantOpen: false,
		},
		{
			name:     "open ended",
			raw:      "75,101,151,*",
			wantVals: []string{"151", "101", "75", "*"},
			wantOpen: true,
		},
		{
			name:     "star only",
			raw:      "*",
			wantVals: []string{"*"},
			wantOpen: true,
		},
		{
			name:      "empty",
			raw:       "",
			wantError: true,
		},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			got, open, err := parsePLValues(tt.raw)
			if tt.wantError {
				if err == nil {
					t.Error("expected error, got nil")
				}
				return
			}
			if err != nil {
				t.Fatalf("unexpected error: %v", err)
			}
			if !reflect.DeepEqual(got, tt.wantVals) {
				t.Errorf("parsePLValues(%q) values = %v; want %v", tt.raw, got, tt.wantVals)
			}
			if open != tt.wantOpen {
				t.Errorf("parsePLValues(%q) open = %v; want %v", tt.raw, open, tt.wantOpen)
			}
		})
	}
}

func TestGetPlaceholdersFromScript(t *testing.T) {
	f, err := os.CreateTemp("", "pltest*.sh")
	if err != nil {
		t.Fatal(err)
	}
	defer os.Remove(f.Name())

	_, _ = f.WriteString(`#!/bin/bash
#PL:star_version:2.7.11b,2.7.11a
#PL:gencode_version:22-24
#PL:read_length:75,101,151,*
#TARGET:grch38/star/{star_version}/gencode{gencode_version}-{read_length}
`)
	f.Close()

	defs, err := GetPlaceholdersFromScript(f.Name())
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}
	if len(defs) != 3 {
		t.Fatalf("expected 3 defs, got %d", len(defs))
	}

	// star_version: sorted desc
	if defs[0].Name != "star_version" {
		t.Errorf("defs[0].Name = %q; want star_version", defs[0].Name)
	}
	if defs[0].Values[0] != "2.7.11b" {
		t.Errorf("defs[0].Values[0] = %q; want 2.7.11b", defs[0].Values[0])
	}

	// gencode_version: range 22-24 sorted desc
	if defs[1].Name != "gencode_version" {
		t.Errorf("defs[1].Name = %q; want gencode_version", defs[1].Name)
	}
	if len(defs[1].Values) != 3 || defs[1].Values[0] != "24" {
		t.Errorf("gencode_version values = %v; want [24 23 22]", defs[1].Values)
	}

	// read_length: open-ended
	if !defs[2].Open {
		t.Error("read_length should be open")
	}
	if defs[2].Values[len(defs[2].Values)-1] != "*" {
		t.Error("last value of open-ended read_length should be *")
	}
}

func TestGetTargetFromScript(t *testing.T) {
	f, err := os.CreateTemp("", "targettest*.sh")
	if err != nil {
		t.Fatal(err)
	}
	defer os.Remove(f.Name())
	_, _ = f.WriteString("#TARGET:grch38/star/{star_version}/gencode{gencode_version}\n")
	f.Close()

	got, err := GetTargetFromScript(f.Name())
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}
	want := "grch38/star/{star_version}/gencode{gencode_version}"
	if got != want {
		t.Errorf("GetTargetFromScript() = %q; want %q", got, want)
	}
}

func TestInterpolateVars(t *testing.T) {
	vars := map[string]string{"star_version": "2.7.11b", "gencode_version": "22"}
	got := InterpolateVars("grch38/star/{star_version}/gencode{gencode_version}-101", vars)
	want := "grch38/star/2.7.11b/gencode22-101"
	if got != want {
		t.Errorf("InterpolateVars() = %q; want %q", got, want)
	}

	// Unknown keys left unchanged
	got2 := InterpolateVars("{unknown}/foo", vars)
	if got2 != "{unknown}/foo" {
		t.Errorf("InterpolateVars unknown key: got %q; want {unknown}/foo", got2)
	}
}

func TestExpandPlaceholders(t *testing.T) {
	defs := []PlaceholderDef{
		{Name: "a", Values: []string{"1", "2"}},
		{Name: "b", Values: []string{"x", "y"}},
	}
	combos := ExpandPlaceholders(defs)
	if len(combos) != 4 {
		t.Fatalf("expected 4 combinations, got %d", len(combos))
	}

	// Collect all unique a,b pairs to verify completeness
	type pair struct{ a, b string }
	got := make(map[pair]bool)
	for _, c := range combos {
		got[pair{c["a"], c["b"]}] = true
	}
	for _, a := range []string{"1", "2"} {
		for _, b := range []string{"x", "y"} {
			if !got[pair{a, b}] {
				t.Errorf("missing combination a=%s b=%s", a, b)
			}
		}
	}
}

func TestExpandPlaceholders_StarSkipped(t *testing.T) {
	defs := []PlaceholderDef{
		{Name: "ver", Values: []string{"1.0", "2.0", "*"}, Open: true},
	}
	combos := ExpandPlaceholders(defs)
	// Should only produce 2 combos (not 3 — "*" is excluded)
	if len(combos) != 2 {
		t.Fatalf("expected 2 combinations, got %d: %v", len(combos), combos)
	}
}

func TestExpandPlaceholders_Empty(t *testing.T) {
	combos := ExpandPlaceholders(nil)
	if len(combos) != 1 {
		t.Fatalf("expected [{}], got %v", combos)
	}
}

func TestDefaultVarsFromPlaceholders(t *testing.T) {
	defs := []PlaceholderDef{
		{Name: "a", Values: []string{"3", "2", "1"}},
		{Name: "b", Values: []string{"x", "y"}},
		{Name: "c", Values: []string{"*"}, Open: true}, // open-only
	}
	vars := DefaultVarsFromPlaceholders(defs)
	if vars["a"] != "3" {
		t.Errorf("default for a = %q; want 3", vars["a"])
	}
	if vars["b"] != "x" {
		t.Errorf("default for b = %q; want x", vars["b"])
	}
	if _, ok := vars["c"]; ok {
		t.Error("open-only placeholder c should not appear in DefaultVars")
	}
}

// helper to ensure a slice of strings is sorted descending (ignores order of equal elements)
func isSortedDesc(vals []string) bool {
	cp := make([]string, len(vals))
	copy(cp, vals)
	sort.Slice(cp, func(i, j int) bool { return naturalVersionGreater(cp[i], cp[j]) })
	return reflect.DeepEqual(vals, cp)
}
