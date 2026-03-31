package build

import (
	"reflect"
	"strings"
	"testing"
)

// ---------------------------------------------------------------------------
// expandTemplates
// ---------------------------------------------------------------------------

func TestExpandTemplates_DepsInterpolated(t *testing.T) {
	scripts := map[string]ScriptInfo{
		"grch38/star-gencode": {
			IsTemplate:     true,
			TargetTemplate: "grch38/star/{star_version}/gencode{gencode_version}-{read_length}",
			PLOrder:        []string{"star_version", "gencode_version", "read_length"},
			Deps:           []string{"grch38/genome/gencode", "grch38/gtf-gencode/{gencode_version}", "star/{star_version}"},
			Whatis:         "STAR index GENCODE{gencode_version} read length {read_length}",
			PL: map[string][]string{
				"star_version":    {"2.7.11b"},
				"gencode_version": {"47"},
				"read_length":     {"101"},
			},
		},
	}

	expandTemplates(scripts)

	entry, ok := scripts["grch38/star/2.7.11b/gencode47-101"]
	if !ok {
		t.Fatal("expanded entry grch38/star/2.7.11b/gencode47-101 not found")
	}

	wantDeps := []string{"grch38/genome/gencode", "grch38/gtf-gencode/47", "star/2.7.11b"}
	if !reflect.DeepEqual(entry.Deps, wantDeps) {
		t.Errorf("Deps = %v, want %v", entry.Deps, wantDeps)
	}

	wantWhatis := "STAR index GENCODE47 read length 101"
	if entry.Whatis != wantWhatis {
		t.Errorf("Whatis = %q, want %q", entry.Whatis, wantWhatis)
	}
}

func TestExpandTemplates_EmptyDeps(t *testing.T) {
	scripts := map[string]ScriptInfo{
		"grch38/gtf-gencode": {
			IsTemplate:     true,
			TargetTemplate: "grch38/gtf-gencode/{gencode_version}",
			PLOrder:        []string{"gencode_version"},
			Deps:           []string{}, // no deps — stored as empty, not nil
			PL:             map[string][]string{"gencode_version": {"47"}},
		},
	}

	expandTemplates(scripts)

	entry, ok := scripts["grch38/gtf-gencode/47"]
	if !ok {
		t.Fatal("expanded entry grch38/gtf-gencode/47 not found")
	}
	// Empty deps list: resolveBuildSource will set base.dependencies = []string{},
	// so parseDependencies correctly skips file re-reading.
	if len(entry.Deps) != 0 {
		t.Errorf("Deps = %v, want empty", entry.Deps)
	}
}

func TestExpandTemplates_OpenEndedSkipped(t *testing.T) {
	// read_length has "*" — only concrete values (101, 151) should be expanded.
	scripts := map[string]ScriptInfo{
		"grch38/star-gencode": {
			IsTemplate:     true,
			TargetTemplate: "grch38/star/{star_version}/gencode{gencode_version}-{read_length}",
			PLOrder:        []string{"star_version", "gencode_version", "read_length"},
			Deps:           []string{"grch38/gtf-gencode/{gencode_version}", "star/{star_version}"},
			PL: map[string][]string{
				"star_version":    {"2.7.11b"},
				"gencode_version": {"47"},
				"read_length":     {"151", "101", "*"},
			},
		},
	}

	expandTemplates(scripts)

	if _, ok := scripts["grch38/star/2.7.11b/gencode47-151"]; !ok {
		t.Error("expected expanded entry for read_length=151")
	}
	if _, ok := scripts["grch38/star/2.7.11b/gencode47-101"]; !ok {
		t.Error("expected expanded entry for read_length=101")
	}
	// "*" itself should never appear as a concrete value in a target name
	for name := range scripts {
		if strings.Contains(name, "*") {
			t.Errorf("unexpected entry with literal '*' in name: %s", name)
		}
	}
}

// ---------------------------------------------------------------------------
// matchTarget
// ---------------------------------------------------------------------------

func TestMatchTarget(t *testing.T) {
	tmpl := "grch38/star/{star_version}/gencode{gencode_version}-{read_length}"

	cases := []struct {
		name    string
		wantOK  bool
		wantVar map[string]string
	}{
		{
			name:   "grch38/star/2.7.11b/gencode47-75",
			wantOK: true,
			wantVar: map[string]string{
				"star_version":    "2.7.11b",
				"gencode_version": "47",
				"read_length":     "75",
			},
		},
		{
			name:   "grch38/star/2.7.11b/gencode47-101",
			wantOK: true,
			wantVar: map[string]string{
				"star_version":    "2.7.11b",
				"gencode_version": "47",
				"read_length":     "101",
			},
		},
		{
			name:   "grch38/star/2.7.10a/gencode45-150",
			wantOK: true,
			wantVar: map[string]string{
				"star_version":    "2.7.10a",
				"gencode_version": "45",
				"read_length":     "150",
			},
		},
		// wrong shape — should not match
		{name: "grch38/salmon/1.10.2/gencode47", wantOK: false},
		{name: "grch38/star/2.7.11b/gencode47", wantOK: false},
		{name: "", wantOK: false},
	}

	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			got, ok := matchTarget(tc.name, tmpl, nil)
			if ok != tc.wantOK {
				t.Fatalf("matchTarget(%q) ok=%v, want %v", tc.name, ok, tc.wantOK)
			}
			if tc.wantOK && !reflect.DeepEqual(got, tc.wantVar) {
				t.Errorf("vars = %v, want %v", got, tc.wantVar)
			}
		})
	}
}

func TestMatchTarget_Salmon(t *testing.T) {
	tmpl := "grch38/salmon/{salmon_version}/gencode{gencode_version}"
	vars, ok := matchTarget("grch38/salmon/1.10.2/gencode47", tmpl, nil)
	if !ok {
		t.Fatal("expected match")
	}
	want := map[string]string{"salmon_version": "1.10.2", "gencode_version": "47"}
	if !reflect.DeepEqual(vars, want) {
		t.Errorf("vars = %v, want %v", vars, want)
	}
}

func TestMatchTarget_CellRanger(t *testing.T) {
	tmpl := "grch38/star/{star_version}/cellranger-{cr_annotation_version}"
	vars, ok := matchTarget("grch38/star/2.7.11b/cellranger-2024-A", tmpl, nil)
	if !ok {
		t.Fatal("expected match")
	}
	want := map[string]string{"star_version": "2.7.11b", "cr_annotation_version": "2024-A"}
	if !reflect.DeepEqual(vars, want) {
		t.Errorf("vars = %v, want %v", vars, want)
	}
}

func TestMatchTarget_NoPlaceholders(t *testing.T) {
	_, ok := matchTarget("grch38/genome/gencode", "grch38/genome/gencode", nil)
	if ok {
		t.Error("template with no placeholders should not match")
	}
}

// ---------------------------------------------------------------------------
// matchTemplateTarget
// ---------------------------------------------------------------------------

func TestMatchTemplateTarget(t *testing.T) {
	scripts := map[string]ScriptInfo{
		"grch38/star-gencode": {
			IsTemplate:     true,
			TargetTemplate: "grch38/star/{star_version}/gencode{gencode_version}-{read_length}",
			PLOrder:        []string{"star_version", "gencode_version", "read_length"},
			Path:           "grch38/star-gencode",
			Whatis:         "STAR GRCh38 GENCODE{gencode_version} index for read length {read_length}",
			Deps:           []string{"grch38/gtf-gencode/{gencode_version}", "star/{star_version}"},
		},
		"grch38/salmon-gencode": {
			IsTemplate:     true,
			TargetTemplate: "grch38/salmon/{salmon_version}/gencode{gencode_version}",
			PLOrder:        []string{"salmon_version", "gencode_version"},
			Path:           "grch38/salmon-gencode",
			Whatis:         "Salmon GRCh38 GENCODE{gencode_version} index",
			Deps:           []string{"grch38/transcript-gencode/{gencode_version}", "salmon/{salmon_version}"},
		},
	}

	t.Run("star custom read_length", func(t *testing.T) {
		info, found := matchTemplateTarget("grch38/star/2.7.11b/gencode47-75", scripts)
		if !found {
			t.Fatal("expected match")
		}
		if info.Name != "grch38/star/2.7.11b/gencode47-75" {
			t.Errorf("Name = %q", info.Name)
		}
		if info.PL["read_length"][0] != "75" {
			t.Errorf("read_length = %q, want 75", info.PL["read_length"])
		}
		if info.PL["star_version"][0] != "2.7.11b" {
			t.Errorf("star_version = %q, want 2.7.11b", info.PL["star_version"])
		}
		wantDeps := []string{"grch38/gtf-gencode/47", "star/2.7.11b"}
		if !reflect.DeepEqual(info.Deps, wantDeps) {
			t.Errorf("Deps = %v, want %v", info.Deps, wantDeps)
		}
		wantWhatis := "STAR GRCh38 GENCODE47 index for read length 75"
		if info.Whatis != wantWhatis {
			t.Errorf("Whatis = %q, want %q", info.Whatis, wantWhatis)
		}
	})

	t.Run("salmon", func(t *testing.T) {
		info, found := matchTemplateTarget("grch38/salmon/1.10.2/gencode45", scripts)
		if !found {
			t.Fatal("expected match")
		}
		if info.PL["salmon_version"][0] != "1.10.2" {
			t.Errorf("salmon_version = %q", info.PL["salmon_version"])
		}
		if info.PL["gencode_version"][0] != "45" {
			t.Errorf("gencode_version = %q", info.PL["gencode_version"])
		}
	})

	t.Run("no match", func(t *testing.T) {
		_, found := matchTemplateTarget("grch38/bowtie2/something", scripts)
		if found {
			t.Error("expected no match")
		}
	})

	t.Run("empty name", func(t *testing.T) {
		_, found := matchTemplateTarget("", scripts)
		if found {
			t.Error("expected no match for empty name")
		}
	})

	t.Run("non-template entries ignored", func(t *testing.T) {
		plain := map[string]ScriptInfo{
			"grch38/genome/gencode": {IsTemplate: false, Path: "grch38/genome/gencode"},
		}
		_, found := matchTemplateTarget("grch38/genome/gencode", plain)
		if found {
			t.Error("non-template entries should not match")
		}
	})

	t.Run("salmon deps interpolated", func(t *testing.T) {
		info, found := matchTemplateTarget("grch38/salmon/1.10.2/gencode45", scripts)
		if !found {
			t.Fatal("expected match")
		}
		wantDeps := []string{"grch38/transcript-gencode/45", "salmon/1.10.2"}
		if !reflect.DeepEqual(info.Deps, wantDeps) {
			t.Errorf("Deps = %v, want %v", info.Deps, wantDeps)
		}
	})

	t.Run("template with no deps returns nil deps", func(t *testing.T) {
		nodeps := map[string]ScriptInfo{
			"grch38/genome-gencode": {
				IsTemplate:     true,
				TargetTemplate: "grch38/genome/{gencode_version}",
				PLOrder:        []string{"gencode_version"},
				Path:           "grch38/genome-gencode",
				Deps:           nil,
			},
		}
		info, found := matchTemplateTarget("grch38/genome/47", nodeps)
		if !found {
			t.Fatal("expected match")
		}
		if info.Deps != nil {
			t.Errorf("Deps should be nil when template has no deps, got %v", info.Deps)
		}
	})
}
