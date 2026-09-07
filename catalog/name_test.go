package catalog

import "testing"

func TestNormalize(t *testing.T) {
	tests := []struct{ in, want string }{
		{"samtools/1.2", "samtools/1.2"},
		{"samtools=1.2", "samtools/1.2"},
		{"samtools@1.2", "samtools/1.2"},
		{"samtools--1.2", "samtools/1.2"},
		{"  samtools@1.2  ", "samtools/1.2"},
		{"samtools", "samtools"},
		{"samtools@1.21>=1.16", "samtools/1.21>=1.16"},
		{"grch38/star/2.7.11b/gencode47-101", "grch38/star/2.7.11b/gencode47-101"},
		{"", ""},
	}
	for _, tt := range tests {
		if got := Normalize(tt.in); got != tt.want {
			t.Errorf("Normalize(%q) = %q, want %q", tt.in, got, tt.want)
		}
	}
}

func TestParseDep(t *testing.T) {
	tests := []struct {
		raw  string
		want Dep
	}{
		{"samtools/1.23.1>=1.10", Dep{Name: "samtools", Version: "1.23.1", Op: ">=", Min: "1.10"}},
		{"star/2.7.11b", Dep{Name: "star", Version: "2.7.11b"}},
		{"star/2.7.11b>2.7.0", Dep{Name: "star", Version: "2.7.11b", Op: ">", Min: "2.7.0"}},
		{"samtools", Dep{Name: "samtools"}},
		// A data dep: the last component is the version, so the name keeps two.
		{"grch38/genome/gencode", Dep{Name: "grch38/genome", Version: "gencode"}},
		// A template dep is unsubstituted until the parent expands.
		{"grch38/gtf-gencode/{gencode_version}", Dep{Name: "grch38/gtf-gencode", Version: "{gencode_version}"}},
	}
	for _, tt := range tests {
		got, err := ParseDep(tt.raw)
		if err != nil {
			t.Errorf("ParseDep(%q): %v", tt.raw, err)
			continue
		}
		if got != tt.want {
			t.Errorf("ParseDep(%q) = %+v, want %+v", tt.raw, got, tt.want)
		}
		if s := got.String(); s != Normalize(tt.raw) {
			t.Errorf("ParseDep(%q).String() = %q, want %q", tt.raw, s, Normalize(tt.raw))
		}
	}
	if _, err := ParseDep("  "); err != ErrEmptyDep {
		t.Errorf("ParseDep(blank) error = %v, want ErrEmptyDep", err)
	}
}

func TestDepSatisfies(t *testing.T) {
	// The preferred version is the implicit upper bound: [1.10, 1.23.1].
	ranged, _ := ParseDep("samtools/1.23.1>=1.10")
	for _, tt := range []struct {
		version string
		want    bool
	}{
		{"1.10", true},
		{"1.21", true},
		{"1.23.1", true},
		{"1.9", false},
		{"1.24", false},
	} {
		if got := ranged.Satisfies(tt.version); got != tt.want {
			t.Errorf("%s.Satisfies(%q) = %v, want %v", ranged, tt.version, got, tt.want)
		}
	}

	// Unconstrained: only the preferred version will do.
	pinned, _ := ParseDep("star/2.7.11b")
	if !pinned.Satisfies("2.7.11b") || pinned.Satisfies("2.7.11a") {
		t.Error("unconstrained dep should admit only its preferred version")
	}
}

func TestDeriveType(t *testing.T) {
	tests := []struct {
		name, target string
		isDef        bool
		declared     string
		want         Type
	}{
		{name: "ubuntu24/build-essential", isDef: true, want: TypeOS},
		{name: "ubuntu24/r", target: "ubuntu24/r/{version}", isDef: true, want: TypeOS},
		{name: "cellranger/9.0.1", want: TypeApp},
		{name: "cytoscape", want: TypeApp},
		{name: "grch38/genome/gencode", want: TypeData},
		// A template's kind comes from its target: one slash as a filename,
		// three as the module path it builds.
		{name: "grch38/star-gencode", target: "grch38/star/{star_version}/gencode{gencode_version}-{read_length}", want: TypeData},
		// #TYPE: overrides app/data only.
		{name: "cellranger/9.0.1", declared: "data", want: TypeData},
		{name: "grch38/genome/gencode", declared: "app", want: TypeApp},
		// #TYPE: has no effect on a .def: every one is os regardless.
		{name: "ubuntu24/build-essential", isDef: true, declared: "app", want: TypeOS},
		{name: "cellranger/9.0.1", declared: "nonsense", want: TypeApp},
	}
	for _, tt := range tests {
		if got := DeriveType(tt.name, tt.target, tt.isDef, tt.declared); got != tt.want {
			t.Errorf("DeriveType(%q, %q, %v, %q) = %q, want %q",
				tt.name, tt.target, tt.isDef, tt.declared, got, tt.want)
		}
	}
}
