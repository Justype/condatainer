package catalog

import (
	"errors"
	"os"
	"path/filepath"
	"slices"
	"strings"
	"testing"
)

func openFake(t *testing.T) Catalog {
	t.Helper()
	cat, err := Open(t.Context(), []Spec{{Name: "local", Base: fakeSource(t)}}, Cache{})
	if err != nil {
		t.Fatal(err)
	}
	return cat
}

func TestLookup(t *testing.T) {
	cat := openFake(t)

	// An index key resolves to itself.
	m, found, err := cat.Lookup(t.Context(), "cellranger/9.0.1")
	if err != nil || !found {
		t.Fatalf("found=%v err=%v", found, err)
	}
	if m.Entry.Kind != KindApp || m.Source.Name != "local" || len(m.Vars) != 0 {
		t.Errorf("match = %+v", m)
	}

	// Spellings of one name are the same name.
	for _, spelling := range []string{"cellranger=9.0.1", "cellranger@9.0.1", " cellranger/9.0.1 "} {
		if _, ok, _ := cat.Lookup(t.Context(), spelling); !ok {
			t.Errorf("Lookup(%q) not found", spelling)
		}
	}

	// A bare template name hands back the template and fills nothing in.
	m, found, err = cat.Lookup(t.Context(), "grch38/star-gencode")
	if err != nil || !found {
		t.Fatalf("template: found=%v err=%v", found, err)
	}
	if !m.Entry.IsTemplate || len(m.Vars) != 0 {
		t.Errorf("bare template returned vars %v, want none", m.Vars)
	}

	// A full target recovers every value.
	m, found, err = cat.Lookup(t.Context(), "grch38/star/2.7.11b/gencode47-101")
	if err != nil || !found {
		t.Fatalf("target: found=%v err=%v", found, err)
	}
	want := map[string]string{"star_version": "2.7.11b", "gencode_version": "47", "read_length": "101"}
	for k, v := range want {
		if m.Vars[k] != v {
			t.Errorf("Vars[%q] = %q, want %q", k, m.Vars[k], v)
		}
	}
}

func TestLookupNotProvided(t *testing.T) {
	cat := openFake(t)
	for _, name := range []string{
		"samtools/1.23.1",              // a conda package: no recipe anywhere
		"grch38/star/2.7.11b",          // partially filled, determines nothing
		"grch38/star/9.9.9/gencode1-1", // outside the #PH: sets
		"",
	} {
		m, found, err := cat.Lookup(t.Context(), name)
		if err != nil {
			t.Errorf("Lookup(%q): %v", name, err)
		}
		if found || m != nil {
			t.Errorf("Lookup(%q) = %+v, want not found", name, m)
		}
	}
}

func TestVersions(t *testing.T) {
	root := fakeSource(t)
	// A second version of one app, and a template with a version axis.
	write := func(rel, body string) {
		p := filepath.Join(root, filepath.FromSlash(rel))
		if err := os.MkdirAll(filepath.Dir(p), 0o755); err != nil {
			t.Fatal(err)
		}
		if err := os.WriteFile(p, []byte(body), 0o644); err != nil {
			t.Fatal(err)
		}
	}
	write("recipes/cellranger/10.0.0", "#DESCRIPTION:newer\n")
	write("recipes/ubuntu24/r.def", "#DESCRIPTION:R {version}\n#TARGET:ubuntu24/r/{version}\n#PH:version:4.5.3,4.6.1,4.4.0\n")

	cat, err := Open(t.Context(), []Spec{{Name: "local", Base: root}}, Cache{})
	if err != nil {
		t.Fatal(err)
	}

	// A plain recipe has one index key per version.
	got, err := cat.Versions(t.Context(), "cellranger")
	if err != nil {
		t.Fatal(err)
	}
	if want := []string{"10.0.0", "9.0.1"}; !slices.Equal(got, want) {
		t.Errorf("cellranger versions = %v, want %v", got, want)
	}

	// A template takes them from #PH:, newest first.
	got, err = cat.Versions(t.Context(), "ubuntu24/r")
	if err != nil {
		t.Fatal(err)
	}
	if want := []string{"4.6.1", "4.5.3", "4.4.0"}; !slices.Equal(got, want) {
		t.Errorf("r versions = %v, want %v", got, want)
	}

	// A multi-axis template has no version axis, so no answer.
	if got, _ := cat.Versions(t.Context(), "grch38/star-gencode"); len(got) != 0 {
		t.Errorf("star-gencode versions = %v, want none", got)
	}
	if got, _ := cat.Versions(t.Context(), "nothing-here"); len(got) != 0 {
		t.Errorf("unknown name versions = %v, want none", got)
	}
}

func TestCatalogOpen(t *testing.T) {
	cat := openFake(t)

	// A non-template comes back as written.
	rec, err := cat.Open(t.Context(), "cellranger/9.0.1", nil)
	if err != nil {
		t.Fatal(err)
	}
	if rec.Name != "cellranger/9.0.1" || rec.Description != "cellranger" {
		t.Errorf("recipe = %+v", rec.Entry)
	}

	// A full target name supplies its own vars.
	rec, err = cat.Open(t.Context(), "grch38/star/2.7.11b/gencode47-101", nil)
	if err != nil {
		t.Fatal(err)
	}
	if rec.Name != "grch38/star/2.7.11b/gencode47-101" || rec.IsTemplate {
		t.Errorf("expanded = %+v", rec.Entry)
	}
	if want := "STAR 2.7.11b index for GENCODE 47"; rec.Description != want {
		t.Errorf("Description = %q, want %q", rec.Description, want)
	}

	// A bare template name needs the caller's choice, and says so without it.
	if _, err := cat.Open(t.Context(), "grch38/star-gencode", nil); err == nil {
		t.Error("opening a bare template without vars should fail")
	}
	rec, err = cat.Open(t.Context(), "grch38/star-gencode", map[string]string{
		"star_version": "2.7.9a", "gencode_version": "49", "read_length": "151"})
	if err != nil {
		t.Fatal(err)
	}
	if want := "grch38/star/2.7.9a/gencode49-151"; rec.Name != want {
		t.Errorf("Name = %q, want %q", rec.Name, want)
	}

	// A name nothing provides is ErrNotProvided, for the caller's fallback.
	_, err = cat.Open(t.Context(), "samtools/1.23.1", nil)
	if !errors.Is(err, ErrNotProvided) {
		t.Errorf("err = %v, want ErrNotProvided", err)
	}
}

func TestReadPath(t *testing.T) {
	cat := openFake(t)
	data, err := cat.ReadPath(t.Context(), cat[0], "recipes/cellranger/9.0.1")
	if err != nil {
		t.Fatal(err)
	}
	if !strings.Contains(string(data), "#DESCRIPTION:cellranger") {
		t.Errorf("ReadPath returned %q", data)
	}
}

// Earlier sources shadow later ones, per entry rather than wholesale.
func TestLookupFirstSourceWins(t *testing.T) {
	lab := t.TempDir()
	p := filepath.Join(lab, "recipes", "cellranger", "9.0.1")
	if err := os.MkdirAll(filepath.Dir(p), 0o755); err != nil {
		t.Fatal(err)
	}
	if err := os.WriteFile(p, []byte("#DESCRIPTION:lab build\n"), 0o644); err != nil {
		t.Fatal(err)
	}

	cat, err := Open(t.Context(), []Spec{
		{Name: "lab", Base: lab},
		{Name: "cnt", Base: fakeSource(t)},
	}, Cache{})
	if err != nil {
		t.Fatal(err)
	}

	m, found, err := cat.Lookup(t.Context(), "cellranger/9.0.1")
	if err != nil || !found {
		t.Fatalf("found=%v err=%v", found, err)
	}
	if m.Source.Name != "lab" || m.Entry.Description != "lab build" {
		t.Errorf("match came from %q (%q), want lab", m.Source.Name, m.Entry.Description)
	}
	// Shadowing is per item: what only the later source has still resolves.
	if _, ok, _ := cat.Lookup(t.Context(), "ubuntu24/base"); !ok {
		t.Error("the later source should still provide what the earlier lacks")
	}
}
