package catalog

import (
	"context"
	"os"
	"path/filepath"
	"testing"
)

// solveSource writes a distro-qualified template (two versions) and a bare
// app (two versions) — enough to exercise every SolveName attempt.
func solveSource(t *testing.T) Catalog {
	t.Helper()
	root := t.TempDir()
	write := func(rel, body string) {
		p := filepath.Join(root, filepath.FromSlash(rel))
		if err := os.MkdirAll(filepath.Dir(p), 0o755); err != nil {
			t.Fatal(err)
		}
		if err := os.WriteFile(p, []byte(body), 0o644); err != nil {
			t.Fatal(err)
		}
	}
	write("recipes/ubuntu24/rstudio-server.def",
		"#DESC:RStudio Server {version}\n#TARGET:ubuntu24/rstudio-server/{version}\n#PH:version:2025.05.0-160,2026.09.0-174\n")
	write("recipes/ubuntu24/xfce4.def", "#DESC:XFCE4\n")
	write("recipes/openjdk/17.0.15", "#DESC:OpenJDK 17.0.15\n")
	write("recipes/openjdk/17.0.18", "#DESC:OpenJDK 17.0.18\n")
	write("recipes/cellranger/9.0.1", "#DESC:cellranger\n")

	cat, err := Open(t.Context(), []Spec{{Name: "local", Base: root}}, Cache{})
	if err != nil {
		t.Fatal(err)
	}
	return cat
}

func TestSolveNameExact(t *testing.T) {
	cat := solveSource(t)

	res, found, err := cat.SolveName(t.Context(), nil, "", "cellranger/9.0.1")
	if err != nil || !found {
		t.Fatalf("found=%v err=%v", found, err)
	}
	if res.Name != "cellranger/9.0.1" || res.Entry == nil {
		t.Errorf("res = %+v, want the exact entry", res)
	}
}

func TestSolveNameBarePicksNewestAvailable(t *testing.T) {
	cat := solveSource(t)

	res, found, err := cat.SolveName(t.Context(), nil, "", "openjdk")
	if err != nil || !found {
		t.Fatalf("found=%v err=%v", found, err)
	}
	if res.Name != "openjdk/17.0.18" || res.Installed != "" {
		t.Errorf("res = %+v, want newest openjdk/17.0.18, nothing installed", res)
	}
}

func TestSolveNamePartialVersionFamilyMatch(t *testing.T) {
	cat := solveSource(t)

	res, found, err := cat.SolveName(t.Context(), nil, "", "openjdk/17")
	if err != nil || !found {
		t.Fatalf("found=%v err=%v", found, err)
	}
	if res.Name != "openjdk/17.0.18" {
		t.Errorf("res.Name = %q, want openjdk/17.0.18", res.Name)
	}
}

func TestSolveNameInstalledBeatsNewer(t *testing.T) {
	cat := solveSource(t)
	have := func(name string) []string {
		if name == "openjdk" {
			return []string{"17.0.15"}
		}
		return nil
	}

	res, found, err := cat.SolveName(t.Context(), have, "", "openjdk/17")
	if err != nil || !found {
		t.Fatalf("found=%v err=%v", found, err)
	}
	if res.Name != "openjdk/17.0.15" || res.Installed != "17.0.15" {
		t.Errorf("res = %+v, want the installed 17.0.15 reused", res)
	}
}

// have answers before the catalog is ever consulted — an installed overlay
// resolves even when the catalog has nothing backing it at all (empty
// source, or the recipe has since moved or dropped), matching the "check
// what's already installed first, no network/catalog call" rule.
func TestSolveNameInstalledResolvesWithEmptyCatalog(t *testing.T) {
	root := t.TempDir()
	cat, err := Open(t.Context(), []Spec{{Name: "empty", Base: root}}, Cache{})
	if err != nil {
		t.Fatal(err)
	}
	have := func(name string) []string {
		if name == "ubuntu24/rstudio-server" {
			return []string{"2026.09.0-174"}
		}
		if name == "openjdk" {
			return []string{"17.0.15"}
		}
		return nil
	}

	res, found, err := cat.SolveName(t.Context(), have, "ubuntu24", "rstudio-server")
	if err != nil || !found {
		t.Fatalf("found=%v err=%v", found, err)
	}
	if res.Name != "ubuntu24/rstudio-server/2026.09.0-174" || res.Installed != "2026.09.0-174" {
		t.Errorf("res = %+v, want the installed version with no catalog backing", res)
	}
	if res.Entry != nil {
		t.Errorf("res.Entry = %+v, want nil — nothing in the catalog backs this", res.Entry)
	}

	res, found, err = cat.SolveName(t.Context(), have, "", "openjdk")
	if err != nil || !found {
		t.Fatalf("found=%v err=%v", found, err)
	}
	if res.Name != "openjdk/17.0.15" || res.Installed != "17.0.15" {
		t.Errorf("res = %+v, want the installed app version with no catalog backing", res)
	}
}

func TestSolveNameDistroPrefixRetry(t *testing.T) {
	cat := solveSource(t)

	// Bare, no distro prefix, no version: only found once retried under distro.
	if _, found, _ := cat.SolveName(t.Context(), nil, "", "rstudio-server"); found {
		t.Fatal("no distro supplied: should not resolve")
	}

	res, found, err := cat.SolveName(t.Context(), nil, "ubuntu24", "rstudio-server")
	if err != nil || !found {
		t.Fatalf("found=%v err=%v", found, err)
	}
	if res.Name != "ubuntu24/rstudio-server/2026.09.0-174" {
		t.Errorf("res.Name = %q, want the newest under the distro prefix", res.Name)
	}
}

func TestSolveNameNoDistroRetryPastOneSlash(t *testing.T) {
	cat := solveSource(t)

	// Already has an internal slash: a distro prefix would add a second, so
	// it must not be tried even when nothing resolves as given.
	_, found, err := cat.SolveName(t.Context(), nil, "ubuntu24", "grch38/genome")
	if err != nil || found {
		t.Errorf("found=%v err=%v, want not found without a distro-prefix retry", found, err)
	}
}

func TestSolveNameNotFound(t *testing.T) {
	cat := solveSource(t)

	for _, raw := range []string{"nothing-here", "openjdk/99", "", "bioconda::star/2.7.11b"} {
		res, found, err := cat.SolveName(t.Context(), nil, "ubuntu24", raw)
		if err != nil || found {
			t.Errorf("SolveName(%q) = %+v, found=%v err=%v, want not found", raw, res, found, err)
		}
	}
}

// A flat, versionless entry directly under a distro — "ubuntu24/xfce4", no
// #PH: axis — must still resolve even though "ubuntu24" is itself a
// registered distro (via this very entry's own os type): the bare-distro
// guard is for a request with no version to satisfy, not for every dep whose
// Name happens to be a distro.
func TestSolveNameFlatKeyUnderRegisteredDistro(t *testing.T) {
	cat := solveSource(t)

	res, found, err := cat.SolveName(t.Context(), nil, "", "ubuntu24/xfce4")
	if err != nil || !found {
		t.Fatalf("found=%v err=%v", found, err)
	}
	if res.Name != "ubuntu24/xfce4" {
		t.Errorf("res.Name = %q, want ubuntu24/xfce4", res.Name)
	}

	// And by the bare-name distro-prefix retry, same as any other bare name.
	res, found, err = cat.SolveName(t.Context(), nil, "ubuntu24", "xfce4")
	if err != nil || !found {
		t.Fatalf("found=%v err=%v", found, err)
	}
	if res.Name != "ubuntu24/xfce4" {
		t.Errorf("res.Name = %q, want ubuntu24/xfce4", res.Name)
	}
}

// dataSource writes two flat, non-templated data entries sharing a prefix —
// alternative sources for the same reference, never versions of one
// another — enough to exercise every "no autofill for data" case.
func dataSource(t *testing.T) Catalog {
	t.Helper()
	root := t.TempDir()
	write := func(rel, body string) {
		p := filepath.Join(root, filepath.FromSlash(rel))
		if err := os.MkdirAll(filepath.Dir(p), 0o755); err != nil {
			t.Fatal(err)
		}
		if err := os.WriteFile(p, []byte(body), 0o644); err != nil {
			t.Fatal(err)
		}
	}
	write("recipes/grch38/genome/ucsc", "#DESC:UCSC genome\n")
	write("recipes/grch38/genome/ensembl", "#DESC:Ensembl genome\n")

	cat, err := Open(t.Context(), []Spec{{Name: "data", Base: root}}, Cache{})
	if err != nil {
		t.Fatal(err)
	}
	return cat
}

// A bare data name must not resolve to whichever alternative source sorts
// highest — "ucsc" and "ensembl" are different data, not different versions
// of the same data, so CompareVersions between them means nothing.
func TestSolveNameBareDataNotFound(t *testing.T) {
	cat := dataSource(t)

	res, found, err := cat.SolveName(t.Context(), nil, "", "grch38/genome")
	if err != nil || found {
		t.Errorf("SolveName(%q) = %+v, found=%v err=%v, want not found", "grch38/genome", res, found, err)
	}
}

// A data name given in full still resolves normally: only a genuine choice
// among data siblings is refused, not the case where exactly one satisfies.
func TestSolveNameExactDataStillResolves(t *testing.T) {
	cat := dataSource(t)

	res, found, err := cat.SolveName(t.Context(), nil, "", "grch38/genome/ucsc")
	if err != nil || !found {
		t.Fatalf("found=%v err=%v", found, err)
	}
	if res.Name != "grch38/genome/ucsc" {
		t.Errorf("res.Name = %q, want grch38/genome/ucsc", res.Name)
	}
}

// A bare distro name is a namespace, not a module: it must not resolve to
// whichever sibling app (xfce4, rstudio-server, ...) happens to sort highest.
func TestSolveNameBareDistroNotFound(t *testing.T) {
	cat := solveSource(t)

	res, found, err := cat.SolveName(t.Context(), nil, "", "ubuntu24")
	if err != nil || found {
		t.Errorf("SolveName(%q) = %+v, found=%v err=%v, want not found", "ubuntu24", res, found, err)
	}
	// Even when it's also the distro being retried under — no self-prefixing.
	res, found, err = cat.SolveName(t.Context(), nil, "ubuntu24", "ubuntu24")
	if err != nil || found {
		t.Errorf("SolveName(%q, distro %q) = %+v, found=%v err=%v, want not found", "ubuntu24", "ubuntu24", res, found, err)
	}
}

func TestSolveInstalledFlatOverlayWithoutRecipe(t *testing.T) {
	installed := map[string][]string{"ubuntu24/xfce4": {"/images/ubuntu24--xfce4.sqf"}}

	name, path, found, err := SolveInstalled(context.Background(), installed, "ubuntu24", "xfce4")
	if err != nil || !found || name != "ubuntu24/xfce4" || path != "/images/ubuntu24--xfce4.sqf" {
		t.Fatalf("got %q %q %v %v", name, path, found, err)
	}
	if _, _, found, _ := SolveInstalled(context.Background(), installed, "ubuntu24", "nothing"); found {
		t.Error("an uninstalled name resolved")
	}
}
