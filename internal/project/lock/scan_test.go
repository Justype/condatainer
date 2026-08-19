package lock

import (
	"errors"
	"os"
	"path/filepath"
	"strings"
	"testing"
)

func write(t *testing.T, root, rel, body string) string {
	t.Helper()
	path := filepath.Join(root, rel)
	if err := os.MkdirAll(filepath.Dir(path), 0o775); err != nil {
		t.Fatal(err)
	}
	if err := os.WriteFile(path, []byte(body), 0o664); err != nil {
		t.Fatal(err)
	}
	return path
}

func scan(t *testing.T, root string) *ScanResult {
	t.Helper()
	result, err := Scan(root, ScanOptions{})
	if err != nil {
		t.Fatal(err)
	}
	return result
}

func keys(result *ScanResult) []string {
	out := make([]string, 0, len(result.Requests))
	for _, request := range result.Requests {
		out = append(out, request.Key)
	}
	return out
}

func find(t *testing.T, result *ScanResult, key string) Request {
	t.Helper()
	for _, request := range result.Requests {
		if request.Key == key {
			return request
		}
	}
	t.Fatalf("no request %q in %v", key, keys(result))
	return Request{}
}

func TestScanFindsDeclarationsAcrossScripts(t *testing.T) {
	root := t.TempDir()
	write(t, root, "analysis.sh", "#!/bin/bash\n#DEP: star/2.7.11b\n#DEP: cutadapt/5.0\nstar --version\n")
	write(t, root, "scripts/align.bash", "#DEP: star/2.7.11b\n#DEP: grch38/genome/gencode49\nalign\n")

	result := scan(t, root)
	want := []string{"cutadapt/5.0", "grch38/genome/gencode49", "star/2.7.11b"}
	got := keys(result)
	if len(got) != len(want) {
		t.Fatalf("requests = %v, want %v", got, want)
	}
	for i := range want {
		if got[i] != want[i] {
			t.Fatalf("requests = %v, want %v (sorted)", got, want)
		}
	}
	star := find(t, result, "star/2.7.11b")
	if len(star.Scripts) != 2 || star.Scripts[0] != "analysis.sh" || star.Scripts[1] != "scripts/align.bash" {
		t.Errorf("star scripts = %v, want both, sorted", star.Scripts)
	}
	if star.Kind != KindName || star.Dep.Name != "star" || star.Dep.Version != "2.7.11b" {
		t.Errorf("star request = %#v", star)
	}
	data := find(t, result, "grch38/genome/gencode49")
	if data.Dep.Name != "grch38/genome" || data.Dep.Version != "gencode49" {
		t.Errorf("a slash-carrying name was mis-split: %#v", data.Dep)
	}
}

// A range is a build-recipe feature. An analysis names one exact version, so a
// constrained declaration is a finding that says what to write instead.
func TestScanRejectsAVersionConstraint(t *testing.T) {
	root := t.TempDir()
	write(t, root, "run.sh", "#DEP: star/2.7.11b>=2.7.0\nrun\n")

	result := scan(t, root)
	if len(result.Requests) != 0 {
		t.Fatalf("requests = %v, want the constrained declaration refused", keys(result))
	}
	if len(result.Findings) != 1 {
		t.Fatalf("findings = %#v, want one", result.Findings)
	}
	if !strings.Contains(result.Findings[0].Reason, "star/2.7.11b") {
		t.Errorf("finding does not name the exact version to write: %q", result.Findings[0].Reason)
	}
	if !strings.Contains(result.Findings[0].Reason, "build recipe") {
		t.Errorf("finding does not say where a range belongs: %q", result.Findings[0].Reason)
	}
}

// Position carries no meaning: a declaration below the first command counts,
// and so does one inside a heredoc that writes another script. The scanner and
// the runtime therefore read a script the same way.
func TestScanReadsDeclarationsAnywhere(t *testing.T) {
	root := t.TempDir()
	write(t, root, "run.sh", "#!/bin/bash\n#DEP: star/2.7.11b\necho hello\n#DEP: cutadapt/5.0\n")

	result := scan(t, root)
	want := []string{"cutadapt/5.0", "star/2.7.11b"}
	if got := keys(result); len(got) != 2 || got[0] != want[0] || got[1] != want[1] {
		t.Fatalf("requests = %v, want %v", got, want)
	}
	if len(result.Findings) != 0 {
		t.Errorf("findings = %#v, want none", result.Findings)
	}
}

func TestScanReadsHeredocDeclarations(t *testing.T) {
	root := t.TempDir()
	write(t, root, "gen.sh", "#!/bin/bash\n#DEP: star/2.7.11b\ncat <<'END' > out.sh\n#DEP: cutadapt/5.0\nEND\n")

	if got := keys(scan(t, root)); len(got) != 2 {
		t.Fatalf("requests = %v, want the heredoc declaration included", got)
	}
}

func TestScanClassifiesPathDeclarations(t *testing.T) {
	root := t.TempDir()
	write(t, root, "run.sh", "#DEP: ./overlays/tool.sqf\n#DEP: env.img  ## unpinned — scratch, rebuilt per machine\nrun\n")

	result := scan(t, root)
	sqf := find(t, result, PathPrefix+"overlays/tool.sqf")
	if sqf.Kind != KindPath || sqf.Path != "overlays/tool.sqf" {
		t.Errorf("sqf request = %#v", sqf)
	}
	if sqf.Unpinned {
		t.Errorf("a .sqf is lockable and must not be marked unpinned by default")
	}
	img := find(t, result, PathPrefix+"env.img")
	if img.Kind != KindWritable {
		t.Errorf("img kind = %q, want writable", img.Kind)
	}
	if !img.Unpinned || img.Reason != "scratch, rebuilt per machine" {
		t.Errorf("unpinned marker = %v, reason = %q", img.Unpinned, img.Reason)
	}
}

// The marker is required; the reason is not.
func TestScanAcceptsABareUnpinnedMarker(t *testing.T) {
	root := t.TempDir()
	write(t, root, "run.sh", "#DEP: env.img  ## unpinned\nrun\n")

	img := find(t, scan(t, root), PathPrefix+"env.img")
	if !img.Unpinned || img.Reason != "" {
		t.Fatalf("bare marker = %#v", img)
	}
}

// A note that is not the marker is just a comment and claims nothing.
func TestScanDoesNotTreatAnyNoteAsTheMarker(t *testing.T) {
	root := t.TempDir()
	write(t, root, "run.sh", "#DEP: env.img  ## the scratch environment\nrun\n")

	if img := find(t, scan(t, root), PathPrefix+"env.img"); img.Unpinned {
		t.Fatalf("an ordinary note was read as the unpinned marker: %#v", img)
	}
}

// The claim is about the artifact, so one script marking it settles it.
func TestScanMergesTheMarkerAcrossScripts(t *testing.T) {
	root := t.TempDir()
	write(t, root, "a.sh", "#DEP: env.img\nrun\n")
	write(t, root, "b.sh", "#DEP: env.img  ## unpinned — shared scratch\nrun\n")

	img := find(t, scan(t, root), PathPrefix+"env.img")
	if !img.Unpinned || img.Reason != "shared scratch" {
		t.Fatalf("marker did not merge: %#v", img)
	}
}

func TestScanSkipsGeneratedAndExcludedDirectories(t *testing.T) {
	root := t.TempDir()
	write(t, root, "keep.sh", "#DEP: star/2.7.11b\nrun\n")
	write(t, root, ".git/hooks/pre-commit.sh", "#DEP: never/1.0\nrun\n")
	write(t, root, "cnt-lock/stale.sh", "#DEP: never/2.0\nrun\n")
	write(t, root, "vendor/dep.sh", "#DEP: never/3.0\nrun\n")

	result, err := Scan(root, ScanOptions{ExcludeDirs: []string{"vendor"}})
	if err != nil {
		t.Fatal(err)
	}
	if got := keys(result); len(got) != 1 || got[0] != "star/2.7.11b" {
		t.Fatalf("requests = %v, want only the kept script", got)
	}
}

func TestScanReadsExtensionlessShebangScriptsOnly(t *testing.T) {
	root := t.TempDir()
	write(t, root, "runner", "#!/usr/bin/env bash\n#DEP: star/2.7.11b\nrun\n")
	write(t, root, "notes", "#DEP: never/1.0\n")
	write(t, root, "data.txt", "#DEP: never/2.0\n")

	if got := keys(scan(t, root)); len(got) != 1 || got[0] != "star/2.7.11b" {
		t.Fatalf("requests = %v", got)
	}
}

// A symlink can point outside the checkout, and a lock describes the checkout.
func TestScanDoesNotFollowSymlinks(t *testing.T) {
	root, outside := t.TempDir(), t.TempDir()
	write(t, outside, "external.sh", "#DEP: never/1.0\nrun\n")
	write(t, root, "keep.sh", "#DEP: star/2.7.11b\nrun\n")
	if err := os.Symlink(filepath.Join(outside, "external.sh"), filepath.Join(root, "linked.sh")); err != nil {
		t.Skipf("symlinks unavailable: %v", err)
	}
	if err := os.Symlink(outside, filepath.Join(root, "linkeddir")); err != nil {
		t.Skipf("symlinks unavailable: %v", err)
	}

	result := scan(t, root)
	if got := keys(result); len(got) != 1 || got[0] != "star/2.7.11b" {
		t.Fatalf("requests = %v, want only the real script", got)
	}
	for _, script := range result.Scripts {
		if script != "keep.sh" {
			t.Errorf("scanned a symlink: %s", script)
		}
	}
}

func TestScanStripsInlineComments(t *testing.T) {
	root := t.TempDir()
	write(t, root, "run.sh", "#DEP: star/2.7.11b # the aligner\nrun\n")

	if got := keys(scan(t, root)); len(got) != 1 || got[0] != "star/2.7.11b" {
		t.Fatalf("requests = %v", got)
	}
}

// A project is where its lock is. Nothing searches a parent, so a subdirectory
// of a project is not a project — that is what keeps a stray cnt-lock/ in $HOME
// from enrolling every script beneath it.
func TestRootAtRequiresTheLockInTheDirectoryItself(t *testing.T) {
	root := t.TempDir()
	if err := os.MkdirAll(filepath.Join(root, DirName), 0o775); err != nil {
		t.Fatal(err)
	}
	got, err := RootAt(root)
	if err != nil {
		t.Fatal(err)
	}
	if resolved, _ := filepath.EvalSymlinks(got); resolved != mustEval(t, root) {
		t.Fatalf("RootAt = %q, want %q", got, root)
	}

	deep := filepath.Join(root, "a", "b")
	if err := os.MkdirAll(deep, 0o775); err != nil {
		t.Fatal(err)
	}
	if _, err := RootAt(deep); !errors.Is(err, ErrNoProject) {
		t.Fatalf("RootAt(subdirectory) = %v, want ErrNoProject", err)
	}
	if _, err := RootAt(t.TempDir()); !errors.Is(err, ErrNoProject) {
		t.Fatalf("RootAt(unrelated) = %v, want ErrNoProject", err)
	}
}

// Naming a project means naming it: an explicit root is never walked past.
func TestRootForDoesNotSearchAboveAnExplicitDirectory(t *testing.T) {
	root := t.TempDir()
	if err := os.MkdirAll(filepath.Join(root, DirName), 0o775); err != nil {
		t.Fatal(err)
	}
	child := filepath.Join(root, "child")
	if err := os.MkdirAll(child, 0o775); err != nil {
		t.Fatal(err)
	}
	got, err := RootFor(child, root)
	if err != nil {
		t.Fatal(err)
	}
	if mustEval(t, got) != mustEval(t, child) {
		t.Fatalf("RootFor = %q, want the explicit %q", got, child)
	}
}

func TestLoadTreatsAMissingLockAsEmpty(t *testing.T) {
	root := t.TempDir()
	if err := os.MkdirAll(filepath.Join(root, DirName), 0o775); err != nil {
		t.Fatal(err)
	}
	l, err := Load(root)
	if err != nil {
		t.Fatal(err)
	}
	if l.SchemaVersion != SchemaVersion || len(l.Selections) != 0 {
		t.Fatalf("Load = %#v", l)
	}
}

func mustEval(t *testing.T, path string) string {
	t.Helper()
	resolved, err := filepath.EvalSymlinks(path)
	if err != nil {
		return path
	}
	return resolved
}

// An external .sqf cannot be a restore destination — restore does not own that
// path — so it is unpinnable like a writable .img rather than lockable.
func TestScanClassifiesExternalPathsAsUnpinnable(t *testing.T) {
	root := t.TempDir()
	write(t, root, "run.sh", "#DEP: /shared/lab/genome.sqf  ## unpinned\n"+
		"#DEP: ../outside/tool.sqf  ## unpinned\n"+
		"#DEP: overlays/inside.sqf\nrun\n")

	result := scan(t, root)
	absolute := find(t, result, PathPrefix+"/shared/lab/genome.sqf")
	if absolute.Kind != KindExternal || absolute.Kind.Pinnable() {
		t.Errorf("absolute path = %#v, want an unpinnable external", absolute)
	}
	climbing := find(t, result, PathPrefix+"../outside/tool.sqf")
	if climbing.Kind != KindExternal {
		t.Errorf("climbing path kind = %q, want external", climbing.Kind)
	}
	inside := find(t, result, PathPrefix+"overlays/inside.sqf")
	if inside.Kind != KindPath || !inside.Kind.Pinnable() {
		t.Errorf("project path = %#v, want a pinnable path", inside)
	}
}

// G13: an unpinnable declaration must say so. Without the marker it is a
// finding, which is what makes `project validate` fail rather than pass while
// the project mounts something unpinned.
func TestScanFlagsAnUndeclaredUnpinnableDeclaration(t *testing.T) {
	root := t.TempDir()
	write(t, root, "run.sh", "#DEP: env.img\n#DEP: /shared/genome.sqf\nrun\n")

	result := scan(t, root)
	if len(result.Findings) != 2 {
		t.Fatalf("findings = %#v, want one per undeclared unpinnable", result.Findings)
	}
	for _, finding := range result.Findings {
		if finding.Line == 0 {
			t.Errorf("finding has no line: %#v", finding)
		}
		if !strings.Contains(finding.Reason, UnpinnedMarker) {
			t.Errorf("finding does not name the marker: %q", finding.Reason)
		}
	}
}

// The marker merges across scripts, so the finding is only decided once every
// script has been read — one script declaring it settles it for the project.
func TestScanDoesNotFlagAnUnpinnableMarkedInAnotherScript(t *testing.T) {
	root := t.TempDir()
	write(t, root, "a.sh", "#DEP: env.img\nrun\n")
	write(t, root, "b.sh", "#DEP: env.img  ## unpinned — scratch\nrun\n")

	if findings := scan(t, root).Findings; len(findings) != 0 {
		t.Fatalf("findings = %#v, want none: b.sh declares it", findings)
	}
}
