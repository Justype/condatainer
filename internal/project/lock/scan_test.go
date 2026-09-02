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
	write(t, root, "run.sh", "#DEP: ./overlays/tool.sqf\n#DEP: env.img\nrun\n")

	result := scan(t, root)
	sqf := find(t, result, PathPrefix+"overlays/tool.sqf")
	if sqf.Kind != KindPath || sqf.Path != "overlays/tool.sqf" {
		t.Errorf("sqf request = %#v", sqf)
	}
	if !sqf.Kind.Pinnable() {
		t.Errorf("a project .sqf must be pinnable")
	}
	img := find(t, result, PathPrefix+"env.img")
	if img.Kind != KindWritable || img.Kind.Pinnable() {
		t.Errorf("img request = %#v, want an unpinnable writable", img)
	}
}

// A `##` note is a comment and never changes what a declaration means. Nothing
// in a script can silence an unpinnable dependency: freezing it is the only way
// to close one.
func TestScanTreatsANoteAsAComment(t *testing.T) {
	root := t.TempDir()
	write(t, root, "run.sh", "#DEP: env.img  ## the scratch environment\nrun\n")

	result := scan(t, root)
	if img := find(t, result, PathPrefix+"env.img"); img.Kind != KindWritable {
		t.Fatalf("img request = %#v", img)
	}
	if len(result.Findings) != 1 {
		t.Fatalf("findings = %#v, want the writable declaration reported", result.Findings)
	}
	if !strings.Contains(result.Findings[0].Reason, "overlay freeze") {
		t.Errorf("the finding does not point at freeze: %q", result.Findings[0].Reason)
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

// Extension only. `run` executes every project script with /bin/bash, so no
// other shell's script could run here, and a file with no extension is not read
// at all — sniffing its shebang cost more than it bought.
func TestScanReadsShellExtensionsOnly(t *testing.T) {
	root := t.TempDir()
	write(t, root, "a.sh", "#!/usr/bin/env bash\n#DEP: star/2.7.11b\nrun\n")
	write(t, root, "b.bash", "#DEP: samtools/1.23.1\n")
	write(t, root, "runner", "#!/usr/bin/env bash\n#DEP: never/1.0\nrun\n")
	write(t, root, "analyze", "#!/home/josh/venv/bin/python\n#DEP: never/2.0\n")
	write(t, root, "notes", "#DEP: never/3.0\n")
	write(t, root, "data.txt", "#DEP: never/4.0\n")
	write(t, root, "job.zsh", "#DEP: never/5.0\n")

	got := keys(scan(t, root))
	want := []string{"samtools/1.23.1", "star/2.7.11b"}
	if len(got) != len(want) {
		t.Fatalf("requests = %v, want %v", got, want)
	}
	for i, key := range want {
		if got[i] != key {
			t.Fatalf("requests = %v, want %v", got, want)
		}
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
	if l.SchemaVersion != SchemaVersion || len(l.Pins) != 0 {
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
	write(t, root, "run.sh", "#DEP: /shared/lab/genome.sqf\n"+
		"#DEP: ../outside/tool.sqf\n"+
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

// An unpinnable declaration is always a finding, which is what makes `project
// validate` fail rather than pass while the project mounts something the lock
// cannot reproduce. Each names what closes it.
func TestScanFlagsEveryUnpinnableDeclaration(t *testing.T) {
	root := t.TempDir()
	write(t, root, "run.sh", "#DEP: env.img\n#DEP: /shared/genome.sqf\nrun\n")

	result := scan(t, root)
	if len(result.Findings) != 2 {
		t.Fatalf("findings = %#v, want one per unpinnable declaration", result.Findings)
	}
	for _, finding := range result.Findings {
		if finding.Line == 0 {
			t.Errorf("finding has no line: %#v", finding)
		}
	}
	if !strings.Contains(result.Findings[0].Reason, "overlay freeze") {
		t.Errorf("the writable finding does not point at freeze: %q", result.Findings[0].Reason)
	}
	if !strings.Contains(result.Findings[1].Reason, "copy it under the project") {
		t.Errorf("the external finding does not say what to do: %q", result.Findings[1].Reason)
	}
}

// Repeating a declaration reports it once: the finding is about the dependency,
// not about each line that names it.
func TestScanReportsOneFindingPerUnpinnableDependency(t *testing.T) {
	root := t.TempDir()
	write(t, root, "a.sh", "#DEP: env.img\nrun\n")
	write(t, root, "b.sh", "#DEP: env.img\nrun\n")

	if findings := scan(t, root).Findings; len(findings) != 1 {
		t.Fatalf("findings = %#v, want one", findings)
	}
}

// A dot-directory carries tool state, not project source. A local conda env,
// .venv, .tox or .snakemake all ship shell scripts of their own, and a #DEP: in
// one of those would become a declaration this project has to pin.
func TestScanSkipsDotDirectories(t *testing.T) {
	root := t.TempDir()
	for _, dir := range []string{".venv/bin", ".snakemake", ".git", "scripts"} {
		if err := os.MkdirAll(filepath.Join(root, dir), 0o775); err != nil {
			t.Fatal(err)
		}
	}
	write := func(rel, dep string) {
		t.Helper()
		body := "#!/usr/bin/env bash\n#DEP: " + dep + "\n"
		if err := os.WriteFile(filepath.Join(root, rel), []byte(body), 0o664); err != nil {
			t.Fatal(err)
		}
	}
	write("scripts/a.sh", "real/1.0")
	write(".venv/bin/activate.sh", "phantom/9.9")
	write(".snakemake/x.sh", "phantom/9.9")
	write(".git/hook.sh", "phantom/9.9")

	result, err := Scan(root, ScanOptions{})
	if err != nil {
		t.Fatal(err)
	}
	if len(result.Requests) != 1 || result.Requests[0].Key != "real/1.0" {
		t.Fatalf("requests = %#v, want only real/1.0", result.Requests)
	}
}

// The root is the project, so it is read even when its own name is hidden.
func TestScanReadsAHiddenProjectRoot(t *testing.T) {
	root := filepath.Join(t.TempDir(), ".hidden-project")
	if err := os.MkdirAll(root, 0o775); err != nil {
		t.Fatal(err)
	}
	body := []byte("#!/usr/bin/env bash\n#DEP: real/1.0\n")
	if err := os.WriteFile(filepath.Join(root, "a.sh"), body, 0o664); err != nil {
		t.Fatal(err)
	}

	result, err := Scan(root, ScanOptions{})
	if err != nil {
		t.Fatal(err)
	}
	if len(result.Requests) != 1 {
		t.Fatalf("requests = %#v, want real/1.0", result.Requests)
	}
}

// A recipe's #DEP: are its artifact's build dependencies, already recorded in
// that artifact's provenance, not declarations the project mounts. $CNT_PREFIX
// is what identifies one, and comments are stripped before the check.
func TestScanSkipsBuildRecipes(t *testing.T) {
	root := t.TempDir()
	write(t, root, "tool.sh", "#DEP: build-only/9.9\ninstall -d \"$CNT_PREFIX/bin\"\n")
	write(t, root, "env.sh", "#DEP: build-only/8.8\ncp x \"${CNT_PREFIX}\"/lib\n")
	write(t, root, "run.sh", "# built into $CNT_PREFIX\n#DEP: real/1.0\n")

	result := scan(t, root)
	if got := keys(result); len(got) != 1 || got[0] != "real/1.0" {
		t.Fatalf("requests = %v, want real/1.0", got)
	}
	if len(result.Scripts) != 1 || result.Scripts[0] != "run.sh" {
		t.Errorf("scripts = %v, want run.sh alone", result.Scripts)
	}
}
