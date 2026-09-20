package build

import (
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/Justype/condatainer/catalog"
	"github.com/Justype/condatainer/internal/image/producer"
)

// Every constructor sites the workspace itself, and they have disagreed before —
// a definition build once got .img from one and .sif from another. These pin the
// full path set each one produces, so a build routes them all through one
// derivation and nothing moves independently.
//
// tmppath_test.go covers the helpers; this covers the constructors' use of them.

// wantPaths is the full path set for one build.
type wantPaths struct {
	baseRoot string
	root     string // producer-private root
	cntDir   string
	target   string
	metaDir  string
	source   string
}

func checkPaths(t *testing.T, b *BuildObject, want wantPaths) {
	t.Helper()
	if b.ws.BaseRoot != want.baseRoot {
		t.Errorf("base root = %q, want %q", b.ws.BaseRoot, want.baseRoot)
	}
	if b.ws.Root != want.root {
		t.Errorf("tmpDir = %q, want %q", b.ws.Root, want.root)
	}
	if b.ws.CntDir != want.cntDir {
		t.Errorf("cntDirPath = %q, want %q", b.ws.CntDir, want.cntDir)
	}
	if b.tgt.Path != want.target {
		t.Errorf("target Path = %q, want %q", b.tgt.Path, want.target)
	}
	if b.ws.MetaDir != want.metaDir {
		t.Errorf("MetaDir = %q, want %q", b.ws.MetaDir, want.metaDir)
	}
	if b.ws.Source != want.source {
		t.Errorf("Source = %q, want %q", b.ws.Source, want.source)
	}
	if wantDir := filepath.Join(filepath.Dir(want.cntDir), "tmp"); b.ws.TmpDir != wantDir {
		t.Errorf("TmpDir = %q, want %q", b.ws.TmpDir, wantDir)
	}
	if wantDir := filepath.Dir(want.cntDir); b.ws.BuildDir != wantDir {
		t.Errorf("BuildDir = %q, want %q", b.ws.BuildDir, wantDir)
	}
}

// expect builds the path set the helpers derive, which is what workspaceFor
// must reproduce.
func expect(name, root, target string, isDef bool) wantPaths {
	ws := workspaceFor(name, root, isDef)
	return wantPaths{
		baseRoot: ws.BaseRoot, root: ws.Root, cntDir: ws.CntDir,
		target: target, metaDir: ws.MetaDir, source: ws.Source,
	}
}

func TestCondaConstructorPaths(t *testing.T) {
	dir := t.TempDir()
	const module = "samtools/1.23.1"

	b, err := NewCondaObjectWithSource(module, "", dir, false)
	if err != nil {
		t.Fatalf("NewCondaObjectWithSource: %v", err)
	}
	checkPaths(t, b, expect(module, b.ws.BaseRoot, filepath.Join(dir, "samtools--1.23.1.sqf"), false))

	if b.spec.Image.Type != catalog.TypeApp {
		t.Errorf("type = %q, want app", b.spec.Image.Type)
	}
}

// A shell build and a definition site their workspaces the same way; only the
// recipe's extension and the sandbox differ.
func TestExternalConstructorPaths(t *testing.T) {
	for _, tc := range []struct {
		name     string
		isDef    bool
		fileName string
	}{
		{"shell script", false, "demo.sh"},
		{"definition", true, "demo.def"},
	} {
		t.Run(tc.name, func(t *testing.T) {
			dir := t.TempDir()
			prefix := filepath.Join(dir, "demo")
			src := filepath.Join(dir, tc.fileName)
			if err := os.WriteFile(src, []byte("#!/usr/bin/env bash\n#DESC:demo\necho hi\n"), 0o644); err != nil {
				t.Fatal(err)
			}

			b, err := FromExternalSource(t.Context(), prefix, src, tc.isDef, dir, false)
			if err != nil {
				t.Fatalf("FromExternalSource: %v", err)
			}
			checkPaths(t, b, expect("demo", b.ws.BaseRoot, prefix+".sqf", tc.isDef))
		})
	}
}

// The one place a workspace is re-sited after construction: it must move every
// path together, never some.
func TestRetargetMovesEveryPath(t *testing.T) {
	dir := t.TempDir()
	const name = "grch38/genome/gencode"

	b := &BuildObject{
		spec: Spec{Image: ImageSpec{Name: name, Type: catalog.TypeApp}},
		ws:   workspaceFor(name, dir, false),
	}

	// data re-sites to the stable shared root rather than fast local scratch.
	b.spec.Image.Type = catalog.TypeData
	b.retargetWorkspace()

	if b.ws.BaseRoot == dir {
		t.Skip("data and app share a tmp root in this environment; nothing to re-site")
	}
	checkPaths(t, b, expect(name, b.ws.BaseRoot, "", false))
}

// The layout itself, spelled out. Everything else in the package derives from
// workspaceFor now, so this is the one place the shape is stated rather than
// computed — a comparison against another derivation would only prove they agree.
func TestWorkspaceForLayout(t *testing.T) {
	for _, tc := range []struct {
		name, module string
		isDef        bool
		wantParent   string // relative to base root
	}{
		{"conda", "samtools/1.23.1", false, "build_samtools_1.23.1"},
		{"definition", "ubuntu24/base", true, "build_ubuntu24_base"},
		{"plain name", "demo", false, "build_demo"},
		{"deep data name", "grch38/star/2.7.11b/gencode47-101", false,
			"build_grch38_star_2.7.11b_gencode47-101"},
	} {
		t.Run(tc.name, func(t *testing.T) {
			root := t.TempDir()
			ws := workspaceFor(tc.module, root, tc.isDef)

			ownerRoot := filepath.Join(root, tc.wantParent, producer.Tag(producer.LocalInfo()))
			buildDir := filepath.Join(ownerRoot, "work")
			sandbox := ""
			if tc.isDef {
				sandbox = filepath.Join(ownerRoot, "rootfs")
			}
			for _, c := range []struct{ got, want, field string }{
				{ws.BaseRoot, root, "BaseRoot"},
				{ws.Root, ownerRoot, "Root"},
				{ws.BuildDir, buildDir, "BuildDir"},
				{ws.CntDir, filepath.Join(buildDir, "cnt"), "CntDir"},
				{ws.TmpDir, filepath.Join(buildDir, "tmp"), "TmpDir"},
				{ws.MetaDir, filepath.Join(buildDir, ".cnt"), "MetaDir"},
				{ws.Source, filepath.Join(ownerRoot, sourceFileName(tc.module, tc.isDef)), "Source"},
				{ws.Sandbox, sandbox, "Sandbox"},
			} {
				if c.got != c.want {
					t.Errorf("%s = %q, want %q", c.field, c.got, c.want)
				}
			}
		})
	}
}

func TestWorkspaceOwnersDoNotShareTemporaryPaths(t *testing.T) {
	root := t.TempDir()
	first := workspaceForOwner("hello/1.0", root, false,
		BuildLockInfo{Runner: "local", Node: "node-a", PID: 10})
	second := workspaceForOwner("hello/1.0", root, false,
		BuildLockInfo{Runner: "local", Node: "node-b", PID: 10})
	for _, pair := range [][2]string{
		{first.Root, second.Root},
		{first.Source, second.Source},
		{first.BuildDir, second.BuildDir},
	} {
		if pair[0] == pair[1] {
			t.Fatalf("two owners share temporary path %q", pair[0])
		}
	}
}

func TestAdoptWorkspaceMovesMaterializedRecipe(t *testing.T) {
	root := t.TempDir()
	b := &BuildObject{
		spec:       Spec{Image: ImageSpec{Name: "hello/1.0"}},
		ws:         workspaceFor("hello/1.0", root, false),
		tempSource: true,
	}
	b.buildSource = b.ws.Source
	if err := os.MkdirAll(filepath.Dir(b.buildSource), 0o755); err != nil {
		t.Fatal(err)
	}
	if err := os.WriteFile(b.buildSource, []byte("echo hello\n"), 0o644); err != nil {
		t.Fatal(err)
	}

	owner := BuildLockInfo{Runner: "slurm", JobID: "42"}
	oldRoot := b.ws.Root
	if err := b.adoptWorkspace(owner); err != nil {
		t.Fatal(err)
	}
	if !strings.Contains(b.ws.Root, "slurm-42") {
		t.Fatalf("adopted root = %q", b.ws.Root)
	}
	if b.buildSource != b.ws.Source {
		t.Fatalf("build source = %q, workspace source = %q", b.buildSource, b.ws.Source)
	}
	if _, err := os.Stat(b.buildSource); err != nil {
		t.Fatalf("materialized recipe was not moved: %v", err)
	}
	if _, err := os.Stat(oldRoot); !os.IsNotExist(err) {
		t.Fatalf("old owner workspace remains: %v", err)
	}
}
