package build

import (
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/Justype/condatainer/catalog"
	"github.com/Justype/condatainer/internal/config"
	"github.com/Justype/condatainer/internal/image/producer"
)

// Every constructor sites the workspace itself, and they have disagreed before —
// a definition build once got .img from one and .sif from another. These pin the
// full path set each one produces, so plan/build-object.md §3 can route them all
// through one derivation and prove nothing moved.
//
// tmppath_test.go covers the helpers; this covers the constructors' use of them.

// wantPaths is the full path set for one build.
type wantPaths struct {
	baseRoot string
	root     string // producer-private root
	cntDir   string
	tmpImg   string // "" when the mode has no scratch image
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
	if b.ws.Overlay != want.tmpImg {
		t.Errorf("tmpOverlayPath = %q, want %q", b.ws.Overlay, want.tmpImg)
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

// expect builds the path set the helpers derive, which is what §3.1 claims
// workspaceFor reproduces.
func expect(name, root, ext, target string) wantPaths {
	ws := workspaceFor(name, root, ext)
	return wantPaths{
		baseRoot: ws.BaseRoot, root: ws.Root, cntDir: ws.CntDir,
		tmpImg: ws.Overlay, target: target, metaDir: ws.MetaDir, source: ws.Source,
	}
}

// The scratch image exists exactly when the mode uses one. This constructor used
// to name a .img unconditionally, so Overlay was set even in directory mode where
// nothing ever created it — the reason UsesImage now decides the mode.
func TestCondaConstructorPaths(t *testing.T) {
	for _, appTmpOverlay := range []bool{false, true} {
		name := map[bool]string{false: "dir mode", true: "ext3 mode"}[appTmpOverlay]
		t.Run(name, func(t *testing.T) {
			prev := config.Global.Build.AppTmpOverlay
			config.Global.Build.AppTmpOverlay = appTmpOverlay
			t.Cleanup(func() { config.Global.Build.AppTmpOverlay = prev })

			dir := t.TempDir()
			const module = "samtools/1.23.1"

			b, err := NewCondaObjectWithSource(module, "", dir, false)
			if err != nil {
				t.Fatalf("NewCondaObjectWithSource: %v", err)
			}
			ext := ""
			if appTmpOverlay {
				ext = ".img"
			}
			checkPaths(t, b, expect(module, b.ws.BaseRoot, ext,
				filepath.Join(dir, "samtools--1.23.1.sqf")))

			if b.ws.UsesImage() != appTmpOverlay {
				t.Errorf("UsesImage = %v, want %v", b.ws.UsesImage(), appTmpOverlay)
			}
			// An app stages inside the image whenever there is one.
			if want := !appTmpOverlay; b.ws.HostPayload() != want {
				t.Errorf("HostPayload = %v, want %v", b.ws.HostPayload(), want)
			}
			if b.spec.Image.Type != catalog.TypeApp {
				t.Errorf("type = %q, want app", b.spec.Image.Type)
			}
		})
	}
}

// The scratch extension is two decisions, not one: a definition always gets a
// SIF, a shell build gets .img only under use_tmp_overlay and otherwise runs in
// dir mode with no scratch image at all. Empty means "dir mode" — the string
// encodes it, which is what plan/build-object.md §3.2 replaces.
func TestExternalConstructorPaths(t *testing.T) {
	for _, tc := range []struct {
		name          string
		isDef         bool
		appTmpOverlay bool
		wantExt       string
		fileName      string
	}{
		{"shell script, dir mode", false, false, "", "demo.sh"},
		{"shell script, ext3 mode", false, true, ".img", "demo.sh"},
		{"definition, dir mode", true, false, ".sif", "demo.def"},
		{"definition, ext3 mode", true, true, ".sif", "demo.def"},
	} {
		t.Run(tc.name, func(t *testing.T) {
			prev := config.Global.Build.AppTmpOverlay
			config.Global.Build.AppTmpOverlay = tc.appTmpOverlay
			t.Cleanup(func() { config.Global.Build.AppTmpOverlay = prev })

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
			checkPaths(t, b, expect("demo", b.ws.BaseRoot, tc.wantExt, prefix+".sqf"))
		})
	}
}

// The one place a workspace is re-sited after construction: it must move every
// path together, never some. What the mode becomes is TestRetargetDropsExt3's
// job; this one checks the set stays internally consistent.
func TestRetargetMovesEveryPath(t *testing.T) {
	withAppTmpOverlayMode(t, true)
	dir := t.TempDir()
	const name = "grch38/genome/gencode"

	b := &BuildObject{
		spec: Spec{Image: ImageSpec{Name: name, Type: catalog.TypeApp}},
		ws:   workspaceFor(name, dir, appExt3ScratchExt(catalog.TypeApp)),
	}

	// data re-sites to the stable shared root rather than fast local scratch.
	b.spec.Image.Type = catalog.TypeData
	b.retargetWorkspace()

	if b.ws.BaseRoot == dir {
		t.Skip("data and app share a tmp root in this environment; nothing to re-site")
	}
	// Data drops the ext3 image, so the whole set is re-derived with no overlay.
	checkPaths(t, b, expect(name, b.ws.BaseRoot, "", ""))
}

// The layout itself, spelled out. Everything else in the package derives from
// workspaceFor now, so this is the one place the shape is stated rather than
// computed — a comparison against another derivation would only prove they agree.
func TestWorkspaceForLayout(t *testing.T) {
	for _, tc := range []struct {
		name, module, ext string
		wantParent        string // relative to base root
	}{
		{"conda", "samtools/1.23.1", ".img", "build_samtools_1.23.1"},
		{"definition", "ubuntu24/base", ".sif", "build_ubuntu24_base"},
		{"dir mode", "demo", "", "build_demo"},
		{"deep data name", "grch38/star/2.7.11b/gencode47-101", ".img",
			"build_grch38_star_2.7.11b_gencode47-101"},
	} {
		t.Run(tc.name, func(t *testing.T) {
			root := t.TempDir()
			ws := workspaceFor(tc.module, root, tc.ext)

			ownerRoot := filepath.Join(root, tc.wantParent, producer.Tag(producer.LocalInfo()))
			buildDir := filepath.Join(ownerRoot, "work")
			overlay := ""
			if tc.ext != "" {
				overlay = filepath.Join(ownerRoot, "rootfs"+tc.ext)
			}
			for _, c := range []struct{ got, want, field string }{
				{ws.BaseRoot, root, "BaseRoot"},
				{ws.Root, ownerRoot, "Root"},
				{ws.BuildDir, buildDir, "BuildDir"},
				{ws.CntDir, filepath.Join(buildDir, "cnt"), "CntDir"},
				{ws.TmpDir, filepath.Join(buildDir, "tmp"), "TmpDir"},
				{ws.MetaDir, filepath.Join(buildDir, ".cnt"), "MetaDir"},
				{ws.Source, filepath.Join(ownerRoot, sourceFileName(tc.module, tc.ext == ".sif")), "Source"},
				{ws.Overlay, overlay, "Overlay"},
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
	first := workspaceForOwner("hello/1.0", root, ".img",
		BuildLockInfo{Runner: "local", Node: "node-a", PID: 10})
	second := workspaceForOwner("hello/1.0", root, ".img",
		BuildLockInfo{Runner: "local", Node: "node-b", PID: 10})
	for _, pair := range [][2]string{
		{first.Root, second.Root},
		{first.Source, second.Source},
		{first.BuildDir, second.BuildDir},
		{first.Overlay, second.Overlay},
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
		ws:         workspaceFor("hello/1.0", root, ""),
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

// withAppTmpOverlayMode sets build.app_tmp_overlay for one test.
func withAppTmpOverlayMode(t *testing.T, on bool) {
	t.Helper()
	prev := config.Global.Build.AppTmpOverlay
	config.Global.Build.AppTmpOverlay = on
	t.Cleanup(func() { config.Global.Build.AppTmpOverlay = prev })
}

// build.app_tmp_overlay is an app-build optimisation: the ext3 overlay keeps a conda
// environment's thousands of small files off the host's inode budget. No other
// type may be put inside one — data stages on the host, and os and base are
// definition builds where apptainer owns the rootfs.
func TestExt3IsAppOnly(t *testing.T) {
	for _, on := range []bool{false, true} {
		t.Run(map[bool]string{false: "use_tmp_overlay off", true: "use_tmp_overlay on"}[on], func(t *testing.T) {
			withAppTmpOverlayMode(t, on)

			wantApp := ""
			if on {
				wantApp = ".img"
			}
			for _, tc := range []struct {
				typ  catalog.Type
				want string
			}{
				{catalog.TypeApp, wantApp},
				{catalog.TypeData, ""}, // never, whatever the config says
				{catalog.TypeOS, ""},
				{catalog.TypeBase, ""},
			} {
				if got := appExt3ScratchExt(tc.typ); got != tc.want {
					t.Errorf("appExt3ScratchExt(%s) = %q, want %q", tc.typ, got, tc.want)
				}
				// The payload lands on the host exactly when there is no image.
				ws := workspaceFor("demo/1.0", t.TempDir(), tc.want)
				if got, want := ws.HostPayload(), tc.want == ""; got != want {
					t.Errorf("HostPayload for %s = %v, want %v", tc.typ, got, want)
				}
			}
		})
	}
}

// The same rule through a constructor, which is where it actually bites.
func TestDataBuildNeverGetsExt3(t *testing.T) {
	withAppTmpOverlayMode(t, true) // the mode that would otherwise hand it an overlay

	dir := t.TempDir()
	src := filepath.Join(dir, "mydata.sh")
	if err := os.WriteFile(src, []byte("#!/usr/bin/env bash\n#TYPE:data\n#DESC:demo\necho hi\n"), 0o644); err != nil {
		t.Fatal(err)
	}

	b, err := FromExternalSource(t.Context(), filepath.Join(dir, "demo"), src, false, dir, false)
	if err != nil {
		t.Fatalf("FromExternalSource: %v", err)
	}
	if b.spec.Image.Type != catalog.TypeData {
		t.Fatalf("type = %q, want data — #TYPE: was ignored", b.spec.Image.Type)
	}
	if b.ws.UsesImage() {
		t.Errorf("data build got an ext3 scratch image at %q", b.ws.Overlay)
	}
	if !b.ws.HostPayload() {
		t.Error("data build does not stage its payload on the host")
	}
}

// A definition is unaffected either way: apptainer owns the rootfs, so the .sif
// is the product and use_tmp_overlay has nothing to say about it.
func TestDefinitionIgnoresExt3Mode(t *testing.T) {
	for _, on := range []bool{false, true} {
		withAppTmpOverlayMode(t, on)
		dir := t.TempDir()
		src := filepath.Join(dir, "demo.def")
		if err := os.WriteFile(src, []byte("Bootstrap: docker\nFrom: alpine:3.19\n"), 0o644); err != nil {
			t.Fatal(err)
		}
		b, err := FromExternalSource(t.Context(), filepath.Join(dir, "demo"), src, true, dir, false)
		if err != nil {
			t.Fatalf("FromExternalSource: %v", err)
		}
		if filepath.Ext(b.ws.Overlay) != ".sif" {
			t.Errorf("use_tmp_overlay=%v gave a definition %q, want a .sif", on, b.ws.Overlay)
		}
		if b.ws.HostPayload() {
			t.Errorf("use_tmp_overlay=%v: definition staged on host; apptainer owns the rootfs", on)
		}
	}
}

// The guess is corrected after the recipe is read. app -> data must drop the
// image, not merely move it: carrying the old extension across was the bug.
func TestRetargetDropsExt3WhenTypeBecomesData(t *testing.T) {
	withAppTmpOverlayMode(t, true)

	dir := t.TempDir()
	const name = "grch38/genome/gencode"
	b := &BuildObject{
		spec: Spec{Image: ImageSpec{Name: name, Type: catalog.TypeApp}},
		ws:   workspaceFor(name, dir, appExt3ScratchExt(catalog.TypeApp)),
	}
	if !b.ws.UsesImage() {
		t.Fatal("fixture should start with an ext3 image")
	}

	b.spec.Image.Type = catalog.TypeData
	b.retargetWorkspace()

	if b.ws.UsesImage() {
		t.Errorf("retarget kept an ext3 image at %q after the type became data", b.ws.Overlay)
	}
	if !b.ws.HostPayload() {
		t.Error("retarget left the payload inside a nonexistent image")
	}
}
