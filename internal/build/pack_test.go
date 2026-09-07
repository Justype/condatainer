package build

import (
	"encoding/json"
	"fmt"
	"os"
	"os/exec"
	"path/filepath"
	"strings"
	"testing"
	"time"

	"github.com/Justype/condatainer/catalog"
	"github.com/Justype/condatainer/internal/artifact/meta"
	"github.com/Justype/condatainer/internal/config"
	"github.com/Justype/condatainer/internal/toolpath"
	"github.com/Justype/condatainer/internal/utils"
)

// withProvisionedLibexec points the scratch tier at a fresh temp directory and
// drops stub bin/apptainer and bin/micromamba there, satisfying libexec.Dir's
// marker check. Mirrors internal/toolpath's own test helper of the same name,
// duplicated here rather than imported: libexec's test helpers are unexported.
func withProvisionedLibexec(t *testing.T) {
	t.Helper()
	scratch := filepath.Join(t.TempDir(), "condatainer")
	t.Setenv("SCRATCH", filepath.Dir(scratch))
	t.Setenv("XDG_DATA_HOME", "")
	t.Setenv("CNT_EXTRA_ROOT", "")
	t.Setenv("CNT_ROOT", "")
	config.InitDataPaths()
	t.Cleanup(config.InitDataPaths)

	bin := filepath.Join(scratch, "libexec", "bin")
	if err := utils.MkdirAllShared(bin); err != nil {
		t.Fatalf("failed to create stub bin dir: %v", err)
	}
	for _, name := range []string{"apptainer", "micromamba"} {
		if err := os.WriteFile(filepath.Join(bin, name), []byte("#!/bin/sh\n"), 0o755); err != nil {
			t.Fatalf("failed to write stub %s: %v", name, err)
		}
	}
}

// newPackObject builds the minimum BuildObject the packer reads: a workspace
// layout and a Spec complete enough to render both metadata documents.
func newPackObject(t *testing.T, typ catalog.Type) *BuildObject {
	t.Helper()
	tmpDir := t.TempDir()
	name := "samtools/1.21"
	b := &BuildObject{
		ws: workspaceFor(name, tmpDir, appExt3ScratchExt(typ), false),
		spec: Spec{
			Image:  ImageSpec{Name: name, Type: typ, Description: "SAMtools", Prefix: meta.Prefix(name, typ)},
			Source: SourceSpec{Conda: &CondaSource{Package: &CondaPackage{Name: "samtools", Version: "1.21"}}},
		},
	}
	return b
}

// withAppTmpOverlay flips the global build mode for one test and restores it. The
// pack settings come along because the zero value is not a legal mksquashfs
// argument, and a test that runs the real command needs them.
func withAppTmpOverlay(t *testing.T, on bool) {
	t.Helper()
	prev := config.Global.Build
	config.Global.Build.AppTmpOverlay = on
	if config.Global.Build.BlockSize == "" {
		config.Global.Build.BlockSize = config.DefaultBlockSize
	}
	if config.Global.Build.DataBlockSize == "" {
		config.Global.Build.DataBlockSize = config.DefaultDataBlockSize
	}
	t.Cleanup(func() { config.Global.Build = prev })
}

// Both documents are staged, since an image carrying only one of them is either
// unmountable or unexplainable.
func TestStageMetadataWritesBothDocuments(t *testing.T) {
	b := newPackObject(t, catalog.TypeApp)

	dir, err := stageMetadata(t.Context(), b)
	if err != nil {
		t.Fatalf("stageMetadata: %v", err)
	}
	if filepath.Base(dir) != meta.DirName {
		t.Errorf("staged into %q, want a directory named %q", dir, meta.DirName)
	}

	runtime, err := meta.MarshalRuntime(b.Runtime())
	if err != nil {
		t.Fatal(err)
	}
	data, err := os.ReadFile(filepath.Join(dir, meta.RuntimeFileName))
	if err != nil {
		t.Fatalf("reading staged %s: %v", meta.RuntimeFileName, err)
	}
	if string(data) != string(runtime) {
		t.Errorf("staged %s = %s, want %s", meta.RuntimeFileName, data, runtime)
	}

	// The manifest is not compared byte for byte: staging stamps build.created,
	// which Manifest() deliberately does not carry. Check the stamp, then clear
	// it and hold the rest to the projection.
	data, err = os.ReadFile(filepath.Join(dir, meta.FileName))
	if err != nil {
		t.Fatalf("reading staged %s: %v", meta.FileName, err)
	}
	var staged meta.Manifest
	if err := json.Unmarshal(data, &staged); err != nil {
		t.Fatalf("staged manifest is not JSON: %v", err)
	}
	if staged.Build.Created.IsZero() {
		t.Error("staging did not stamp build.created")
	} else if since := time.Since(staged.Build.Created); since < 0 || since > time.Minute {
		t.Errorf("build.created = %s, %s away from now", staged.Build.Created, since)
	}
	staged.Build.Created = time.Time{}

	want, err := meta.MarshalManifest(b.Manifest())
	if err != nil {
		t.Fatal(err)
	}
	got, err := meta.MarshalManifest(staged)
	if err != nil {
		t.Fatal(err)
	}
	if string(got) != string(want) {
		t.Errorf("staged %s = %s, want %s", meta.FileName, got, want)
	}
}

// Metadata that would not survive a read back must stop the build while there is
// still no image, not after one is installed and read back as broken.
func TestStageMetadataRejectsInvalidMetadata(t *testing.T) {
	b := newPackObject(t, catalog.TypeApp)
	b.spec.Image.Prefix = "" // an app with no prefix cannot be resolved at load time

	if _, err := stageMetadata(t.Context(), b); err == nil {
		t.Fatal("staged metadata that does not validate")
	}
	if _, err := os.Stat(b.ws.MetaDir); err == nil {
		t.Error("invalid metadata left a staged directory behind")
	}
}

// The staged directory has to reach mksquashfs still called .cnt, in both modes,
// because mksquashfs names an archive root after its source's basename.
func TestSquashfsSourcesCarryMetaDirName(t *testing.T) {
	tests := []struct {
		name      string
		appTmpOvl bool
		sourceDir string
		ext3Base  string
	}{
		{name: "dir mode", appTmpOvl: false, sourceDir: "/host/build/cnt"},
		{name: "ext3 payload in image", appTmpOvl: true, ext3Base: "/mnt/upper"},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			withAppTmpOverlay(t, tt.appTmpOvl)
			b := newPackObject(t, catalog.TypeApp)
			metaDir := b.ws.MetaDir

			sources, keepAsDirectory := packSources(b, tt.sourceDir, metaDir, tt.ext3Base)

			wantFirst := tt.sourceDir
			if tt.appTmpOvl {
				wantFirst = filepath.Join(tt.ext3Base, "cnt")
			}
			if len(sources) != 2 || sources[0] != wantFirst || sources[1] != metaDir {
				t.Errorf("sources = %v, want [%q %q]", sources, wantFirst, metaDir)
			}
			if !keepAsDirectory {
				t.Error("keepAsDirectory = false, want true")
			}
			if filepath.Base(metaDir) != meta.DirName {
				t.Errorf("metadata source %q does not end in %q", metaDir, meta.DirName)
			}
		})
	}
}

// An empty metaDir packs the payload alone, which is what keeps the packer
// usable for anything that has no metadata to add.
func TestSquashfsWithoutMetaDirPacksPayloadOnly(t *testing.T) {
	withAppTmpOverlay(t, false)
	b := newPackObject(t, catalog.TypeApp)

	sources, keepAsDirectory := packSources(b, b.ws.CntDir, "", "")

	if len(sources) != 1 || sources[0] != b.ws.CntDir {
		t.Errorf("sources = %v, want [%q]", sources, b.ws.CntDir)
	}
	if !keepAsDirectory {
		t.Error("keepAsDirectory = false, want true")
	}
}

// The two-source layout is only correct if mksquashfs actually puts both
// documents under /.cnt. Generating the command and asserting on its text cannot
// show that, so this runs the real thing and reads them back.
func TestPackedImageMetadataIsReadable(t *testing.T) {
	if _, err := exec.LookPath("mksquashfs"); err != nil {
		t.Skip("mksquashfs not available")
	}
	withAppTmpOverlay(t, false)

	b := newPackObject(t, catalog.TypeApp)
	payload := filepath.Join(b.ws.CntDir, b.spec.Image.Name, "bin")
	if err := os.MkdirAll(payload, 0o755); err != nil {
		t.Fatal(err)
	}
	if err := os.WriteFile(filepath.Join(payload, "samtools"), []byte("#!/bin/sh\n"), 0o755); err != nil {
		t.Fatal(err)
	}

	metaDir, err := stageMetadata(t.Context(), b)
	if err != nil {
		t.Fatalf("stageMetadata: %v", err)
	}

	out := filepath.Join(t.TempDir(), "samtools--1.21.sqf")
	if err := createSquashfs(t.Context(), b, false, b.ws.CntDir, metaDir, out); err != nil {
		t.Fatalf("createSquashfs: %v", err)
	}

	manifest, err := meta.ReadManifest(out)
	if err != nil {
		t.Fatalf("meta.ReadManifest on the packed image: %v", err)
	}
	if manifest.Name != b.spec.Image.Name || manifest.Type != catalog.TypeApp {
		t.Errorf("manifest = %+v, want name %q type %q", manifest, b.spec.Image.Name, catalog.TypeApp)
	}
	if manifest.BuildType != BuildTypeConda {
		t.Errorf("build_type = %q, want %q", manifest.BuildType, BuildTypeConda)
	}

	rt, err := meta.ReadRuntime(out)
	if err != nil {
		t.Fatalf("meta.ReadRuntime on the packed image: %v", err)
	}
	if rt.Prefix != meta.Prefix(b.spec.Image.Name, catalog.TypeApp) {
		t.Errorf("prefix = %q", rt.Prefix)
	}

	// The payload has to survive beside the metadata, not be replaced by it.
	list, err := exec.CommandContext(t.Context(), "unsquashfs", "-l", out).Output()
	if err != nil {
		t.Fatalf("unsquashfs -l: %v", err)
	}
	for _, want := range []string{
		"squashfs-root/cnt/samtools/1.21/bin/samtools",
		"squashfs-root" + meta.Path,
		"squashfs-root" + meta.RuntimePath,
	} {
		if !strings.Contains(string(list), want) {
			t.Errorf("packed image is missing %s:\n%s", want, list)
		}
	}
}

// The build's own scratch directory sits beside the payload and the metadata, so
// the packer must not sweep it into the image.
func TestPackedImageExcludesBuildScratch(t *testing.T) {
	if _, err := exec.LookPath("mksquashfs"); err != nil {
		t.Skip("mksquashfs not available")
	}
	withAppTmpOverlay(t, false)

	b := newPackObject(t, catalog.TypeApp)
	if err := os.MkdirAll(filepath.Join(b.ws.CntDir, b.spec.Image.Name), 0o755); err != nil {
		t.Fatal(err)
	}
	scratch := b.ws.TmpDir
	if err := os.MkdirAll(scratch, 0o755); err != nil {
		t.Fatal(err)
	}
	if err := os.WriteFile(filepath.Join(scratch, "huge.pkg"), []byte("x"), 0o644); err != nil {
		t.Fatal(err)
	}

	metaDir, err := stageMetadata(t.Context(), b)
	if err != nil {
		t.Fatal(err)
	}
	out := filepath.Join(t.TempDir(), "out.sqf")
	if err := createSquashfs(t.Context(), b, false, b.ws.CntDir, metaDir, out); err != nil {
		t.Fatalf("createSquashfs: %v", err)
	}

	list, err := exec.CommandContext(t.Context(), "unsquashfs", "-l", out).Output()
	if err != nil {
		t.Fatal(err)
	}
	if strings.Contains(string(list), "huge.pkg") {
		t.Errorf("build scratch leaked into the image:\n%s", list)
	}
}

// TestPackedImageFromScratchOverlayReadsThroughFuse2fs exercises the ext3-mode
// pack path end to end: a real .img, populated via debugfs the way an actual
// conda install leaves it, packed by mounting it with fuse2fs — no container,
// no Apptainer, the same apptainer-free mechanism internal/image/freeze uses.
func TestPackedImageFromScratchOverlayReadsThroughFuse2fs(t *testing.T) {
	for _, name := range []string{"mksquashfs", "unsquashfs", "mke2fs", "debugfs", "fuse2fs", "unshare"} {
		if _, err := exec.LookPath(name); err != nil {
			t.Skipf("%s not available", name)
		}
	}
	withAppTmpOverlay(t, true)
	config.Global.Build.AppTmpOverlaySizeMB = 64

	b := newPackObject(t, catalog.TypeApp)
	if err := b.CreateTmpOverlay(t.Context(), false); err != nil {
		t.Fatalf("CreateTmpOverlay: %v", err)
	}

	hostFile := filepath.Join(t.TempDir(), "samtools")
	if err := os.WriteFile(hostFile, []byte("#!/bin/sh\n"), 0o755); err != nil {
		t.Fatal(err)
	}
	debugfsPath, err := toolpath.Resolve("debugfs")
	if err != nil {
		t.Fatal(err)
	}
	script := "cd upper\nmkdir cnt\ncd cnt\nmkdir samtools\ncd samtools\nmkdir 1.21\ncd 1.21\nmkdir bin\ncd bin\n" +
		fmt.Sprintf("write %s samtools\nquit\n", hostFile)
	cmd := exec.CommandContext(t.Context(), debugfsPath, "-w", b.ws.Overlay)
	cmd.Stdin = strings.NewReader(script)
	if out, err := cmd.CombinedOutput(); err != nil {
		t.Fatalf("debugfs inject: %v\n%s", err, out)
	}

	metaDir, err := stageMetadata(t.Context(), b)
	if err != nil {
		t.Fatalf("stageMetadata: %v", err)
	}

	out := filepath.Join(t.TempDir(), "samtools--1.21.sqf")
	if err := createSquashfs(t.Context(), b, false, "/cnt", metaDir, out); err != nil {
		t.Fatalf("createSquashfs: %v", err)
	}

	list, err := exec.CommandContext(t.Context(), "unsquashfs", "-l", out).Output()
	if err != nil {
		t.Fatalf("unsquashfs -l: %v", err)
	}
	for _, want := range []string{
		"squashfs-root/cnt/samtools/1.21/bin/samtools",
		"squashfs-root" + meta.Path,
	} {
		if !strings.Contains(string(list), want) {
			t.Errorf("packed image is missing %s:\n%s", want, list)
		}
	}
}

func TestSquashfsShowsProgressWithoutFinalStatistics(t *testing.T) {
	script := squashfsScript("mksquashfs", []string{"/cnt"}, "/images/out.sqf", 2, "128k", "-comp zstd", true)
	if !strings.Contains(script, " -quiet ") {
		t.Fatalf("mksquashfs command does not suppress final statistics:\n%s", script)
	}
	if strings.Contains(script, "-no-progress") {
		t.Fatalf("mksquashfs command suppresses progress:\n%s", script)
	}
}

// The conda install phase must not know where the image lands: that is what lets
// a cancelled install leave nothing next to the installed images.
func TestCondaInstallDoesNotPackOrTouchTarget(t *testing.T) {
	withAppTmpOverlay(t, false)
	withProvisionedLibexec(t)
	b := newPackObject(t, catalog.TypeApp)
	b.buildType = BuildTypeConda
	b.packageName, b.packageVersion = "samtools", "1.21"
	b.tgt = targetFor(filepath.Join(t.TempDir(), "images", "samtools--1.21.sqf"))

	opts, err := b.condaInstallExecOpts()
	if err != nil {
		t.Fatalf("condaInstallExecOpts: %v", err)
	}

	script := strings.Join(opts.Command, " ")
	if strings.Contains(script, "mksquashfs") {
		t.Errorf("install phase still packs:\n%s", script)
	}
	if !strings.Contains(script, "micromamba create") {
		t.Errorf("install phase does not install:\n%s", script)
	}
	if strings.Contains(strings.Join(opts.BindPaths, " "), filepath.Dir(b.tgt.Path)) {
		t.Errorf("install phase binds the images dir: %v", opts.BindPaths)
	}
}

// Conda and script must agree on where the payload is, or one writes to the
// image while the other packs the host directory.
func TestCondaInstallPayloadMatchesPackSource(t *testing.T) {
	for _, appTmpOverlay := range []bool{false, true} {
		t.Run(map[bool]string{false: "dir mode", true: "ext3 mode"}[appTmpOverlay], func(t *testing.T) {
			withAppTmpOverlay(t, appTmpOverlay)
			withProvisionedLibexec(t)
			b := newPackObject(t, catalog.TypeApp)
			b.buildType = BuildTypeConda
			b.packageName, b.packageVersion = "samtools", "1.21"

			opts, err := b.condaInstallExecOpts()
			if err != nil {
				t.Fatal(err)
			}
			installsToHost := strings.Contains(strings.Join(opts.BindPaths, " "), b.ws.CntDir+":/cnt")

			if installsToHost != b.ws.HostPayload() {
				t.Errorf("install writes to host = %v but packOutput reads host = %v",
					installsToHost, b.ws.HostPayload())
			}
		})
	}
}
