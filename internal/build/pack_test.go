package build

import (
	"os"
	"os/exec"
	"path/filepath"
	"strings"
	"testing"

	"github.com/Justype/condatainer/catalog"
	"github.com/Justype/condatainer/internal/config"
	"github.com/Justype/condatainer/internal/image/meta"
)

// newPackObject builds the minimum BuildObject the packer reads: a workspace
// layout and a Spec complete enough to render a manifest.
func newPackObject(t *testing.T, typ catalog.Type) *BuildObject {
	t.Helper()
	tmpDir := t.TempDir()
	name := "samtools/1.21"
	b := &BuildObject{
		ws: workspaceFor(name, tmpDir, appExt3ScratchExt(typ)),
		spec: Spec{
			Image:   ImageSpec{Name: name, Type: typ, Description: "SAMtools"},
			Source:  SourceSpec{Conda: &CondaSource{Package: &CondaPackage{Name: "samtools", Version: "1.21"}}},
			Runtime: meta.Runtime{Prefix: meta.Prefix(name, typ)},
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

func TestStageMetadataWritesManifest(t *testing.T) {
	b := newPackObject(t, catalog.TypeApp)

	dir, err := stageMetadata(t.Context(), b)
	if err != nil {
		t.Fatalf("stageMetadata: %v", err)
	}
	if filepath.Base(dir) != meta.DirName {
		t.Errorf("staged into %q, want a directory named %q", dir, meta.DirName)
	}

	data, err := os.ReadFile(filepath.Join(dir, meta.FileName))
	if err != nil {
		t.Fatalf("reading staged manifest: %v", err)
	}
	want, err := meta.Marshal(b.Manifest())
	if err != nil {
		t.Fatal(err)
	}
	if string(data) != string(want) {
		t.Errorf("staged manifest = %s, want %s", data, want)
	}
}

// A manifest that would not survive meta.Read must stop the build while there is
// still no image, not after one is installed and read back as broken.
func TestStageMetadataRejectsInvalidManifest(t *testing.T) {
	b := newPackObject(t, catalog.TypeApp)
	b.spec.Runtime.Prefix = "" // an app with no prefix cannot be resolved at load time
	b.spec.Image.Name = ""

	if _, err := stageMetadata(t.Context(), b); err == nil {
		t.Fatal("staged a manifest that does not validate")
	}
	if _, err := os.Stat(b.ws.MetaDir); err == nil {
		t.Error("invalid manifest left a staged directory behind")
	}
}

// The staged directory has to reach mksquashfs still called .cnt, in both modes,
// because mksquashfs names an archive root after its source's basename.
func TestSquashfsSourcesCarryMetaDirName(t *testing.T) {
	tests := []struct {
		name         string
		appTmpOvl    bool
		sourceDir    string
		wantMetaArg  string
		wantBindHost bool
	}{
		{name: "dir mode", appTmpOvl: false, sourceDir: "/host/build/cnt", wantMetaArg: "", wantBindHost: true},
		{name: "ext3 payload in image", appTmpOvl: true, sourceDir: "/cnt", wantMetaArg: metaMountPath},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			withAppTmpOverlay(t, tt.appTmpOvl)
			b := newPackObject(t, catalog.TypeApp)
			metaDir := b.ws.MetaDir

			script, _, binds := buildSquashfsOpts(b, false, tt.sourceDir, metaDir, "/images/out.sqf")

			wantArg := tt.wantMetaArg
			if wantArg == "" {
				wantArg = metaDir
			}
			if !strings.Contains(script, tt.sourceDir+" "+wantArg+" ") {
				t.Errorf("sources are not %q then %q:\n%s", tt.sourceDir, wantArg, script)
			}
			if filepath.Base(wantArg) != meta.DirName {
				t.Errorf("metadata source %q does not end in %q", wantArg, meta.DirName)
			}

			joined := strings.Join(binds, " ")
			switch {
			case tt.wantBindHost && !strings.Contains(joined, metaDir):
				t.Errorf("host metadata dir not bound: %v", binds)
			case !tt.wantBindHost && !strings.Contains(joined, metaDir+":"+metaMountPath):
				t.Errorf("metadata dir not bound to %s: %v", metaMountPath, binds)
			}
		})
	}
}

// An empty metaDir packs the payload alone, which is what keeps the packer
// usable for anything that has no manifest to add.
func TestSquashfsWithoutMetaDirPacksPayloadOnly(t *testing.T) {
	withAppTmpOverlay(t, false)
	b := newPackObject(t, catalog.TypeApp)

	script, _, binds := buildSquashfsOpts(b, false, b.ws.CntDir, "", "/images/out.sqf")

	if strings.Contains(script, meta.DirName) {
		t.Errorf("packed a metadata source that was not staged:\n%s", script)
	}
	if strings.Contains(strings.Join(binds, " "), meta.DirName) {
		t.Errorf("bound a metadata dir that was not staged: %v", binds)
	}
}

// The two-source layout is only correct if mksquashfs actually puts the manifest
// at meta.Path. Generating the command and asserting on its text cannot show
// that, so this runs the real thing and reads it back through meta.Read.
func TestPackedImageManifestIsReadable(t *testing.T) {
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
	script, _, _ := buildSquashfsOpts(b, false, b.ws.CntDir, metaDir, out)

	// The script is plain bash and the container only supplies mksquashfs, so
	// running it directly tests the same command the packer would issue.
	cmd := exec.CommandContext(t.Context(), "/bin/bash", "-c", script)
	if output, err := cmd.CombinedOutput(); err != nil {
		t.Fatalf("mksquashfs failed: %v\n%s", err, output)
	}

	got, err := meta.Read(out)
	if err != nil {
		t.Fatalf("meta.Read on the packed image: %v", err)
	}
	if got.Name != b.spec.Image.Name || got.Type != catalog.TypeApp {
		t.Errorf("manifest = %+v, want name %q type %q", got, b.spec.Image.Name, catalog.TypeApp)
	}
	if got.BuildType != BuildTypeConda.String() {
		t.Errorf("build_type = %q, want %q", got.BuildType, BuildTypeConda)
	}
	if got.Runtime.Prefix != meta.Prefix(b.spec.Image.Name, catalog.TypeApp) {
		t.Errorf("prefix = %q", got.Runtime.Prefix)
	}

	// The payload has to survive beside the manifest, not be replaced by it.
	list, err := exec.CommandContext(t.Context(), "unsquashfs", "-l", out).Output()
	if err != nil {
		t.Fatalf("unsquashfs -l: %v", err)
	}
	for _, want := range []string{"squashfs-root/cnt/samtools/1.21/bin/samtools", "squashfs-root" + meta.Path} {
		if !strings.Contains(string(list), want) {
			t.Errorf("packed image is missing %s:\n%s", want, list)
		}
	}
}

// The build's own scratch directory sits beside the payload and the manifest, so
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
	script, _, _ := buildSquashfsOpts(b, false, b.ws.CntDir, metaDir, out)
	if output, err := exec.CommandContext(t.Context(), "/bin/bash", "-c", script).CombinedOutput(); err != nil {
		t.Fatalf("mksquashfs failed: %v\n%s", err, output)
	}

	list, err := exec.CommandContext(t.Context(), "unsquashfs", "-l", out).Output()
	if err != nil {
		t.Fatal(err)
	}
	if strings.Contains(string(list), "huge.pkg") {
		t.Errorf("build scratch leaked into the image:\n%s", list)
	}
}

// The conda install phase must not know where the image lands: that is what lets
// a cancelled install leave nothing next to the installed images.
func TestCondaInstallDoesNotPackOrTouchTarget(t *testing.T) {
	withAppTmpOverlay(t, false)
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
