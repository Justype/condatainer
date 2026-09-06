package conda

import (
	"os"
	"os/exec"
	"path/filepath"
	"testing"
)

func requireSquashfsTools(t *testing.T) {
	t.Helper()
	for _, bin := range []string{"mksquashfs", "unsquashfs"} {
		if _, err := exec.LookPath(bin); err != nil {
			t.Skipf("%s not available", bin)
		}
	}
}

// packSqf builds a .sqf whose root is dir.
func packSqf(t *testing.T, dir, out string) {
	t.Helper()
	cmd := exec.Command("mksquashfs", dir, out, "-no-progress", "-noappend", "-quiet", "-no-xattrs")
	if output, err := cmd.CombinedOutput(); err != nil {
		t.Fatalf("mksquashfs: %v\n%s", err, output)
	}
}

func writeFile(t *testing.T, path, content string) {
	t.Helper()
	if err := os.MkdirAll(filepath.Dir(path), 0o755); err != nil {
		t.Fatal(err)
	}
	if err := os.WriteFile(path, []byte(content), 0o644); err != nil {
		t.Fatal(err)
	}
}

func TestListCondaPackagesSqfReadsCondaMeta(t *testing.T) {
	requireSquashfsTools(t)
	root := t.TempDir()
	writeFile(t, filepath.Join(root, "cnt_env/conda-meta/numpy-1.24.0-py311h1234567_0.json"), "{}")
	writeFile(t, filepath.Join(root, "cnt_env/conda-meta/history"), "")
	sqf := filepath.Join(t.TempDir(), "env.sqf")
	packSqf(t, root, sqf)

	pkgs, err := ListCondaPackagesSqf(sqf)
	if err != nil {
		t.Fatalf("ListCondaPackagesSqf: %v", err)
	}
	if pkgs["numpy"] != "1.24.0" {
		t.Fatalf("pkgs = %v, want numpy=1.24.0", pkgs)
	}
}

func TestListCondaPackagesSqfNonSqfPath(t *testing.T) {
	pkgs, err := ListCondaPackagesSqf("/tmp/does-not-matter.img")
	if err != nil || pkgs != nil {
		t.Fatalf("ListCondaPackagesSqf(non-sqf) = (%v, %v), want (nil, nil)", pkgs, err)
	}
}

// The .img's own history is the continuation of the snapshot's, so merging
// reads them as one log: channels union in order, and the .img's own later
// spec line wins for a name reinstalled after the freeze.
func TestReadCondaInfoMergedConcatenatesAsOneLog(t *testing.T) {
	requireSquashfsTools(t)

	snapRoot := t.TempDir()
	writeFile(t, filepath.Join(snapRoot, "cnt_env/conda-meta/history"),
		"+https://conda.anaconda.org/conda-forge/linux-64::numpy-1.24.0-py311h1234567_0\n"+
			"# update specs: [\"numpy\"]\n")
	snapSqf := filepath.Join(t.TempDir(), "env.sqf")
	packSqf(t, snapRoot, snapSqf)

	imgRoot := t.TempDir()
	writeFile(t, filepath.Join(imgRoot, "cnt_env/conda-meta/history"),
		"+https://conda.anaconda.org/bioconda/linux-64::numpy-1.26.0-py311h1234567_0\n"+
			"# update specs: [\"numpy=1.26\"]\n")
	// A second .sqf standing in for the .img: image.ReadFile only branches on
	// extension, and this test exercises the merge/parse logic, not the ext3
	// read path (already covered by ListCondaPackages's own tests).
	imgSqf := filepath.Join(t.TempDir(), "img-as.sqf")
	packSqf(t, imgRoot, imgSqf)

	info := ReadCondaInfoMerged(snapSqf, imgSqf, "/cnt_env")
	if info == nil {
		t.Fatal("ReadCondaInfoMerged returned nil")
	}
	if len(info.Channels) != 2 || info.Channels[0] != "conda-forge" || info.Channels[1] != "bioconda" {
		t.Fatalf("channels = %v, want [conda-forge bioconda]", info.Channels)
	}
	if len(info.Specs) != 1 || info.Specs[0] != "numpy=1.26" {
		t.Fatalf("specs = %v, want [numpy=1.26] (the .img's own later spec wins)", info.Specs)
	}
}

func TestReadCondaInfoMergedNilWhenNeitherHasHistory(t *testing.T) {
	requireSquashfsTools(t)
	empty := filepath.Join(t.TempDir(), "empty.sqf")
	packSqf(t, t.TempDir(), empty)

	if info := ReadCondaInfoMerged(empty, empty, "/cnt_env"); info != nil {
		t.Fatalf("info = %v, want nil", info)
	}
}
