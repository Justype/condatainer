package config

import "testing"

func TestNormalizeCompressArgs(t *testing.T) {
	cases := map[string]string{
		"gzip":        "-comp gzip",
		"lz4":         "-comp lz4",
		"zstd":        "-comp zstd -Xcompression-level 14",
		"zstd-fast":   "-comp zstd -Xcompression-level 3",
		"zstd-medium": "-comp zstd -Xcompression-level 8",
		"zstd-high":   "-comp zstd -Xcompression-level 19",
		// unknown values should be returned unchanged
		"-comp foo": "-comp foo",
		"":          "",
	}

	for input, want := range cases {
		got := NormalizeCompressArgs(input)
		if got != want {
			t.Errorf("NormalizeCompressArgs(%q) = %q; want %q", input, got, want)
		}
	}
}

func TestIsValidDistro(t *testing.T) {
	for _, d := range GetAvailableDistros() {
		if !IsValidDistro(d) {
			t.Errorf("expected distro %q to be valid", d)
		}
	}
	if IsValidDistro("not-a-distro") {
		t.Errorf("unexpectedly accepted invalid distro")
	}
}

func TestArgsForCompressAndNames(t *testing.T) {
	names := CompressNames()
	if len(names) == 0 {
		t.Fatalf("no compress names returned")
	}

	for _, name := range names {
		args := ArgsForCompress(name)
		if args == name {
			t.Errorf("ArgsForCompress(%q) returned input unchanged", name)
		}
	}
	// unknown name should be returned unchanged
	if got := ArgsForCompress("foo"); got != "foo" {
		t.Errorf("ArgsForCompress of unknown name changed value: %q", got)
	}
}

func TestParseLatestModuleFromAvailOutput(t *testing.T) {
	output := `
/opt/modulefiles:
apptainer/1.3.5
apptainer/1.4.0(default)
apptainer/1.2.9
`

	module, version := parseLatestModuleFromAvailOutput("apptainer", output)
	if module != "apptainer/1.4.0" {
		t.Fatalf("unexpected module: got %q", module)
	}
	if version != "1.4.0" {
		t.Fatalf("unexpected version: got %q", version)
	}
}

func TestParseLatestModuleFromAvailOutputWithMixedNoise(t *testing.T) {
	output := `
Lmod has detected the following error: ...
----------------------------------- /opt/modulefiles -----------------------------------
singularity/3.10
singularity/3.11*
otherpkg/1.0
`

	module, version := parseLatestModuleFromAvailOutput("singularity", output)
	if module != "singularity/3.11" {
		t.Fatalf("unexpected module: got %q", module)
	}
	if version != "3.11" {
		t.Fatalf("unexpected version: got %q", version)
	}
}

func TestCompareModuleVersion(t *testing.T) {
	tests := []struct {
		a    string
		b    string
		want int
	}{
		{a: "1.4.0", b: "1.3.9", want: 1},
		{a: "1.4", b: "1.4.0", want: -1},
		{a: "3.11", b: "3.11", want: 0},
		{a: "", b: "3.11", want: -1},
		{a: "3.11", b: "", want: 1},
		{a: "2024a", b: "2023b", want: 1},
	}

	for _, tt := range tests {
		got := compareModuleVersion(tt.a, tt.b)
		switch {
		case tt.want > 0 && got <= 0:
			t.Fatalf("compareModuleVersion(%q, %q) = %d; want > 0", tt.a, tt.b, got)
		case tt.want < 0 && got >= 0:
			t.Fatalf("compareModuleVersion(%q, %q) = %d; want < 0", tt.a, tt.b, got)
		case tt.want == 0 && got != 0:
			t.Fatalf("compareModuleVersion(%q, %q) = %d; want 0", tt.a, tt.b, got)
		}
	}
}
