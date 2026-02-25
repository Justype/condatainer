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
