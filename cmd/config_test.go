package cmd

import (
    "testing"

    "github.com/Justype/condatainer/internal/config"
)

func TestConfigValueCompletion(t *testing.T) {
    opts := configValueCompletion("build.compress_args")
    expected := config.CompressNames()
    for _, e := range expected {
        found := false
        for _, o := range opts {
            if o == e {
                found = true
                break
            }
        }
        if !found {
            t.Errorf("expected completion option %q not present", e)
        }
    }
}

func TestConfigSetNormalizeCompressArgs(t *testing.T) {
    // exercise the normalization path that would be executed when the
    // configSetCmd.Run logic runs
    tests := map[string]string{
        "gzip":       "-comp gzip",
        "lz4":        "-comp lz4",
        "zstd":       "-comp zstd -Xcompression-level 14",
        "zstd-fast":  "-comp zstd -Xcompression-level 3",
        "zstd-medium": "-comp zstd -Xcompression-level 8",
        "zstd-high":  "-comp zstd -Xcompression-level 19",
        // these should be left untouched
        "-comp foo": "-comp foo",
    }

    for input, want := range tests {
        got := config.NormalizeCompressArgs(input)
        if got != want {
            t.Errorf("NormalizeCompressArgs(%q) = %q, want %q", input, got, want)
        }
    }
}
