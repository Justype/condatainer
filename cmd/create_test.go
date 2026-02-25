package cmd

import (
    "testing"

    "github.com/Justype/condatainer/internal/config"
)

// ensure that every compress option declared in config is registered as a
// flag on the create command.  This guards against drift when new options are
// added.
func TestCreateFlagsForCompressOptions(t *testing.T) {
    for _, opt := range config.CompressOptions {
        if createCmd.Flags().Lookup(opt.Name) == nil {
            t.Errorf("create command missing flag for compression option %q", opt.Name)
        }
    }
}

func TestCompressArgsFromFlags(t *testing.T) {
    // helper to build a map with all flags set to the given boolean value
    makeMap := func(setName string) map[string]*bool {
        m := make(map[string]*bool)
        for _, opt := range config.CompressOptions {
            v := false
            if opt.Name == setName {
                v = true
            }
            m[opt.Name] = &v
        }
        return m
    }

    // no flag => empty string, no error
    if got, err := compressArgsFromFlags(makeMap("")); err != nil {
        t.Fatalf("unexpected error for no flag: %v", err)
    } else if got != "" {
        t.Errorf("expected empty result for no flag, got %q", got)
    }

    // each individual option returns the appropriate args
    for _, opt := range config.CompressOptions {
        m := makeMap(opt.Name)
        if got, err := compressArgsFromFlags(m); err != nil {
            t.Errorf("unexpected error for option %q: %v", opt.Name, err)
        } else if got != opt.Args {
            t.Errorf("compressArgsFromFlags(%q) = %q, want %q", opt.Name, got, opt.Args)
        }
    }

    // multiple options should error
    m := makeMap("")
    if len(config.CompressOptions) >= 2 {
        // set first two
        names := []string{config.CompressOptions[0].Name, config.CompressOptions[1].Name}
        for _, n := range names {
            v := true
            m[n] = &v
        }
        if _, err := compressArgsFromFlags(m); err == nil {
            t.Errorf("expected error when multiple compression flags set")
        }
    }
}
