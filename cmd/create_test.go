package cmd

import (
	"strings"
	"testing"

	"github.com/Justype/condatainer/catalog"
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

// --name picks the target; --file must not quietly replace it with a prefix
// derived from the filename. The mode dispatch tests --prefix first, so a
// derived prefix used to make the --name branch unreachable.
func TestDerivePrefixFromFile(t *testing.T) {
	cases := []struct {
		name             string
		file, prefix, nm string
		want             string
	}{
		{"file alone derives a prefix", "environment.yml", "", "", "environment"},
		{"name wins over the filename", "environment.yml", "", "myenv", ""},
		{"explicit prefix is kept", "environment.yml", "/images/x", "", ""},
		{"no file, nothing to derive", "", "", "", ""},
		{"path keeps its directory", "/tmp/envs/build.sh", "", "", "/tmp/envs/build"},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			if got := derivePrefixFromFile(tc.file, tc.prefix, tc.nm); got != tc.want {
				t.Errorf("derivePrefixFromFile(%q, %q, %q) = %q, want %q",
					tc.file, tc.prefix, tc.nm, got, tc.want)
			}
		})
	}
}

// Restoring a project from a lockfile recreates the names the catalog uses, and
// a data image is several levels deep. The name has to survive the round trip
// through the filename, which is what the depth limit used to be guarding.
func TestNormalizedTargetNameKeepsDepth(t *testing.T) {
	prev := createName
	t.Cleanup(func() { createName = prev })

	for _, want := range []string{
		"myenv",
		"samtools/1.23.1",
		"grch38/star/2.7.11b/gencode47-101",
	} {
		createName = want
		got := normalizedTargetName()
		if got != want {
			t.Errorf("normalizedTargetName() = %q, want %q", got, want)
		}
		// / becomes -- on the way to a filename, and back on the way in.
		roundTrip := catalog.Normalize(strings.ReplaceAll(got, "/", "--"))
		if roundTrip != want {
			t.Errorf("round trip through filename = %q, want %q", roundTrip, want)
		}
	}
}
