package cmd

import (
	"strings"
	"testing"

	"github.com/fatih/color"

	"github.com/Justype/condatainer/internal/artifact/meta"
	"github.com/Justype/condatainer/internal/store"
)

func candidate(name, digest string, layout store.Layout) store.Candidate {
	return store.Candidate{
		Name:     name,
		Layout:   layout,
		Identity: meta.KeyRef{Scheme: "conda-explicit-v1", SHA256: digest},
	}
}

// A flat artifact answers to its plain name, so its identity is informational
// and greyed; a store entry can only be addressed by one, so it stays legible.
func TestStyleStoreAddressGreysOnlyAFlatIdentity(t *testing.T) {
	previous := color.NoColor
	color.NoColor = false
	defer func() { color.NoColor = previous }()

	digest := strings.Repeat("a", 64)
	flat := styleStoreAddress(candidate("cutadapt/5.0", digest, store.LayoutFlat))
	stored := styleStoreAddress(candidate("cutadapt/5.0", digest, store.LayoutStored))

	const grey = "\x1b[90m"
	if !strings.Contains(flat, grey+"@") {
		t.Errorf("a flat identity is not greyed: %q", flat)
	}
	if strings.Contains(stored, grey+"@") {
		t.Errorf("a store identity was greyed: %q", stored)
	}
	// Both remain pasteable: the address survives whatever styling is applied.
	for _, rendered := range []string{flat, stored} {
		if !strings.Contains(stripANSI(rendered), "cutadapt/5.0@"+digest[:store.DefaultPrefixChars]) {
			t.Errorf("address lost in styling: %q", rendered)
		}
	}
}

func stripANSI(s string) string {
	var out strings.Builder
	for i := 0; i < len(s); i++ {
		if s[i] == '\x1b' {
			for i < len(s) && s[i] != 'm' {
				i++
			}
			continue
		}
		out.WriteByte(s[i])
	}
	return out.String()
}

// Splitting on the first @ is what lets a full scheme-backed key, which carries
// its own @, be given inline.
func TestSplitStoreAddress(t *testing.T) {
	digest := strings.Repeat("b", 64)
	for _, tc := range []struct {
		name, arg, flag, wantName string
		wantErr                   bool
	}{
		{name: "inline prefix", arg: "star/2.7.11b@bbbbbbbbbbbb", wantName: "star/2.7.11b"},
		{name: "flag only", arg: "star/2.7.11b", flag: "bbbbbbbbbbbb", wantName: "star/2.7.11b"},
		{name: "complete key inline", arg: "star/2.7.11b@identity-v1@sha256:" + digest, wantName: "star/2.7.11b"},
		{name: "agreeing duplicate", arg: "star/2.7.11b@bbbbbbbbbbbb", flag: "bbbbbbbbbbbb", wantName: "star/2.7.11b"},
		{name: "conflict", arg: "star/2.7.11b@bbbbbbbbbbbb", flag: "cccccccccccc", wantErr: true},
		{name: "neither", arg: "star/2.7.11b", wantErr: true},
		{name: "empty inline", arg: "star/2.7.11b@", wantErr: true},
	} {
		t.Run(tc.name, func(t *testing.T) {
			name, query, err := splitStoreAddress(tc.arg, tc.flag)
			if tc.wantErr {
				if err == nil {
					t.Fatalf("accepted %q with --identity %q", tc.arg, tc.flag)
				}
				return
			}
			if err != nil {
				t.Fatal(err)
			}
			if name != tc.wantName {
				t.Errorf("name = %q, want %q", name, tc.wantName)
			}
			if query.SHA256 == "" {
				t.Error("no digest parsed")
			}
		})
	}
}
