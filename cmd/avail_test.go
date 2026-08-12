package cmd

import (
	"strings"
	"testing"
)

func TestFormatPackageLinePutsDescriptionOnSecondLine(t *testing.T) {
	pkg := PackageInfo{
		Name:        "ubuntu24/simple-os",
		IsContainer: true,
		Source:      "test",
		Description: "Small OS overlay providing jq",
	}

	got := formatPackageLine(pkg, true)
	if !strings.Contains(got, "\n  "+pkg.Description) {
		t.Fatalf("formatPackageLine = %q, want an indented description line", got)
	}
	if strings.Contains(strings.SplitN(got, "\n", 2)[0], pkg.Description) {
		t.Fatalf("description remained on the metadata line: %q", got)
	}
}

func TestFormatPackageLineCanHideDescription(t *testing.T) {
	pkg := PackageInfo{Name: "hello/1.0", Description: "distinct description"}
	if got := formatPackageLine(pkg, false); strings.Contains(got, pkg.Description) {
		t.Fatalf("formatPackageLine = %q, want no description", got)
	}
}

func TestFormatDescriptionWrapsAndIndents(t *testing.T) {
	got := formatDescription("one two three four", 2, 12)
	want := "  one two\n  three four"
	if got != want {
		t.Fatalf("formatDescription = %q, want %q", got, want)
	}
}

func TestFormatListDescriptionUsesSecondLine(t *testing.T) {
	got := formatListDescription("hello/1.0", "one two three four", 13)
	want := "hello/1.0\n  one two\n  three four"
	if got != want {
		t.Fatalf("formatListDescription = %q, want %q", got, want)
	}
}
