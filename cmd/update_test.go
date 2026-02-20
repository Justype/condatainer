package cmd

import "testing"

func TestCompareVersions(t *testing.T) {
	cases := []struct {
		v1, v2 string
		want   int
	}{
		{"v1.2.3", "v1.2.3", 0},
		{"1.2.3", "v1.2.3", 0},
		{"v1.2.3-alpha", "v1.2.3", -1},
		{"v1.2.3", "v1.2.3-alpha", 1},
		{"v1.2.3-alpha", "v1.2.3-beta", -1},
		{"v1.2.3+001", "v1.2.3+002", -1},
		{"v1.2.4", "v1.2.3+999", 1},
	}
	for _, c := range cases {
		got := compareVersions(c.v1, c.v2)
		if got != c.want {
			t.Errorf("compareVersions(%q,%q) = %d; want %d", c.v1, c.v2, got, c.want)
		}
	}
}

func TestChangeHelpers(t *testing.T) {
	// major change detection
	if !isMajorChange("v1.2.3", "v2.0.0") {
		t.Error("expected major change from 1.x to 2.x")
	}
	if isMajorChange("v1.2.3", "v1.3.0") {
		t.Error("did not expect major change when major versions equal")
	}

	// base image update (major/minor) logic
	if !isMinorChange("v1.2.3", "v1.3.0") {
		t.Error("expected base image update when minor version increases")
	}
	if isMinorChange("v1.2.3", "v1.2.5") {
		t.Error("should not update base image for patch-only change")
	}
	if !isMinorChange("v1.2.3", "v2.0.0") {
		t.Error("expected base image update when major version increases")
	}

	// helper extraction
	if maj, err := getMajorNumber("v3.4.5"); err != nil || maj != 3 {
		t.Errorf("getMajorNumber returned %d,%v; want 3,nil", maj, err)
	}
	if mm := getMajorMinor("1.2.3"); mm != "v1.2" {
		t.Errorf("getMajorMinor returned %q; want \"v1.2\"", mm)
	}
}
