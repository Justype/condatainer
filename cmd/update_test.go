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
