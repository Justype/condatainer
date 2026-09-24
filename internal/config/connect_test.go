package config

import "testing"

func TestParseConnect(t *testing.T) {
	for in, want := range map[string]string{"auto": "auto", " SSH ": "ssh", "Scheduler": "scheduler", "direct": "direct"} {
		if got, ok := ParseConnect(in); !ok || got != want {
			t.Errorf("ParseConnect(%q) = %q, %v; want %q", in, got, ok, want)
		}
	}
	for _, in := range []string{"", "true", "bind_all", "srun"} {
		if _, ok := ParseConnect(in); ok {
			t.Errorf("ParseConnect(%q) accepted", in)
		}
	}
}
