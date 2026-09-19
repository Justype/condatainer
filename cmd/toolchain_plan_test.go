package cmd

import (
	"errors"
	"slices"
	"strings"
	"testing"
)

func TestPlanToolchain(t *testing.T) {
	none := errors.New("not found")
	cases := []struct {
		name  string
		probe systemProbe
		want  []string
		note  string // a substring the report must contain
	}{
		{"a complete host needs nothing", systemProbe{Apptainer: "1.5.3", Squashfs: true, Squashfuse: true, Fuse2fs: true}, nil, "1.5.3 ok"},
		{"no usable apptainer installs apptainer alone", systemProbe{ApptainerErr: none, Fuse2fs: true}, []string{"apptainer"}, "install apptainer"},
		{"apptainer already brings the squashfs tools", systemProbe{ApptainerErr: none}, []string{"apptainer"}, "fuse2fs: not found"},
		{"a missing squashfs-tools with a good apptainer", systemProbe{Apptainer: "1.5.3", Squashfuse: true, Fuse2fs: true}, []string{"squashfs-tools"}, "install squashfs-tools"},
		{"a missing squashfuse", systemProbe{Apptainer: "1.5.3", Squashfs: true, Fuse2fs: true}, []string{"squashfuse"}, "install squashfuse"},
		{"both missing", systemProbe{Apptainer: "1.5.3", Fuse2fs: true}, []string{"squashfs-tools", "squashfuse"}, ""},
		{"fuse2fs is only reported", systemProbe{Apptainer: "1.5.3", Squashfs: true, Squashfuse: true}, nil, "fuse2fs: not found"},
	}
	for _, c := range cases {
		install, report, _ := planToolchain(c.probe)
		if !slices.Equal(install, c.want) {
			t.Errorf("%s: install = %v, want %v", c.name, install, c.want)
		}
		if !strings.Contains(strings.Join(report, "\n"), c.note) {
			t.Errorf("%s: report %q lacks %q", c.name, report, c.note)
		}
	}
}

func TestPlanToolchainWarnsAboutAMissingSystemApptainer(t *testing.T) {
	none := errors.New("not found")
	ok := systemProbe{Apptainer: "installed in libexec", Squashfs: true, Squashfuse: true, Fuse2fs: true}

	// libexec's apptainer serves exec/run, but os overlays need the system one.
	missing := ok
	missing.SystemApptainerErr = none
	if _, _, warnings := planToolchain(missing); len(warnings) != 1 || !strings.Contains(warnings[0], "os overlays") {
		t.Errorf("warnings = %q, want one naming os overlays", warnings)
	}

	if _, _, warnings := planToolchain(ok); len(warnings) != 0 {
		t.Errorf("warned with a system apptainer present: %q", warnings)
	}
}

func TestLibexecHostError(t *testing.T) {
	if err := libexecHostError(true); !errors.Is(err, errLibexecInContainer) {
		t.Errorf("libexecHostError(true) = %v, want errLibexecInContainer", err)
	}
	if err := libexecHostError(false); err != nil {
		t.Errorf("libexecHostError(false) = %v, want nil", err)
	}
}
