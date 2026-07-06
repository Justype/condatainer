package cmd

import "testing"

func TestParsePostScriptHelperFlags(t *testing.T) {
	flags, rest, err := parsePostScriptHelperFlags([]string{
		"-m", "16g",
		"--cpus=6",
		"-t12h",
		"--gpu", "a100:1",
		"-b", "base.sif",
		"-eenv.img",
		"-o", "one.sqfs",
		"--overlay=two.sqfs",
		"-w", "work/dir",
		"--new",
		"-p", "11908",
	})
	if err != nil {
		t.Fatalf("parsePostScriptHelperFlags returned error: %v", err)
	}
	if !flags.memSet || flags.mem != "16g" {
		t.Fatalf("mem flag = (%v, %q), want (true, 16g)", flags.memSet, flags.mem)
	}
	if !flags.cpusSet || flags.cpus != 6 {
		t.Fatalf("cpus flag = (%v, %d), want (true, 6)", flags.cpusSet, flags.cpus)
	}
	if !flags.timeSet || flags.time != "12h" {
		t.Fatalf("time flag = (%v, %q), want (true, 12h)", flags.timeSet, flags.time)
	}
	if !flags.gpuSet || flags.gpu != "a100:1" {
		t.Fatalf("gpu flag = (%v, %q), want (true, a100:1)", flags.gpuSet, flags.gpu)
	}
	if !flags.baseSet || flags.base != "base.sif" {
		t.Fatalf("base flag = (%v, %q), want (true, base.sif)", flags.baseSet, flags.base)
	}
	if !flags.envSet || flags.env != "env.img" {
		t.Fatalf("env flag = (%v, %q), want (true, env.img)", flags.envSet, flags.env)
	}
	if !flags.cwdSet || flags.cwd != "work/dir" {
		t.Fatalf("cwd flag = (%v, %q), want (true, work/dir)", flags.cwdSet, flags.cwd)
	}
	if !flags.newSet {
		t.Fatalf("new flag = %v, want true", flags.newSet)
	}
	if len(flags.overlays) != 2 || flags.overlays[0] != "one.sqfs" || flags.overlays[1] != "two.sqfs" {
		t.Fatalf("overlays = %#v, want [one.sqfs two.sqfs]", flags.overlays)
	}
	if len(rest) != 2 || rest[0] != "-p" || rest[1] != "11908" {
		t.Fatalf("rest = %#v, want [-p 11908]", rest)
	}
}

func TestParsePostScriptHelperFlagsPassesThroughAfterDoubleDash(t *testing.T) {
	_, rest, err := parsePostScriptHelperFlags([]string{"-m", "16g", "--", "--mem", "helper-value"})
	if err != nil {
		t.Fatalf("parsePostScriptHelperFlags returned error: %v", err)
	}
	if len(rest) != 3 || rest[0] != "--" || rest[1] != "--mem" || rest[2] != "helper-value" {
		t.Fatalf("rest = %#v, want [-- --mem helper-value]", rest)
	}
}

func TestParsePostScriptHelperFlagsRequiresValues(t *testing.T) {
	if _, _, err := parsePostScriptHelperFlags([]string{"-m"}); err == nil {
		t.Fatal("expected missing -m value to return an error")
	}
	if _, _, err := parsePostScriptHelperFlags([]string{"--cpus=0"}); err == nil {
		t.Fatal("expected invalid --cpus value to return an error")
	}
}
