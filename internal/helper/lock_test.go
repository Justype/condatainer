package helper

import (
	"errors"
	"testing"

	"github.com/Justype/condatainer/internal/image/producer"
)

func TestHeadlessLiveness(t *testing.T) {
	t.Setenv("CNT_HELPER_ID", "t-1")
	t.Setenv("CNT_HELPER_STATE_DIR", t.TempDir())

	if live, _ := HeadlessLiveness("t-1"); live != LivenessUnknown {
		t.Fatalf("no lock file: got %v, want unknown", live)
	}

	info := producer.Info{Runner: "local", Node: producer.ShortHostname(), PID: 4242}
	lock, err := HoldLock(LockFilePath("t-1"), info)
	if err != nil {
		t.Fatal(err)
	}
	live, got := HeadlessLiveness("t-1")
	if live != LivenessAlive || got.PID != 4242 {
		t.Fatalf("held: got %v %+v, want alive with the recorded owner", live, got)
	}

	lock.Close()
	if live, got = HeadlessLiveness("t-1"); live != LivenessGone || got.PID != 4242 {
		t.Fatalf("released: got %v %+v, want gone with the recorded owner", live, got)
	}
}

func TestKillHeadlessProcessRefusesOtherHost(t *testing.T) {
	t.Setenv("CNT_HELPER_ID", "t-2")
	t.Setenv("CNT_HELPER_STATE_DIR", t.TempDir())

	lock, err := HoldLock(LockFilePath("t-2"), producer.Info{Runner: "local", Node: "elsewhere", PID: 4242})
	if err != nil {
		t.Fatal(err)
	}
	defer lock.Close()
	if err := KillHeadlessProcess("t-2"); !errors.Is(err, ErrOtherHost) {
		t.Fatalf("got %v, want ErrOtherHost", err)
	}
}
