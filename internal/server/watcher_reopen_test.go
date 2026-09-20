package server

import (
	"testing"
	"time"
)

func TestReopenInterval(t *testing.T) {
	for failures, want := range map[int]time.Duration{
		0: 30 * time.Second, 1: 30 * time.Second, 2: time.Minute,
		3: 2 * time.Minute, 4: 4 * time.Minute, 5: 5 * time.Minute, 50: 5 * time.Minute,
	} {
		if got := reopenInterval(failures); got != want {
			t.Errorf("reopenInterval(%d) = %v, want %v", failures, got, want)
		}
	}
}
