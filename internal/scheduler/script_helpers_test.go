package scheduler

import (
	"testing"
	"time"
)

func TestFormatHMSTime(t *testing.T) {
	tests := []struct {
		input time.Duration
		want  string
	}{
		{2*time.Hour + 30*time.Minute, "02:30:00"},
		{time.Hour + 30*time.Second, "01:00:30"},
		{90 * time.Minute, "01:30:00"},
		{168 * time.Hour, "168:00:00"},
		{0, "00:00:00"},
	}

	for _, tt := range tests {
		t.Run(tt.want, func(t *testing.T) {
			got := formatHMSTime(tt.input)
			if got != tt.want {
				t.Errorf("formatHMSTime(%v) = %q; want %q", tt.input, got, tt.want)
			}
		})
	}
}

