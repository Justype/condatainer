package scheduler

import (
	"testing"
	"time"
)

func TestParseMemoryMB(t *testing.T) {
	tests := []struct {
		input  string
		wantMB int64
	}{
		{"8G", 8 * 1024},
		{"8GB", 8 * 1024},
		{"1024M", 1024},
		{"1024MB", 1024},
		{"4096K", 4},
		{"4096KB", 4},
		{"1T", 1024 * 1024},
		{"1TB", 1024 * 1024},
	}

	for _, tt := range tests {
		t.Run(tt.input, func(t *testing.T) {
			mb, err := parseMemoryMB(tt.input)
			if err != nil {
				t.Errorf("parseMemoryMB(%q) error: %v", tt.input, err)
				return
			}
			if mb != tt.wantMB {
				t.Errorf("parseMemoryMB(%q) = %d MB; want %d MB", tt.input, mb, tt.wantMB)
			}
		})
	}
}

func TestParseHMSTime(t *testing.T) {
	tests := []struct {
		name    string
		input   string
		wantDur time.Duration
		wantErr bool
	}{
		{"HH:MM:SS", "02:30:00", 2*time.Hour + 30*time.Minute, false},
		{"HH:MM", "10:30", 10*time.Hour + 30*time.Minute, false},
		{"minutes only", "90", 90 * time.Minute, false},
		{"with seconds", "01:00:30", time.Hour + 30*time.Second, false},
		{"empty string", "", 0, false},
		{"zero", "00:00", 0, false},
		{"large via HH:MM:SS", "168:00:00", 168 * time.Hour, false},
		{"large via HH:MM", "168:00", 168 * time.Hour, false},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			dur, err := parseHMSTime(tt.input)
			if tt.wantErr && err == nil {
				t.Error("Expected error, got nil")
			}
			if !tt.wantErr && err != nil {
				t.Errorf("Unexpected error: %v", err)
			}
			if dur != tt.wantDur {
				t.Errorf("parseHMSTime(%q) = %v; want %v", tt.input, dur, tt.wantDur)
			}
		})
	}
}

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
