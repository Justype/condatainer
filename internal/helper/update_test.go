package helper

import (
	"errors"
	"os"
	"testing"
)

func TestHelperIndexNotFound(t *testing.T) {
	tests := []struct {
		name string
		err  error
		want bool
	}{
		{name: "HTTP 404", err: errors.New("catalog: https://example.test/index/helpers.json: 404 Not Found"), want: true},
		{name: "local absence", err: &os.PathError{Op: "open", Path: "index/helpers.json", Err: os.ErrNotExist}, want: true},
		{name: "HTTP 403", err: errors.New("catalog: https://example.test/index/helpers.json: 403 Forbidden")},
		{name: "server error", err: errors.New("catalog: https://example.test/index/helpers.json: 500 Internal Server Error")},
		{name: "malformed index", err: errors.New("failed to parse helper index: invalid character")},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			if got := helperIndexNotFound(tt.err); got != tt.want {
				t.Errorf("helperIndexNotFound(%q) = %t, want %t", tt.err, got, tt.want)
			}
		})
	}
}
