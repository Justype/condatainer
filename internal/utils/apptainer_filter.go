package utils

import (
	"bytes"
	"io"
	"strings"
)

// Apptainer log-line prefixes, for selecting which lines to drop.
const (
	ApptainerWarning = "WARNING:"
	ApptainerInfo    = "INFO:"
	ApptainerDebug   = "DEBUG:"
	ApptainerVerbose = "VERBOSE:"
	ApptainerError   = "ERROR:"
	ApptainerFatal   = "FATAL:"
)

// ApptainerNonError lists the non-error log prefixes; dropping these keeps only
// ERROR:/FATAL: and non-Apptainer output.
var ApptainerNonError = []string{ApptainerWarning, ApptainerInfo, ApptainerDebug, ApptainerVerbose}

// ApptainerFilter is an io.Writer that drops Apptainer log lines whose prefix is
// in drop (e.g. noisy bind-mount warnings) and forwards the rest. Partial lines
// are buffered until a newline; call Flush to emit a trailing unterminated line.
type ApptainerFilter struct {
	w    io.Writer
	buf  bytes.Buffer
	drop []string
}

// NewApptainerFilter wraps w, dropping lines that start with any of dropPrefixes.
func NewApptainerFilter(w io.Writer, dropPrefixes ...string) *ApptainerFilter {
	return &ApptainerFilter{w: w, drop: dropPrefixes}
}

// Write buffers p and flushes every complete line, dropping selected lines.
func (f *ApptainerFilter) Write(p []byte) (int, error) {
	f.buf.Write(p)
	for {
		i := bytes.IndexByte(f.buf.Bytes(), '\n')
		if i < 0 {
			break
		}
		line := f.buf.Next(i + 1) // includes the newline
		if err := f.emit(line); err != nil {
			return len(p), err
		}
	}
	return len(p), nil
}

// Flush emits any buffered line that has no trailing newline.
func (f *ApptainerFilter) Flush() error {
	if f.buf.Len() == 0 {
		return nil
	}
	line := f.buf.Next(f.buf.Len())
	return f.emit(line)
}

// emit writes a line unless its prefix is in drop.
func (f *ApptainerFilter) emit(line []byte) error {
	trimmed := strings.TrimSpace(string(line))
	for _, prefix := range f.drop {
		if strings.HasPrefix(trimmed, prefix) {
			return nil
		}
	}
	_, err := f.w.Write(line)
	return err
}
