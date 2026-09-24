package proxy

import (
	"bufio"
	"bytes"
	"fmt"
	"io"
	"strconv"
	"strings"
)

// execMaxPayload bounds one data frame; senders stay at execChunk.
const execMaxPayload = 1 << 20

// frame is one message between the dashboard and the relay. Control frames are
// one text line "op id [arg]"; a data frame ("D") is the line "D id len"
// followed by len raw bytes and a newline. The trailing newline is what makes
// srun deliver the frame instead of holding it back.
type frame struct {
	op, id, arg string
	data        []byte
}

// writeFrame writes f with a single Write, so frames from one writer never interleave.
func writeFrame(w io.Writer, f frame) error {
	var b bytes.Buffer
	b.WriteString(f.op)
	if f.id != "" {
		b.WriteString(" " + f.id)
	}
	if f.op == "D" {
		b.WriteString(" " + strconv.Itoa(len(f.data)) + "\n")
		b.Write(f.data)
	} else if f.arg != "" {
		b.WriteString(" " + f.arg)
	}
	b.WriteByte('\n')
	_, err := w.Write(b.Bytes())
	return err
}

// readFrame returns the next frame. A line that is not a frame (a warning
// printed on the same stream) comes back with an op no caller handles; a
// malformed "D" line is dropped. It returns an error only when the stream fails.
func readFrame(r *bufio.Reader) (frame, error) {
	for {
		line, err := r.ReadString('\n')
		if err != nil {
			return frame{}, err
		}
		p := strings.Fields(line)
		if len(p) == 0 {
			continue
		}
		f := frame{op: p[0]}
		if len(p) > 1 {
			f.id = p[1]
		}
		if len(p) > 2 {
			f.arg = p[2]
		}
		if f.op != "D" {
			return f, nil
		}
		n, err := strconv.Atoi(f.arg)
		if err != nil || n < 0 || n > execMaxPayload || f.id == "" {
			continue
		}
		f.data = make([]byte, n)
		if _, err := io.ReadFull(r, f.data); err != nil {
			return frame{}, err
		}
		if c, err := r.ReadByte(); err != nil {
			return frame{}, err
		} else if c != '\n' {
			return frame{}, fmt.Errorf("exec frame: missing terminator after %d bytes", n)
		}
		return f, nil
	}
}
