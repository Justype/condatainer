package proxy

import (
	"bufio"
	"context"
	"errors"
	"fmt"
	"io"
	"net"
	"os/exec"
	"strconv"
	"strings"
	"sync"
	"sync/atomic"
	"time"
)

const (
	execChunk        = 16384
	execQueueLen     = 256
	execReadyTimeout = 30 * time.Second
	execStderrKeep   = 2048
)

var errExecEnded = errors.New("exec tunnel ended")

// DialViaExec runs cmd, which must end up running RunExecRelay on the node of a
// running job, and returns a DialFunc that reaches loopback ports there through
// it. done closes when the command exits.
func DialViaExec(cmd []string) (DialFunc, func(), <-chan struct{}, error) {
	if len(cmd) == 0 {
		return nil, nil, nil, errors.New("exec tunnel: empty command")
	}
	c := exec.Command(cmd[0], cmd[1:]...)
	stderr := &tailBuffer{}
	c.Stderr = stderr
	in, err := c.StdinPipe()
	if err != nil {
		return nil, nil, nil, err
	}
	out, err := c.StdoutPipe()
	if err != nil {
		return nil, nil, nil, err
	}
	if err := c.Start(); err != nil {
		return nil, nil, nil, fmt.Errorf("exec tunnel: %w", err)
	}

	m := &execMux{
		stdin: in,
		conns: make(map[string]*execConn),
		ready: make(chan struct{}),
		done:  make(chan struct{}),
	}
	go m.readLoop(c, bufio.NewReaderSize(out, 64<<10))
	stop := func() {
		in.Close() //nolint:errcheck
		_ = c.Process.Kill()
	}

	select {
	case <-m.ready:
	case <-m.done:
		return nil, nil, nil, fmt.Errorf("exec tunnel: relay exited: %s", stderr.String())
	case <-time.After(execReadyTimeout):
		stop()
		return nil, nil, nil, fmt.Errorf("exec tunnel: relay not ready after %s: %s", execReadyTimeout, stderr.String())
	}
	return m.dial, stop, m.done, nil
}

// execConn is one proxied connection: srv is the end the mux drives, and the
// caller holds the other end of the net.Pipe.
type execConn struct {
	id  string
	srv net.Conn
	ack chan bool   // true: connected, false: refused
	out chan []byte // node -> caller; closed by the read loop only
}

type execMux struct {
	stdin io.Writer
	wmu   sync.Mutex // one frame per write
	mu    sync.Mutex // guards conns
	conns map[string]*execConn
	next  atomic.Uint64
	ready chan struct{}
	done  chan struct{}
}

func (m *execMux) send(f frame) {
	m.wmu.Lock()
	defer m.wmu.Unlock()
	_ = writeFrame(m.stdin, f)
}

func (m *execMux) lookup(id string) *execConn {
	m.mu.Lock()
	defer m.mu.Unlock()
	return m.conns[id]
}

func (m *execMux) forget(c *execConn) {
	m.mu.Lock()
	delete(m.conns, c.id)
	m.mu.Unlock()
}

// readLoop delivers frames from the node until the relay exits, then fails every
// open connection. It is the only sender on, and closer of, each conn's out.
func (m *execMux) readLoop(cmd *exec.Cmd, r *bufio.Reader) {
	var readyOnce sync.Once
	for {
		f, err := readFrame(r)
		if err != nil {
			break
		}
		if f.op == "R" {
			readyOnce.Do(func() { close(m.ready) })
			continue
		}
		c := m.lookup(f.id)
		if c == nil {
			continue
		}
		switch f.op {
		case "A":
			c.ack <- true
		case "D":
			c.out <- f.data
		case "C":
			select {
			case c.ack <- false:
			default:
			}
			m.forget(c)
			close(c.out)
		}
	}
	_ = cmd.Wait()
	close(m.done)
	m.mu.Lock()
	open := m.conns
	m.conns = make(map[string]*execConn)
	m.mu.Unlock()
	for _, c := range open {
		select {
		case c.ack <- false:
		default:
		}
		close(c.out)
	}
}

func (m *execMux) dial(ctx context.Context, _, target string) (net.Conn, error) {
	host, portStr, err := net.SplitHostPort(target)
	if err != nil {
		return nil, err
	}
	if host != "127.0.0.1" && host != "localhost" && host != "::1" {
		return nil, fmt.Errorf("exec tunnel reaches only loopback on the node, not %q", host)
	}
	if _, err := strconv.ParseUint(portStr, 10, 16); err != nil {
		return nil, fmt.Errorf("invalid port %q: %w", portStr, err)
	}
	select {
	case <-m.done:
		return nil, errExecEnded
	default:
	}

	cli, srv := net.Pipe()
	c := &execConn{
		id:  strconv.FormatUint(m.next.Add(1), 10),
		srv: srv,
		ack: make(chan bool, 1),
		out: make(chan []byte, execQueueLen),
	}
	m.mu.Lock()
	m.conns[c.id] = c
	m.mu.Unlock()
	m.send(frame{op: "O", id: c.id, arg: portStr})

	fail := func(err error) (net.Conn, error) {
		m.send(frame{op: "C", id: c.id})
		m.forget(c)
		cli.Close() //nolint:errcheck
		srv.Close() //nolint:errcheck
		return nil, err
	}
	select {
	case ok := <-c.ack:
		if !ok {
			return fail(fmt.Errorf("connection to %s refused on the node", target))
		}
	case <-ctx.Done():
		return fail(ctx.Err())
	case <-m.done:
		return fail(errExecEnded)
	}
	go c.writeLoop()
	go m.upLoop(c)
	return cli, nil
}

// writeLoop copies node data to the caller, then closes the pipe once the node
// has closed the connection.
func (c *execConn) writeLoop() {
	for b := range c.out {
		if _, err := c.srv.Write(b); err != nil {
			for range c.out { // caller is gone: discard until the node closes
			}
			break
		}
	}
	c.srv.Close() //nolint:errcheck
}

// upLoop copies the caller's data to the node until the caller closes.
func (m *execMux) upLoop(c *execConn) {
	buf := make([]byte, execChunk)
	for {
		n, err := c.srv.Read(buf)
		if n > 0 {
			m.send(frame{op: "D", id: c.id, data: buf[:n]})
		}
		if err != nil {
			break
		}
	}
	m.send(frame{op: "C", id: c.id}) // the node answers "C" once the service closes
}

// tailBuffer keeps the last execStderrKeep bytes written to it.
type tailBuffer struct {
	mu  sync.Mutex
	buf []byte
}

func (t *tailBuffer) Write(p []byte) (int, error) {
	t.mu.Lock()
	defer t.mu.Unlock()
	t.buf = append(t.buf, p...)
	if len(t.buf) > execStderrKeep {
		t.buf = t.buf[len(t.buf)-execStderrKeep:]
	}
	return len(p), nil
}

func (t *tailBuffer) String() string {
	t.mu.Lock()
	defer t.mu.Unlock()
	return strings.TrimSpace(string(t.buf))
}
