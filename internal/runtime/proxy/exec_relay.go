package proxy

import (
	"bufio"
	"io"
	"net"
	"sync"
	"time"
)

type relayConn struct {
	c    net.Conn
	in   chan []byte   // caller -> service; closed by the frame loop on "C"
	gone chan struct{} // closed when the service side has ended
}

// RunExecRelay is the node side of DialViaExec: it reads frames from in, keeps
// one loopback connection per id, and writes the service's data to out. It
// returns when in closes.
func RunExecRelay(in io.Reader, out io.Writer) error {
	var wmu sync.Mutex
	send := func(f frame) {
		wmu.Lock()
		defer wmu.Unlock()
		_ = writeFrame(out, f)
	}
	var mu sync.Mutex
	conns := map[string]*relayConn{}
	lookup := func(id string) *relayConn {
		mu.Lock()
		defer mu.Unlock()
		return conns[id]
	}

	send(frame{op: "R"})
	r := bufio.NewReaderSize(in, 64<<10)
	for {
		f, err := readFrame(r)
		if err != nil {
			break
		}
		id := f.id
		if id == "" {
			continue
		}
		switch f.op {
		case "O":
			nc, err := net.DialTimeout("tcp", net.JoinHostPort("127.0.0.1", f.arg), 10*time.Second)
			if err != nil {
				send(frame{op: "C", id: id})
				continue
			}
			rc := &relayConn{c: nc, in: make(chan []byte, execQueueLen), gone: make(chan struct{})}
			mu.Lock()
			conns[id] = rc
			mu.Unlock()
			send(frame{op: "A", id: id})
			go rc.writeLoop()
			go func() {
				buf := make([]byte, execChunk)
				for {
					n, err := nc.Read(buf)
					if n > 0 {
						send(frame{op: "D", id: id, data: buf[:n]})
					}
					if err != nil {
						break
					}
				}
				send(frame{op: "C", id: id})
				mu.Lock()
				if conns[id] == rc {
					delete(conns, id)
				}
				mu.Unlock()
				close(rc.gone)
				nc.Close() //nolint:errcheck
			}()
		case "D":
			rc := lookup(id)
			if rc == nil {
				continue
			}
			select {
			case rc.in <- f.data:
			case <-rc.gone:
			}
		case "C":
			mu.Lock()
			rc := conns[id]
			delete(conns, id)
			mu.Unlock()
			if rc != nil {
				close(rc.in)
			}
		}
	}

	mu.Lock()
	open := conns
	conns = map[string]*relayConn{}
	mu.Unlock()
	for _, rc := range open {
		close(rc.in)
		rc.c.Close() //nolint:errcheck
	}
	return nil
}

// writeLoop copies the caller's data to the service, then half-closes so the
// service sees the caller's end of the stream.
func (rc *relayConn) writeLoop() {
	for {
		select {
		case b, ok := <-rc.in:
			if !ok {
				if tc, ok := rc.c.(*net.TCPConn); ok {
					tc.CloseWrite() //nolint:errcheck
				}
				return
			}
			if _, err := rc.c.Write(b); err != nil {
				return
			}
		case <-rc.gone:
			return
		}
	}
}
