package proxy

import (
	"bufio"
	"context"
	"fmt"
	"io"
	"net"
	"net/http"
	"sync"
)

var brPool = sync.Pool{
	New: func() any { return bufio.NewReader(nil) },
}

type closeWriter interface{ CloseWrite() error }

// DialFunc dials a network connection to addr. Compatible with
// net.Dialer.DialContext and http.Transport.DialContext.
type DialFunc func(ctx context.Context, network, addr string) (net.Conn, error)

// RunProxy listens on listenAddr and serves both SOCKS5 and HTTP CONNECT,
// forwarding connections through dial. Stops when ctx is cancelled.
func RunProxy(ctx context.Context, listenAddr string, dial DialFunc) error {
	ln, err := net.Listen("tcp", listenAddr)
	if err != nil {
		return fmt.Errorf("proxy listen on %s: %w", listenAddr, err)
	}
	go func() {
		<-ctx.Done()
		ln.Close()
	}()
	for {
		conn, err := ln.Accept()
		if err != nil {
			select {
			case <-ctx.Done():
				return nil
			default:
				return err
			}
		}
		go serveConn(ctx, conn, dial)
	}
}

// serveConn peeks the first byte to distinguish SOCKS5 (0x05) from HTTP CONNECT.
func serveConn(ctx context.Context, conn net.Conn, dial DialFunc) {
	defer conn.Close()
	br := brPool.Get().(*bufio.Reader)
	br.Reset(conn)
	defer func() { br.Reset(nil); brPool.Put(br) }()
	b, err := br.Peek(1)
	if err != nil {
		return
	}
	if b[0] == 0x05 {
		handleSOCKS5(ctx, br, conn, dial)
	} else {
		handleHTTPConnect(ctx, br, conn, dial)
	}
}

// handleSOCKS5 implements RFC 1928 SOCKS5: no-auth negotiation + CONNECT command.
func handleSOCKS5(ctx context.Context, br *bufio.Reader, conn net.Conn, dial DialFunc) {
	// Greeting: {ver=5, nMethods, methods...} → reply {5, 0} (no auth)
	hdr := make([]byte, 2)
	if _, err := io.ReadFull(br, hdr); err != nil || hdr[0] != 0x05 {
		return
	}
	methods := make([]byte, hdr[1])
	if _, err := io.ReadFull(br, methods); err != nil {
		return
	}
	conn.Write([]byte{0x05, 0x00}) //nolint:errcheck

	// Request: {ver, cmd=1, rsv, atyp, addr..., port(2)}
	req := make([]byte, 4)
	if _, err := io.ReadFull(br, req); err != nil || req[0] != 0x05 || req[1] != 0x01 {
		conn.Write([]byte{0x05, 0x07, 0x00, 0x01, 0, 0, 0, 0, 0, 0}) //nolint:errcheck
		return
	}

	var host string
	switch req[3] {
	case 0x01: // IPv4
		addr := make([]byte, 4)
		if _, err := io.ReadFull(br, addr); err != nil {
			return
		}
		host = net.IP(addr).String()
	case 0x03: // domain
		lenB := make([]byte, 1)
		if _, err := io.ReadFull(br, lenB); err != nil {
			return
		}
		name := make([]byte, lenB[0])
		if _, err := io.ReadFull(br, name); err != nil {
			return
		}
		host = string(name)
	case 0x04: // IPv6
		addr := make([]byte, 16)
		if _, err := io.ReadFull(br, addr); err != nil {
			return
		}
		host = "[" + net.IP(addr).String() + "]"
	default:
		conn.Write([]byte{0x05, 0x08, 0x00, 0x01, 0, 0, 0, 0, 0, 0}) //nolint:errcheck
		return
	}

	portB := make([]byte, 2)
	if _, err := io.ReadFull(br, portB); err != nil {
		return
	}
	target := fmt.Sprintf("%s:%d", host, int(portB[0])<<8|int(portB[1]))

	upstream, err := dial(ctx, "tcp", target)
	if err != nil {
		conn.Write([]byte{0x05, 0x04, 0x00, 0x01, 0, 0, 0, 0, 0, 0}) //nolint:errcheck
		return
	}
	defer upstream.Close()

	conn.Write([]byte{0x05, 0x00, 0x00, 0x01, 0, 0, 0, 0, 0, 0}) //nolint:errcheck
	relay(conn, upstream)
}

// handleHTTPConnect parses an HTTP CONNECT request and tunnels through dial.
func handleHTTPConnect(ctx context.Context, br *bufio.Reader, conn net.Conn, dial DialFunc) {
	req, err := http.ReadRequest(br)
	if err != nil {
		return
	}
	if req.Method != http.MethodConnect {
		fmt.Fprintf(conn, "HTTP/1.1 405 Method Not Allowed\r\n\r\n") //nolint:errcheck
		return
	}

	upstream, err := dial(ctx, "tcp", req.Host)
	if err != nil {
		fmt.Fprintf(conn, "HTTP/1.1 502 Bad Gateway\r\n\r\n") //nolint:errcheck
		return
	}
	defer upstream.Close()

	fmt.Fprintf(conn, "HTTP/1.1 200 Connection Established\r\n\r\n") //nolint:errcheck
	relay(conn, upstream)
}

func relay(a, b net.Conn) {
	done := make(chan struct{}, 2)
	go func() {
		io.Copy(b, a) //nolint:errcheck
		if cw, ok := b.(closeWriter); ok {
			cw.CloseWrite() //nolint:errcheck
		}
		done <- struct{}{}
	}()
	go func() {
		io.Copy(a, b) //nolint:errcheck
		if cw, ok := a.(closeWriter); ok {
			cw.CloseWrite() //nolint:errcheck
		}
		done <- struct{}{}
	}()
	<-done
	<-done
}
