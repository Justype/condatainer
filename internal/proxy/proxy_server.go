package proxy

import (
	"bufio"
	"context"
	"encoding/binary"
	"fmt"
	"io"
	"net"
	"net/http"
	"strconv"
)

// RunProxy listens on listenAddr and serves both SOCKS5 and HTTP CONNECT,
// forwarding connections through the upstream SOCKS5 proxy at socks5Addr.
// Stops when ctx is cancelled.
func RunProxy(ctx context.Context, listenAddr, socks5Addr string) error {
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
		go serveConn(conn, socks5Addr)
	}
}

// serveConn peeks the first byte to distinguish SOCKS5 (0x05) from HTTP CONNECT.
func serveConn(conn net.Conn, socks5Addr string) {
	defer conn.Close()
	br := bufio.NewReader(conn)
	b, err := br.Peek(1)
	if err != nil {
		return
	}
	if b[0] == 0x05 {
		handleSOCKS5(br, conn, socks5Addr)
	} else {
		handleHTTPConnect(br, conn, socks5Addr)
	}
}

// handleSOCKS5 implements RFC 1928 SOCKS5: no-auth negotiation + CONNECT command.
// Dials the target via the upstream SOCKS5 proxy.
func handleSOCKS5(br *bufio.Reader, conn net.Conn, socks5Addr string) {
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
		host = net.IP(addr).String()
	default:
		conn.Write([]byte{0x05, 0x08, 0x00, 0x01, 0, 0, 0, 0, 0, 0}) //nolint:errcheck
		return
	}

	portB := make([]byte, 2)
	if _, err := io.ReadFull(br, portB); err != nil {
		return
	}
	port := binary.BigEndian.Uint16(portB)
	target := fmt.Sprintf("%s:%d", host, port)

	upstream, err := dialViaSocks5(socks5Addr, target)
	if err != nil {
		conn.Write([]byte{0x05, 0x04, 0x00, 0x01, 0, 0, 0, 0, 0, 0}) //nolint:errcheck
		return
	}
	defer upstream.Close()

	conn.Write([]byte{0x05, 0x00, 0x00, 0x01, 0, 0, 0, 0, 0, 0}) //nolint:errcheck
	relay(conn, upstream)
}

// handleHTTPConnect parses an HTTP CONNECT request and tunnels through the upstream SOCKS5.
func handleHTTPConnect(br *bufio.Reader, conn net.Conn, socks5Addr string) {
	req, err := http.ReadRequest(br)
	if err != nil {
		return
	}
	if req.Method != http.MethodConnect {
		fmt.Fprintf(conn, "HTTP/1.1 405 Method Not Allowed\r\n\r\n") //nolint:errcheck
		return
	}

	upstream, err := dialViaSocks5(socks5Addr, req.Host)
	if err != nil {
		fmt.Fprintf(conn, "HTTP/1.1 502 Bad Gateway\r\n\r\n") //nolint:errcheck
		return
	}
	defer upstream.Close()

	fmt.Fprintf(conn, "HTTP/1.1 200 Connection Established\r\n\r\n") //nolint:errcheck
	relay(conn, upstream)
}

// dialViaSocks5 connects to target ("host:port") through the SOCKS5 proxy at socks5Addr.
func dialViaSocks5(socks5Addr, target string) (net.Conn, error) {
	host, portStr, err := net.SplitHostPort(target)
	if err != nil {
		return nil, err
	}
	p, err := strconv.ParseUint(portStr, 10, 16)
	if err != nil {
		return nil, fmt.Errorf("invalid port %q: %w", portStr, err)
	}
	port := uint16(p)

	conn, err := net.Dial("tcp", socks5Addr)
	if err != nil {
		return nil, err
	}

	// Greeting
	if _, err := conn.Write([]byte{0x05, 0x01, 0x00}); err != nil {
		conn.Close()
		return nil, err
	}
	resp := make([]byte, 2)
	if _, err := io.ReadFull(conn, resp); err != nil || resp[0] != 0x05 || resp[1] != 0x00 {
		conn.Close()
		return nil, fmt.Errorf("socks5 handshake failed")
	}

	// CONNECT request with domain name
	req := []byte{0x05, 0x01, 0x00, 0x03, byte(len(host))}
	req = append(req, []byte(host)...)
	req = append(req, byte(port>>8), byte(port))
	if _, err := conn.Write(req); err != nil {
		conn.Close()
		return nil, err
	}

	// Response header (4 bytes minimum)
	hdr := make([]byte, 4)
	if _, err := io.ReadFull(conn, hdr); err != nil || hdr[1] != 0x00 {
		conn.Close()
		return nil, fmt.Errorf("socks5 CONNECT failed: status %d", hdr[1])
	}
	// Skip bound address in response
	switch hdr[3] {
	case 0x01:
		io.ReadFull(conn, make([]byte, 6)) //nolint:errcheck
	case 0x03:
		l := make([]byte, 1)
		io.ReadFull(conn, l)                         //nolint:errcheck
		io.ReadFull(conn, make([]byte, int(l[0])+2)) //nolint:errcheck
	case 0x04:
		io.ReadFull(conn, make([]byte, 18)) //nolint:errcheck
	}
	return conn, nil
}

func relay(a, b io.ReadWriter) {
	done := make(chan struct{}, 2)
	go func() { io.Copy(b, a); done <- struct{}{} }() //nolint:errcheck
	go func() { io.Copy(a, b); done <- struct{}{} }() //nolint:errcheck
	<-done
	<-done
}
