package proxy

import (
	"bufio"
	"bytes"
	"context"
	"crypto/rand"
	"io"
	"net"
	"os"
	"strings"
	"sync"
	"testing"
	"time"
)

// The test binary doubles as the relay command.
func TestMain(m *testing.M) {
	if os.Getenv("CNT_TEST_EXEC_RELAY") == "1" {
		if RunExecRelay(os.Stdin, os.Stdout) != nil {
			os.Exit(1)
		}
		os.Exit(0)
	}
	os.Exit(m.Run())
}

func echoServer(t *testing.T) string {
	t.Helper()
	ln, err := net.Listen("tcp", "127.0.0.1:0")
	if err != nil {
		t.Fatal(err)
	}
	t.Cleanup(func() { ln.Close() })
	go func() {
		for {
			c, err := ln.Accept()
			if err != nil {
				return
			}
			go func() { defer c.Close(); io.Copy(c, c) }() //nolint:errcheck
		}
	}()
	return ln.Addr().String()
}

func TestDialViaExec(t *testing.T) {
	addr := echoServer(t)
	dial, stop, done, err := DialViaExec([]string{"env", "CNT_TEST_EXEC_RELAY=1", os.Args[0]})
	if err != nil {
		t.Fatal(err)
	}
	defer stop()
	ctx, cancel := context.WithTimeout(context.Background(), 20*time.Second)
	defer cancel()

	// Binary frames with no newline echo back promptly, request after request.
	c, err := dial(ctx, "tcp", addr)
	if err != nil {
		t.Fatal(err)
	}
	for i := 0; i < 50; i++ {
		msg := make([]byte, 20)
		rand.Read(msg)                                 //nolint:errcheck
		c.SetDeadline(time.Now().Add(5 * time.Second)) //nolint:errcheck
		if _, err := c.Write(msg); err != nil {
			t.Fatal(err)
		}
		got := make([]byte, len(msg))
		if _, err := io.ReadFull(c, got); err != nil || !bytes.Equal(got, msg) {
			t.Fatalf("round trip %d: err=%v equal=%v", i, err, bytes.Equal(got, msg))
		}
	}
	c.Close()

	// Concurrent bulk transfers share the one relay.
	var wg sync.WaitGroup
	for i := 0; i < 4; i++ {
		wg.Add(1)
		go func() {
			defer wg.Done()
			c, err := dial(ctx, "tcp", addr)
			if err != nil {
				t.Error(err)
				return
			}
			defer c.Close()
			data := make([]byte, 200000)
			rand.Read(data)  //nolint:errcheck
			go c.Write(data) //nolint:errcheck
			got := make([]byte, len(data))
			c.SetDeadline(time.Now().Add(15 * time.Second)) //nolint:errcheck
			if _, err := io.ReadFull(c, got); err != nil || !bytes.Equal(got, data) {
				t.Errorf("bulk echo: err=%v equal=%v", err, bytes.Equal(got, data))
			}
		}()
	}
	wg.Wait()

	// A port nothing listens on is refused, not hung.
	ln, _ := net.Listen("tcp", "127.0.0.1:0")
	closed := ln.Addr().String()
	ln.Close()
	if _, err := dial(ctx, "tcp", closed); err == nil {
		t.Fatal("dial to a closed port succeeded")
	}
	if _, err := dial(ctx, "tcp", "10.0.0.1:80"); err == nil {
		t.Fatal("dial to a non-loopback host succeeded")
	}

	stop()
	select {
	case <-done:
	case <-time.After(5 * time.Second):
		t.Fatal("done not closed after stop")
	}
	if _, err := dial(ctx, "tcp", addr); err == nil {
		t.Fatal("dial after stop succeeded")
	}
}

// A relay that cannot start is reported with what the command printed.
func TestDialViaExecReportsRelayFailure(t *testing.T) {
	_, _, _, err := DialViaExec([]string{"sh", "-c", "echo boom >&2; exit 3"})
	if err == nil || !strings.Contains(err.Error(), "boom") {
		t.Fatalf("err = %v, want the command's message", err)
	}
}

// Data frames carry any bytes, including newlines and NULs; a stray line surfaces as an unhandled op.
func TestExecFrameRoundTrip(t *testing.T) {
	payload := []byte("a\nb\x00\r\n\xff\n")
	var buf bytes.Buffer
	buf.WriteString("warning: not a frame\n")
	for _, f := range []frame{{op: "O", id: "1", arg: "8080"}, {op: "D", id: "1", data: payload}, {op: "C", id: "1"}} {
		if err := writeFrame(&buf, f); err != nil {
			t.Fatal(err)
		}
	}
	r := bufio.NewReader(&buf)
	f, _ := readFrame(r)
	if f.op != "warning:" { // a non-frame line parses as an unknown op and is ignored by callers
		t.Fatalf("first = %+v", f)
	}
	for _, want := range []frame{{op: "O", id: "1", arg: "8080"}, {op: "D", id: "1", arg: "8", data: payload}, {op: "C", id: "1"}} {
		got, err := readFrame(r)
		if err != nil || got.op != want.op || got.id != want.id || got.arg != want.arg || !bytes.Equal(got.data, want.data) {
			t.Fatalf("got %+v, %v; want %+v", got, err, want)
		}
	}
}
