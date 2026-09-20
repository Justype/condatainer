package cmd

import (
	"context"
	"net"
	"testing"
	"time"
)

func TestWaitForPort(t *testing.T) {
	ln, err := net.Listen("tcp", "127.0.0.1:0")
	if err != nil {
		t.Fatal(err)
	}
	port := ln.Addr().(*net.TCPAddr).Port

	if !waitForPort(context.Background(), port, 2*time.Second) {
		t.Fatal("did not see a listening port")
	}

	ln.Close()
	if waitForPort(context.Background(), port, 400*time.Millisecond) {
		t.Fatal("reported a closed port as open")
	}
}
