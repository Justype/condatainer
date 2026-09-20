package cmd

import (
	"context"
	"net"
	"strings"
	"testing"
	"time"

	"github.com/spf13/cobra"
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

func TestServerControlRefusedInsideContainer(t *testing.T) {
	t.Setenv("IN_CONDATAINER", "1")
	for _, c := range []*cobra.Command{serverStartCmd, serverStopCmd, serverRestartCmd} {
		err := requireLoginNode(c)
		if err == nil || !strings.Contains(err.Error(), "cannot run inside a container") {
			t.Errorf("server %s inside a container: err = %v", c.Name(), err)
		}
	}
}
