// Package server implements the condatainer dashboard HTTP server.
// It runs on the login node, provides a web dashboard, reverse-proxies to
// helper services on compute nodes via SSH tunnels, and streams SSE logs.
package server

import (
	"context"
	"encoding/json"
	"fmt"
	"log/slog"
	"net"
	"net/http"
	"os"
	"os/signal"
	"path/filepath"
	"sync"
	"syscall"
	"time"

	"github.com/Justype/condatainer/internal/config"
	"github.com/Justype/condatainer/internal/logging"
	"github.com/Justype/condatainer/internal/utils"
)

// taskEntry tracks an in-progress background task (overlay creation, etc.).
type taskEntry struct {
	broker *SSEBroker
	cancel context.CancelFunc
}

// srv is the central server instance shared by all handlers.
type srv struct {
	port      int
	startTime time.Time
	ctx       context.Context
	brokers   *brokerRegistry
	proxies   *proxyRegistry
	watcher   *watcher
	tasks     sync.Map // taskID → *taskEntry
}

// RunDaemon starts the HTTP server, writes the PID file, signals readiness via
// reportPipe (if non-nil), then blocks until SIGTERM/SIGINT.
// If watchPID > 0, the server also exits when that process dies (used in
// non-daemon mode to tie the server lifetime to the parent shell).
func RunDaemon(port int, reportPipe *os.File, watchPID int) error {
	report := func(msg string) {
		if reportPipe == nil {
			return
		}
		if msg != "" {
			fmt.Fprint(reportPipe, msg) //nolint:errcheck
		}
		reportPipe.Close()
		reportPipe = nil
	}

	ln, err := net.Listen("tcp", fmt.Sprintf("127.0.0.1:%d", port))
	if err != nil {
		report(fmt.Sprintf("listen on port %d: %v", port, err))
		return err
	}
	actualPort := ln.Addr().(*net.TCPAddr).Port
	ctx, cancel := context.WithCancel(context.Background())
	defer cancel()
	ctx = logging.WithLogger(ctx, slog.New(slog.NewTextHandler(os.Stderr, &slog.HandlerOptions{
		Level: slog.LevelDebug,
	})))

	logger := logging.FromContext(ctx)
	s := &srv{
		port:      actualPort,
		startTime: time.Now(),
		ctx:       ctx,
		brokers:   newBrokerRegistry(),
		proxies:   newProxyRegistry(logger),
	}
	s.watcher = newWatcher(s)

	mux := http.NewServeMux()
	registerRoutes(mux, s)
	httpServer := &http.Server{Handler: s.hostDispatch(mux)}

	// Write PID file.
	host, _ := os.Hostname()
	pidFile := config.GetServerPidFilePath()
	if pidFile != "" {
		if err := writeState(pidFile, config.ServerState{
			Host:   host,
			Port:   actualPort,
			PID:    os.Getpid(),
			Daemon: watchPID == 0,
		}); err != nil {
			report(fmt.Sprintf("writing PID file: %v", err))
			return err
		}
		defer os.Remove(pidFile)
	}

	// Signal parent: ready.
	report("")

	// Non-daemon: exit when the parent shell dies.
	if watchPID > 0 {
		go func() {
			t := time.NewTicker(5 * time.Second)
			defer t.Stop()
			for {
				select {
				case <-ctx.Done():
					return
				case <-t.C:
					p, err := os.FindProcess(watchPID)
					if err != nil || p.Signal(syscall.Signal(0)) != nil {
						logger.Debug("server: parent shell gone, shutting down", "pid", watchPID)
						syscall.Kill(os.Getpid(), syscall.SIGTERM) //nolint:errcheck
						return
					}
				}
			}
		}()
	}

	// Start NFS watcher.
	go s.watcher.run(ctx)

	// Handle shutdown signals.
	sig := make(chan os.Signal, 1)
	signal.Notify(sig, syscall.SIGTERM, syscall.SIGINT)
	go func() {
		<-sig
		cancel()
		shutCtx, shutCancel := context.WithTimeout(context.Background(), 5*time.Second)
		defer shutCancel()
		httpServer.Shutdown(shutCtx) //nolint:errcheck
	}()

	logging.FromContext(ctx).Debug("server: listening", "addr", fmt.Sprintf("http://127.0.0.1:%d", actualPort))
	if err := httpServer.Serve(ln); err != nil && err != http.ErrServerClosed {
		return err
	}
	return nil
}

// writeState writes the server State JSON to the PID file.
func writeState(path string, s config.ServerState) error {
	data, err := json.Marshal(s)
	if err != nil {
		return err
	}
	if !utils.EnsureWritableDir(filepath.Dir(path)) {
		return fmt.Errorf("cannot create state directory: %s", filepath.Dir(path))
	}
	return os.WriteFile(path, data, utils.PermFile)
}

// ReadState reads and parses the server State from the PID file.
func ReadState(pidFile string) (*config.ServerState, error) {
	data, err := os.ReadFile(pidFile)
	if err != nil {
		return nil, err
	}
	var s config.ServerState
	if err := json.Unmarshal(data, &s); err != nil {
		return nil, err
	}
	return &s, nil
}

// IsAlive returns true if the server described by the PID file is running.
func IsAlive(pidFile string) bool {
	s, err := ReadState(pidFile)
	if err != nil || s.PID == 0 {
		return false
	}
	p, err := os.FindProcess(s.PID)
	if err != nil {
		return false
	}
	return p.Signal(syscall.Signal(0)) == nil
}
