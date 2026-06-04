package server

import (
	"encoding/json"
	"fmt"
	"net/http"
	"sync"
)

// SSEBroker fan-outs SSE messages to subscribed HTTP clients.
// One broker per helper ID; the watcher pushes messages into it.
type SSEBroker struct {
	mu      sync.Mutex
	clients map[chan []byte]struct{}
}

func newSSEBroker() *SSEBroker {
	return &SSEBroker{clients: make(map[chan []byte]struct{})}
}

func (b *SSEBroker) subscribe() chan []byte {
	ch := make(chan []byte, 32)
	b.mu.Lock()
	b.clients[ch] = struct{}{}
	b.mu.Unlock()
	return ch
}

func (b *SSEBroker) unsubscribe(ch chan []byte) {
	b.mu.Lock()
	delete(b.clients, ch)
	b.mu.Unlock()
}

func (b *SSEBroker) publish(data []byte) {
	b.mu.Lock()
	defer b.mu.Unlock()
	for ch := range b.clients {
		select {
		case ch <- data:
		default: // drop if client is slow
		}
	}
}

// ServeSSE writes SSE events to the client until r.Context() is done.
func (b *SSEBroker) ServeSSE(w http.ResponseWriter, r *http.Request) {
	flusher, ok := w.(http.Flusher)
	if !ok {
		http.Error(w, "streaming not supported", http.StatusInternalServerError)
		return
	}
	w.Header().Set("Content-Type", "text/event-stream")
	w.Header().Set("Cache-Control", "no-cache")
	w.Header().Set("Connection", "keep-alive")
	w.Header().Set("X-Accel-Buffering", "no")

	ch := b.subscribe()
	defer b.unsubscribe(ch)

	for {
		select {
		case <-r.Context().Done():
			return
		case data := <-ch:
			fmt.Fprintf(w, "data: %s\n\n", data) //nolint:errcheck
			flusher.Flush()
		}
	}
}

// brokerRegistry holds one SSEBroker per helper ID.
type brokerRegistry struct {
	mu      sync.Mutex
	brokers map[string]*SSEBroker
}

func newBrokerRegistry() *brokerRegistry {
	return &brokerRegistry{brokers: make(map[string]*SSEBroker)}
}

func (r *brokerRegistry) get(id string) *SSEBroker {
	r.mu.Lock()
	defer r.mu.Unlock()
	if b, ok := r.brokers[id]; ok {
		return b
	}
	b := newSSEBroker()
	r.brokers[id] = b
	return b
}

func (r *brokerRegistry) publish(id string, level, text string) {
	data, _ := json.Marshal(map[string]string{"level": level, "text": text})
	r.get(id).publish(data)
}

// brokerWriter wraps an SSEBroker as an io.Writer.
// Each Write call broadcasts the raw chunk as a JSON task output message.
type brokerWriter struct{ b *SSEBroker }

func (w *brokerWriter) Write(p []byte) (int, error) {
	data, _ := json.Marshal(map[string]string{"t": "out", "d": string(p)})
	w.b.publish(data)
	return len(p), nil
}

// broadcastDone sends a terminal "done" event to all SSE subscribers of b.
func broadcastDone(b *SSEBroker, err error) {
	msg := map[string]interface{}{"t": "done", "ok": err == nil}
	if err != nil {
		msg["err"] = err.Error()
	}
	data, _ := json.Marshal(msg)
	b.publish(data)
}
