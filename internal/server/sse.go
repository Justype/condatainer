package server

import (
	"context"
	"encoding/json"
	"fmt"
	"net/http"
	"sync"
)

// SSEBroker fan-outs SSE messages to subscribed HTTP clients.
// One broker per helper ID; the watcher pushes messages into it.
// finished/finalMsg cache the terminal "done" event (see publishFinal) so
// a client that connects after the task finished still gets the result.
type SSEBroker struct {
	mu       sync.Mutex
	clients  map[chan []byte]struct{}
	finished bool
	finalMsg []byte
}

func newSSEBroker() *SSEBroker {
	return &SSEBroker{clients: make(map[chan []byte]struct{})}
}

// subscribe registers ch and, atomically with that registration, reports
// whether the task already finished along with its cached final message.
func (b *SSEBroker) subscribe() (ch chan []byte, finalMsg []byte, finished bool) {
	ch = make(chan []byte, 32)
	b.mu.Lock()
	defer b.mu.Unlock()
	b.clients[ch] = struct{}{}
	return ch, b.finalMsg, b.finished
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

// publishFinal is like publish, but also caches data as the terminal
// message so later subscribers immediately learn the task is done.
func (b *SSEBroker) publishFinal(data []byte) {
	b.mu.Lock()
	defer b.mu.Unlock()
	b.finished = true
	b.finalMsg = data
	for ch := range b.clients {
		select {
		case ch <- data:
		default:
		}
	}
}

// ServeSSE writes SSE events to the client until r.Context() is done, or
// immediately writes the cached final event and returns if the task had
// already finished before this client connected.
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

	ch, finalMsg, finished := b.subscribe()
	defer b.unsubscribe(ch)

	if finished {
		fmt.Fprintf(w, "data: %s\n\n", finalMsg) //nolint:errcheck
		flusher.Flush()
		return
	}

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
	b.publishFinal(data)
}

// broadcastCancelled sends a terminal "done" event marked cancelled:true,
// so the client can show "Cancelled" rather than an error.
func broadcastCancelled(b *SSEBroker) {
	data, _ := json.Marshal(map[string]interface{}{"t": "done", "ok": false, "cancelled": true})
	b.publishFinal(data)
}

// broadcastResult sends a task's terminal "done" event, reporting it as a
// cancellation rather than an error when the failure came from ctx being
// cancelled.
func broadcastResult(b *SSEBroker, ctx context.Context, err error) {
	if err != nil && ctx.Err() != nil {
		broadcastCancelled(b)
		return
	}
	broadcastDone(b, err)
}

// broadcastProgress sends a {"t":"progress","current":n,"total":n} event
// to all SSE subscribers of b.
func broadcastProgress(b *SSEBroker, current, total int) {
	data, _ := json.Marshal(map[string]interface{}{"t": "progress", "current": current, "total": total})
	b.publish(data)
}
