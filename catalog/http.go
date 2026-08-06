package catalog

import (
	"bytes"
	"compress/gzip"
	"context"
	"encoding/json"
	"fmt"
	"io"
	"net/http"
	"strings"
	"time"
)

// maxIndexSize bounds what a fetch will read, so a wrong URL cannot pull an
// unbounded body into memory.
const maxIndexSize = 32 << 20

// httpBackend reads a source over HTTP, where there is no directory listing and
// no cheap way to read a hundred recipes — hence the generated index.
type httpBackend struct {
	src    *Source
	cache  Cache
	client *http.Client
}

func (h *httpBackend) read(ctx context.Context, path string) ([]byte, error) {
	if data, fresh := h.cache.get(h.src.Base, path); fresh {
		return data, nil
	}

	data, err := h.fetch(ctx, path)
	if err == nil {
		_ = h.cache.put(h.src.Base, path, data)
		return data, nil
	}

	// Expired with no route out is the normal state of a compute node, not an
	// error. Serve what is cached and say that it happened.
	if stale, _ := h.cache.get(h.src.Base, path); stale != nil {
		h.src.Stale = true
		return stale, nil
	}
	return nil, err
}

func (h *httpBackend) fetch(ctx context.Context, path string) ([]byte, error) {
	url := strings.TrimSuffix(h.src.Base, "/") + "/" + path
	req, err := http.NewRequestWithContext(ctx, http.MethodGet, url, nil)
	if err != nil {
		return nil, err
	}
	resp, err := h.httpClient().Do(req)
	if err != nil {
		return nil, err
	}
	defer resp.Body.Close()
	if resp.StatusCode != http.StatusOK {
		return nil, fmt.Errorf("catalog: %s: %s", url, resp.Status)
	}
	return io.ReadAll(io.LimitReader(resp.Body, maxIndexSize))
}

func (h *httpBackend) httpClient() *http.Client {
	if h.client != nil {
		return h.client
	}
	return &http.Client{Timeout: 60 * time.Second}
}

// entries reads the generated index, preferring the compressed copy.
func (h *httpBackend) entries(ctx context.Context) (map[string]*Entry, error) {
	data, err := h.readIndex(ctx, indexDir+"/recipes.json")
	if err != nil {
		return nil, err
	}
	var raw map[string]*Entry
	if err := json.Unmarshal(data, &raw); err != nil {
		return nil, fmt.Errorf("catalog: %s index: %w", h.src.Name, err)
	}
	for name, e := range raw {
		e.Name = name
	}
	return raw, nil
}

// readIndex tries path.gz first, since that is what a collection publishes for
// a fetch, and falls back to the plain file.
func (h *httpBackend) readIndex(ctx context.Context, path string) ([]byte, error) {
	gz, gzErr := h.read(ctx, path+".gz")
	if gzErr == nil {
		if data, err := gunzip(gz); err == nil {
			return data, nil
		}
	}
	data, err := h.read(ctx, path)
	if err != nil && gzErr != nil {
		return nil, gzErr
	}
	return data, err
}

func gunzip(data []byte) ([]byte, error) {
	zr, err := gzip.NewReader(bytes.NewReader(data))
	if err != nil {
		return nil, err
	}
	defer zr.Close()
	return io.ReadAll(io.LimitReader(zr, maxIndexSize))
}
