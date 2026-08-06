package catalog

import (
	"net/http"
	"net/http/httptest"
	"os"
	"path/filepath"
	"testing"
	"time"
)

func TestCacheRoundTrip(t *testing.T) {
	c := Cache{Dir: t.TempDir(), TTL: time.Minute}
	if err := c.put("https://a.invalid", "index/recipes.json", []byte("hello")); err != nil {
		t.Fatal(err)
	}

	data, fresh := c.get("https://a.invalid", "index/recipes.json")
	if string(data) != "hello" || !fresh {
		t.Errorf("get = %q fresh=%v", data, fresh)
	}
	// Two bases cannot collide.
	if data, _ := c.get("https://b.invalid", "index/recipes.json"); data != nil {
		t.Errorf("a different base read %q", data)
	}
	// A disabled cache stores nothing and reports nothing.
	var off Cache
	if err := off.put("https://a.invalid", "x", []byte("y")); err != nil {
		t.Fatal(err)
	}
	if data, _ := off.get("https://a.invalid", "x"); data != nil {
		t.Error("disabled cache returned data")
	}
}

func TestCacheExpiry(t *testing.T) {
	c := Cache{Dir: t.TempDir(), TTL: time.Hour}
	if err := c.put("base", "p", []byte("v")); err != nil {
		t.Fatal(err)
	}
	old := time.Now().Add(-2 * time.Hour)
	if err := os.Chtimes(c.file("base", "p"), old, old); err != nil {
		t.Fatal(err)
	}
	// Expired bytes are still returned, so a failed fetch can serve them.
	data, fresh := c.get("base", "p")
	if string(data) != "v" || fresh {
		t.Errorf("get = %q fresh=%v, want the stale copy", data, fresh)
	}
}

// An expired cache with no route out is the normal state of a compute node.
// The cached copy is served and the source says it happened.
func TestStaleBeatsFailing(t *testing.T) {
	root := fakeSource(t)
	srv := httptest.NewServer(http.FileServer(http.Dir(root)))
	cache := Cache{Dir: t.TempDir(), TTL: time.Hour}

	cat, err := Open(t.Context(), []Spec{{Name: "remote", Base: srv.URL}}, cache)
	if err != nil {
		t.Fatal(err)
	}
	if len(cat.Entries(t.Context())) == 0 {
		t.Fatal("a fresh fetch returned nothing")
	}
	if cat[0].Stale {
		t.Error("a fresh fetch should not be stale")
	}

	// Expire everything cached for this source, then cut the network.
	dir := filepath.Dir(cache.file(srv.URL, "index/recipes.json"))
	old := time.Now().Add(-2 * time.Hour)
	entries, err := os.ReadDir(dir)
	if err != nil {
		t.Fatal(err)
	}
	for _, e := range entries {
		if err := os.Chtimes(filepath.Join(dir, e.Name()), old, old); err != nil {
			t.Fatal(err)
		}
	}
	srv.Close()

	offline, err := Open(t.Context(), []Spec{{Name: "remote", Base: srv.URL}}, cache)
	if err != nil {
		t.Fatal(err)
	}
	got := offline.Entries(t.Context())
	if _, ok := got["cellranger/9.0.1"]; !ok {
		t.Error("cached index not served")
	}
	if !offline[0].Stale {
		t.Error("Stale not set after serving an expired copy")
	}

	// With no cache at all, the source offers nothing and records why. Not an
	// error: the catalog is empty, which callers already handle, and the reason
	// lives on the source for one report to cover every source at once.
	bare, err := Open(t.Context(), []Spec{{Name: "remote", Base: srv.URL}}, Cache{})
	if err != nil {
		t.Fatal(err)
	}
	if got := bare.Entries(t.Context()); len(got) != 0 {
		t.Errorf("Entries = %v, want empty", got)
	}
	if bare[0].Err == nil {
		t.Error("an unreachable, uncached source should carry its error")
	}
}
