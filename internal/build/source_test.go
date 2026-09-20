package build

import (
	"bytes"
	"compress/gzip"
	"context"
	"crypto/sha256"
	"encoding/hex"
	"errors"
	"log/slog"
	"net/http"
	"net/http/httptest"
	"os"
	"path/filepath"
	"strings"
	"testing"
	"time"

	"github.com/Justype/condatainer/catalog"
	"github.com/Justype/condatainer/internal/artifact/meta"
)

func sumOf(text string) string {
	sum := sha256.Sum256([]byte(text))
	return hex.EncodeToString(sum[:])
}

func TestFetchSourcesRecordsWhatArrived(t *testing.T) {
	const gtf, genome = "annotation bytes", "assembly bytes"
	srv := httptest.NewServer(http.HandlerFunc(func(w http.ResponseWriter, r *http.Request) {
		switch r.URL.Path {
		case "/gtf":
			w.Write([]byte(gtf)) //nolint:errcheck
		case "/genome":
			w.Write([]byte(genome)) //nolint:errcheck
		default:
			w.WriteHeader(http.StatusNotFound)
		}
	}))
	defer srv.Close()

	dir := t.TempDir()
	got, err := fetchSources(context.Background(), []catalog.SourceURL{
		{Name: "gtf", URL: srv.URL + "/gtf"},
		{Name: "genome", URL: srv.URL + "/genome"},
	}, filepath.Join(dir, "sources"))
	if err != nil {
		t.Fatal(err)
	}

	want := []meta.SourceFile{
		{Name: "gtf", SHA256: sumOf(gtf)},
		{Name: "genome", SHA256: sumOf(genome)},
	}
	if len(got) != len(want) {
		t.Fatalf("fetched %d sources, want %d", len(got), len(want))
	}
	for i := range want {
		if got[i] != want[i] {
			t.Errorf("source %d = %+v, want %+v", i, got[i], want[i])
		}
	}

	// The body reads the file, not the hash, so it has to actually be there.
	data, err := os.ReadFile(filepath.Join(dir, "sources", "gtf"))
	if err != nil {
		t.Fatal(err)
	}
	if string(data) != gtf {
		t.Errorf("file holds %q, want %q", data, gtf)
	}
}

// A truncated download must fail rather than be hashed: its digest would be a
// perfectly good digest for bytes nobody wants.
func TestShortReadFailsAndLeavesNothingBehind(t *testing.T) {
	srv := httptest.NewServer(http.HandlerFunc(func(w http.ResponseWriter, r *http.Request) {
		w.Header().Set("Content-Length", "100")
		w.Write([]byte("only twenty bytes...")) //nolint:errcheck
	}))
	defer srv.Close()

	dir := t.TempDir()
	_, err := fetchSources(context.Background(),
		[]catalog.SourceURL{{Name: "gtf", URL: srv.URL}}, dir)
	if err == nil {
		t.Fatal("a short read was accepted")
	}

	entries, readErr := os.ReadDir(dir)
	if readErr != nil {
		t.Fatal(readErr)
	}
	for _, e := range entries {
		t.Errorf("left %q behind; a partial download must be removed", e.Name())
	}
}

func TestSourceEnvNamesWhatArrived(t *testing.T) {
	got := sourceEnvFor([]meta.SourceFile{{Name: "gtf"}, {Name: "genome"}})
	want := []string{"CNT_SRC_gtf=/cnt_src/gtf", "CNT_SRC_genome=/cnt_src/genome"}
	if strings.Join(got, " ") != strings.Join(want, " ") {
		t.Errorf("env = %v, want %v", got, want)
	}
}

// A vendor download link carries an auth token in its query. Nothing that
// reaches a log may hold one.
func TestHostOfKeepsOnlyTheHost(t *testing.T) {
	for _, tc := range []struct{ in, want string }{
		{"https://cf.10xgenomics.com/rel/cellranger-9.0.1.tar.gz?Expires=1&Signature=SECRET", "cf.10xgenomics.com"},
		{"https://ftp.ebi.ac.uk/pub/gencode/v49.gtf.gz", "ftp.ebi.ac.uk"},
		{"example.invalid", "example.invalid"},
	} {
		if got := hostOf(tc.in); got != tc.want {
			t.Errorf("hostOf(%q) = %q, want %q", tc.in, got, tc.want)
		}
		if strings.Contains(hostOf(tc.in), "SECRET") {
			t.Errorf("hostOf(%q) leaked the query", tc.in)
		}
	}
}

// A compute node with no route out reaches the internet through the tunnel that
// sets http_proxy, so a transport that ignored the environment would fail exactly
// where the proxy exists to help. A bare http.Transport has no Proxy, which is
// what this guards.
//
// Not driven through a live proxy: ProxyFromEnvironment reads the environment
// once per process, so a test cannot set it after another has already fetched.
// A server that labels a .gz file Content-Encoding: gzip must not have it
// decoded: the recipe and the recorded digest are about the published bytes.
func TestFetchKeepsTheBytesTheServerSent(t *testing.T) {
	var packed bytes.Buffer
	zw := gzip.NewWriter(&packed)
	zw.Write([]byte("annotation bytes")) //nolint:errcheck
	zw.Close()
	srv := httptest.NewServer(http.HandlerFunc(func(w http.ResponseWriter, r *http.Request) {
		w.Header().Set("Content-Encoding", "gzip")
		w.Write(packed.Bytes()) //nolint:errcheck
	}))
	defer srv.Close()

	dir := t.TempDir()
	got, err := fetchSources(context.Background(), []catalog.SourceURL{{Name: "gtf", URL: srv.URL}}, dir)
	if err != nil {
		t.Fatal(err)
	}
	if want := sumOf(packed.String()); got[0].SHA256 != want {
		t.Errorf("digest = %s, want %s (the gzip bytes as sent)", got[0].SHA256, want)
	}
}

// A connection that goes quiet mid-body must fail rather than hang the build.
func TestStalledTransferFailsAndLeavesNothingBehind(t *testing.T) {
	old := sourceIdleTimeout
	sourceIdleTimeout = 100 * time.Millisecond
	defer func() { sourceIdleTimeout = old }()

	release := make(chan struct{})
	srv := httptest.NewServer(http.HandlerFunc(func(w http.ResponseWriter, r *http.Request) {
		w.Header().Set("Content-Length", "100")
		w.Write([]byte("partial")) //nolint:errcheck
		w.(http.Flusher).Flush()
		<-release
	}))
	defer srv.Close()
	defer close(release)

	dir := t.TempDir()
	_, err := fetchSources(context.Background(), []catalog.SourceURL{{Name: "gtf", URL: srv.URL}}, dir)
	if !errors.Is(err, errSourceStalled) {
		t.Fatalf("error = %v, want a stalled transfer", err)
	}
	if _, statErr := os.Stat(filepath.Join(dir, "gtf.part")); statErr == nil {
		t.Error("a partial file was left behind")
	}
}

func TestSourceTransportReadsTheProxyEnvironment(t *testing.T) {
	if newSourceTransport().Proxy == nil {
		t.Error("the transport has no Proxy, so http_proxy and https_proxy are ignored")
	}
}

func TestResolveAskedTakesAnswersInDeclarationOrder(t *testing.T) {
	got, err := resolveAsked([]catalog.SourceURL{
		{Name: "gtf", URL: "https://example.invalid/a"},
		{Name: "crx", Prompt: "first"},
		{Name: "idx", Prompt: "second"},
	}, []string{"https://vendor.invalid/one", "https://vendor.invalid/two"})
	if err != nil {
		t.Fatal(err)
	}
	if got[0].URL != "https://example.invalid/a" {
		t.Errorf("a declared source was altered: %+v", got[0])
	}
	if got[1].URL != "https://vendor.invalid/one" || got[2].URL != "https://vendor.invalid/two" {
		t.Errorf("answers mapped wrongly: %+v", got)
	}
}

// --yes supplies empty answers. That must be a refusal, not a fetch of "".
// Neither refusal may echo the answer: a pasted link is a per-user secret.
func TestResolveAskedRefusesWithoutEchoingTheAnswer(t *testing.T) {
	ask := []catalog.SourceURL{{Name: "crx", Prompt: "link"}}
	for _, tc := range []struct {
		name    string
		answers []string
	}{
		{"empty", []string{""}},
		{"none supplied", nil},
		{"not http", []string{"ftp://vendor.invalid/x?Signature=SECRET"}},
	} {
		t.Run(tc.name, func(t *testing.T) {
			_, err := resolveAsked(ask, tc.answers)
			if err == nil {
				t.Fatal("accepted an unusable answer")
			}
			if strings.Contains(err.Error(), "SECRET") {
				t.Errorf("error echoed the answer: %v", err)
			}
		})
	}
}

// The first #INPUT: answers reach the recipe and the rest go to the fetch.
func TestAnswersSplitBetweenRecipeAndFetch(t *testing.T) {
	b := &BuildObject{inputAnswers: []string{"y", "https://vendor.invalid/x"}, recipeInputs: 1}
	if got := b.recipeAnswers(); len(got) != 1 || got[0] != "y" {
		t.Errorf("recipe answers = %v, want the #INPUT: one only", got)
	}
	if got := b.sourceAnswers(); len(got) != 1 || got[0] != "https://vendor.invalid/x" {
		t.Errorf("source answers = %v, want the ask: one only", got)
	}
	none := &BuildObject{}
	if len(none.recipeAnswers()) != 0 || len(none.sourceAnswers()) != 0 {
		t.Error("a build with no prompts produced answers")
	}
}

// net/http puts the full request URL in its errors, and a pasted vendor link
// carries an auth token. A failed fetch must not print it.
func TestFetchErrorDoesNotLeakTheLink(t *testing.T) {
	srv := httptest.NewServer(http.HandlerFunc(func(w http.ResponseWriter, r *http.Request) {}))
	url := srv.URL + "/f?Signature=SECRET"
	srv.Close() // nothing listens now, so the dial fails

	_, err := fetchSources(context.Background(),
		[]catalog.SourceURL{{Name: "crx", URL: url}}, t.TempDir())
	if err == nil {
		t.Fatal("expected the fetch to fail")
	}
	if strings.Contains(err.Error(), "SECRET") {
		t.Errorf("error leaked the link: %v", err)
	}
}

func TestLinkMentionsVersionReadsSeparatorsAsOne(t *testing.T) {
	for _, tc := range []struct {
		link string
		want bool
	}{
		{"https://cf.example/cellranger-9.0.1.tar.gz", true},
		{"https://cf.example/cellranger_9_0_1.tar.gz", true},
		{"https://cf.example/cellranger-9-0-1.tar.gz", true},
		{"https://cf.example/cellranger-8.0.1.tar.gz", false},
		{"https://cf.example/download?token=abc", false},
	} {
		if got := linkMentionsVersion(tc.link, "9.0.1"); got != tc.want {
			t.Errorf("linkMentionsVersion(%q) = %v, want %v", tc.link, got, tc.want)
		}
	}
}

// A name whose last component is not a version has nothing to look for.
func TestVersionOfIgnoresWhatIsNotAVersion(t *testing.T) {
	for name, want := range map[string]string{
		"cellranger/9.0.1":                  "9.0.1",
		"grch38/star/2.7.11b/gencode49-101": "gencode49-101",
		"samtools":                          "",
		"grch38/genome/gencode":             "",
	} {
		if got := versionOf(name); got != want {
			t.Errorf("versionOf(%q) = %q, want %q", name, got, want)
		}
	}
}

// The warning names the artifact and version and never the link, which carries
// an auth token.
func TestVersionWarningDoesNotEchoTheLink(t *testing.T) {
	var out bytes.Buffer
	log := slog.New(slog.NewTextHandler(&out, nil))

	warnIfLinkLacksVersion(log, "cellranger/9.0.1", "https://cf.example/x?Signature=SECRET")
	if !strings.Contains(out.String(), "9.0.1") {
		t.Errorf("no warning for a link without the version: %q", out.String())
	}
	if strings.Contains(out.String(), "SECRET") {
		t.Errorf("the warning echoed the link: %q", out.String())
	}

	out.Reset()
	warnIfLinkLacksVersion(log, "cellranger/9.0.1", "https://cf.example/cellranger-9.0.1.tar.gz")
	if out.Len() != 0 {
		t.Errorf("warned about a link that names the version: %q", out.String())
	}
}
