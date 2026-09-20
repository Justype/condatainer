package build

import (
	"bytes"
	"context"
	"crypto/sha256"
	"encoding/hex"
	"errors"
	"fmt"
	"io"
	"log/slog"
	"net/http"
	"net/url"
	"os"
	"path/filepath"
	"strings"
	"time"

	"github.com/Justype/condatainer/catalog"
	"github.com/Justype/condatainer/internal/artifact/meta"
	"github.com/Justype/condatainer/internal/logging"
	"github.com/Justype/condatainer/internal/utils"
)

// SourcePath is where fetched #SOURCE: inputs appear inside the container. It is
// bound read-only: a source is an input recorded by its digest, so a recipe that
// edited one in place would leave the record describing bytes that are gone.
const SourcePath = "/cnt_src"

// sourceEnvPrefix is how a recipe names a fetched source: $CNT_SRC_<name>.
const sourceEnvPrefix = "CNT_SRC_"

// sourceHeaderTimeout bounds waiting for a server to answer, not the transfer.
// A genome FASTA legitimately takes an hour; a server that has said nothing in
// two minutes is not going to.
const sourceHeaderTimeout = 2 * time.Minute

// sourceIdleTimeout is how long a transfer may go without a byte before it
// fails, so a connection that died silently does not hang the build. Not a
// bound on the transfer: a slow, steady download never trips it. A variable so
// a test can shorten it.
var sourceIdleTimeout = 5 * time.Minute

var errSourceStalled = errors.New("no data received")

// declaredSources reads the #SOURCE: lines this build will fetch, with template
// placeholders already substituted.
//
// Parsed from the embedded recipe rather than carried from resolution, for the
// reason recipeKeyArtifact exists: one spelling of the build's inputs. A second
// copy is a second answer waiting to disagree.
func (b *BuildObject) declaredSources() ([]catalog.SourceURL, error) {
	file, ok := b.spec.Source.RecipeFile()
	if !ok || b.spec.Source.Script == nil {
		return nil, nil
	}
	recipe, err := catalog.ParseRecipe(meta.RecipeFileName, bytes.NewReader(file.Data))
	if err != nil {
		return nil, fmt.Errorf("cannot read the recipe's sources: %w", err)
	}
	if len(b.spec.Source.Placeholders) > 0 {
		recipe, err = catalog.Expand(recipe, b.spec.Source.Placeholders)
		if err != nil {
			return nil, fmt.Errorf("cannot apply placeholders to the recipe's sources: %w", err)
		}
	}
	return recipe.Sources, nil
}

// fetchDeclaredSources downloads what the recipe declared and records the
// digests on the spec, where the manifest and both key paths read them.
//
// Called before the identity is predicted, because source digests are part of
// it. Nothing is downloaded to answer a question a published artifact could have
// answered: the prebuilt check runs first and decides on equivalence, which no
// source enters.
func (b *BuildObject) fetchDeclaredSources(ctx context.Context) error {
	sources, err := b.declaredSources()
	if err != nil || len(sources) == 0 {
		return err
	}
	sources, err = resolveAsked(sources, b.sourceAnswers())
	if err != nil {
		return err
	}
	if err := ensureWorkspaceRoot(b); err != nil {
		return err
	}
	fetched, err := fetchSources(ctx, sources, b.sourceDir())
	if err != nil {
		return err
	}
	b.spec.Source.Fetched = fetched
	return nil
}

// resolveAsked gives each ask: source the link its answer supplied, taking
// answers in declaration order.
//
// An empty or non-http answer is refused, and the refusal does not echo it: an
// answered link is a per-user secret. --yes supplies empty answers, which is
// this refusal rather than a fetch of "".
func resolveAsked(sources []catalog.SourceURL, answers []string) ([]catalog.SourceURL, error) {
	out := make([]catalog.SourceURL, len(sources))
	next := 0
	for i, src := range sources {
		out[i] = src
		if src.Prompt == "" {
			continue
		}
		if next >= len(answers) {
			return nil, fmt.Errorf("source %s needs a download link, but none was supplied", src.Name)
		}
		link := strings.TrimSpace(answers[next])
		next++
		if !strings.HasPrefix(link, "http://") && !strings.HasPrefix(link, "https://") {
			return nil, fmt.Errorf("source %s needs an http or https download link", src.Name)
		}
		out[i].URL = link
	}
	return out, nil
}

// warnIfLinkLacksVersion flags a pasted link that never mentions the artifact's
// version, before a multi-gigabyte download rather than after it.
//
// A warning and never a refusal: vendor URL layouts change, some carry the
// version only inside a token, and nothing can tell a wrong link from an
// unfamiliar one. The link is not echoed, since it carries an auth token.
func warnIfLinkLacksVersion(log *slog.Logger, name, link string) {
	version := versionOf(name)
	if version == "" || linkMentionsVersion(link, version) {
		return
	}
	log.Warn(fmt.Sprintf("the link does not mention version %s of %s; check it is the right download", version, name))
}

// versionOf is the last component of a name, or "" when that is not a version:
// a name with one component, or a last component holding no digit.
func versionOf(name string) string {
	i := strings.LastIndex(name, "/")
	if i < 0 {
		return ""
	}
	version := name[i+1:]
	if !strings.ContainsAny(version, "0123456789") {
		return ""
	}
	return version
}

// linkMentionsVersion reports whether link contains version, reading . _ and -
// as the same separator: 9.0.1, 9_0_1 and 9-0-1 are one version.
func linkMentionsVersion(link, version string) bool {
	norm := strings.NewReplacer("_", ".", "-", ".")
	return strings.Contains(norm.Replace(strings.ToLower(link)), norm.Replace(strings.ToLower(version)))
}

// sourceDir is where fetched sources are staged on the host, before being bound
// read-only into the container at SourcePath.
func (b *BuildObject) sourceDir() string { return filepath.Join(b.ws.Root, "sources") }

// fetchSources downloads every declared source into dir and returns what
// arrived, in declaration order.
//
// Sources are fetched one at a time. There are rarely more than three, and
// concurrent transfers would put several progress readers on one terminal line.
func fetchSources(ctx context.Context, sources []catalog.SourceURL, dir string) ([]meta.SourceFile, error) {
	if len(sources) == 0 {
		return nil, nil
	}
	if err := utils.MkdirAllShared(dir); err != nil {
		return nil, fmt.Errorf("cannot create source directory: %w", err)
	}

	out := make([]meta.SourceFile, 0, len(sources))
	for _, src := range sources {
		sum, err := fetchSource(ctx, src, filepath.Join(dir, src.Name))
		if err != nil {
			return nil, err
		}
		out = append(out, meta.SourceFile{Name: src.Name, SHA256: sum})
	}
	return out, nil
}

// fetchSource downloads one source to destPath and returns the SHA-256 of what
// arrived.
func fetchSource(ctx context.Context, src catalog.SourceURL, destPath string) (string, error) {
	logging.FromContext(ctx).Info(
		fmt.Sprintf("Fetching source %s from %s", src.Name, hostOf(src.URL)), "kind", "note")
	sum, err := downloadHashed(ctx, src.URL, destPath, src.Name)
	if err != nil {
		return "", fmt.Errorf("cannot fetch source %s from %s: %w", src.Name, hostOf(src.URL), err)
	}
	logging.FromContext(ctx).Info("source ready", "source", src.Name, "sha256", sum[:12])
	return sum, nil
}

// newSourceTransport reads the proxy from the environment: a compute node with
// no route out reaches the internet through the tunnel that sets http_proxy and
// https_proxy, and a bare Transport ignores both. Compression is off so the
// bytes saved are the bytes the server sent, even for a .gz labelled
// Content-Encoding: gzip.
func newSourceTransport() *http.Transport {
	return &http.Transport{
		Proxy:                 http.ProxyFromEnvironment,
		DisableCompression:    true,
		ResponseHeaderTimeout: sourceHeaderTimeout,
	}
}

// withoutURL drops the request URL from a transport error. net/http embeds it,
// and an answered link carries an auth token that must reach no log or error;
// the host is already named by the caller.
func withoutURL(err error) error {
	var urlErr *url.Error
	if errors.As(err, &urlErr) {
		return urlErr.Err
	}
	return err
}

// downloadHashed streams url to destPath, hashing on the way past, and returns
// the hex SHA-256.
//
// A short read fails rather than being recorded. The digest of a truncated file
// is a perfectly good digest, so recording one would turn a dropped connection
// into a permanent identity for bytes nobody wants — worse than the build
// failing on a half-written input, which at least says so.
func downloadHashed(ctx context.Context, url, destPath, label string) (string, error) {
	ctx, cancel := context.WithCancelCause(ctx)
	defer cancel(nil)
	idle := time.AfterFunc(sourceIdleTimeout, func() { cancel(errSourceStalled) })
	defer idle.Stop()

	client := &http.Client{Transport: newSourceTransport()}
	req, err := http.NewRequestWithContext(ctx, http.MethodGet, url, nil)
	if err != nil {
		return "", err
	}
	resp, err := client.Do(req)
	if err != nil {
		return "", stalledOr(ctx, withoutURL(err))
	}
	defer resp.Body.Close()
	if resp.StatusCode != http.StatusOK {
		return "", fmt.Errorf("HTTP %d", resp.StatusCode)
	}

	tmpPath := destPath + ".part"
	file, err := utils.CreateFileWritable(tmpPath)
	if err != nil {
		return "", err
	}

	digest := sha256.New()
	body := &idleReader{r: resp.Body, timer: idle}
	written, err := io.Copy(io.MultiWriter(file, digest), newSourceProgress(ctx, label, resp.ContentLength, body))
	closeErr := file.Close()
	if err == nil {
		err = closeErr
	}
	err = stalledOr(ctx, err)
	if err == nil && resp.ContentLength >= 0 && written != resp.ContentLength {
		err = fmt.Errorf("received %d bytes, server said %d", written, resp.ContentLength)
	}
	if err != nil {
		os.Remove(tmpPath) //nolint:errcheck — the partial is being discarded anyway
		return "", err
	}

	if err := os.Rename(tmpPath, destPath); err != nil {
		os.Remove(tmpPath) //nolint:errcheck
		return "", err
	}
	utils.ShareWithParentGroup(destPath)
	return hex.EncodeToString(digest.Sum(nil)), nil
}

// idleReader restarts the idle timer whenever bytes arrive.
type idleReader struct {
	r     io.Reader
	timer *time.Timer
}

func (i *idleReader) Read(b []byte) (int, error) {
	n, err := i.r.Read(b)
	if n > 0 {
		i.timer.Reset(sourceIdleTimeout)
	}
	return n, err
}

// stalledOr names an idle timeout as such, and otherwise returns err.
func stalledOr(ctx context.Context, err error) error {
	if err != nil && errors.Is(context.Cause(ctx), errSourceStalled) {
		return fmt.Errorf("%w for %s", errSourceStalled, sourceIdleTimeout)
	}
	return err
}

// Report every sourceProgressStep bytes or sourceProgressInterval, whichever
// comes first, matching what a registry pull draws.
const (
	sourceProgressStep     = int64(64 << 20)
	sourceProgressInterval = 2 * time.Second
)

// sourceProgress counts one source's bytes onto the interactive progress line.
// No lock: sources are fetched one at a time and io.Copy reads on one goroutine.
type sourceProgress struct {
	log        *slog.Logger
	label      string
	total      int64
	reader     io.Reader
	done       int64
	nextReport int64
	lastReport time.Time
	lastDone   int64
	start      time.Time
}

func newSourceProgress(ctx context.Context, label string, total int64, r io.Reader) io.Reader {
	return &sourceProgress{
		log: logging.FromContext(ctx), label: label, total: total, reader: r,
		nextReport: sourceProgressStep, lastReport: time.Now(), start: time.Now(),
	}
}

func (p *sourceProgress) Read(b []byte) (int, error) {
	n, err := p.reader.Read(b)
	if n > 0 {
		p.done += int64(n)
		now := time.Now()
		if p.done >= p.nextReport || now.Sub(p.lastReport) >= sourceProgressInterval {
			p.report(false)
			p.nextReport, p.lastReport, p.lastDone = p.done+sourceProgressStep, now, p.done
		}
	}
	if err == io.EOF {
		p.report(true)
	}
	return n, err
}

// report draws one line: what has arrived, out of what the server advertised,
// and the speed since the previous line (the average, on the last one). A server
// that advertised no length reports no total.
func (p *sourceProgress) report(final bool) {
	if p.log == nil {
		return
	}
	args := []any{"kind", "progress", "done", utils.FormatSize(p.done)}
	if p.total >= 0 {
		args = append(args, "total", utils.FormatSize(p.total))
	}
	bytes, since := p.done-p.lastDone, p.lastReport
	if final {
		bytes, since = p.done, p.start
	}
	if elapsed := time.Since(since).Seconds(); elapsed > 0 {
		args = append(args, "speed", utils.FormatSize(int64(float64(bytes)/elapsed))+"/s")
	}
	args = append(args, "final", final, "last", final)
	p.log.Info("Downloading "+p.label, args...)
}

// sourceEnvFor renders the $CNT_SRC_<name> settings a recipe body reads, from
// what actually arrived rather than what was declared — a name reaches the
// recipe only once there is a file behind it.
func sourceEnvFor(fetched []meta.SourceFile) []string {
	out := make([]string, 0, len(fetched))
	for _, f := range fetched {
		out = append(out, sourceEnvPrefix+f.Name+"="+SourcePath+"/"+f.Name)
	}
	return out
}

// hostOf renders a URL for a log line: the host alone. A vendor download link
// carries an auth token in its query, and nothing that reaches a log, a manifest
// or a key may hold one.
func hostOf(rawURL string) string {
	rest := rawURL
	if _, after, ok := strings.Cut(rest, "://"); ok {
		rest = after
	}
	if host, _, ok := strings.Cut(rest, "/"); ok {
		return host
	}
	return rest
}
