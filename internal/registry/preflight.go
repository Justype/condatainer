package registry

import (
	"context"
	"encoding/json"
	"errors"
	"fmt"
	"io"
	"net/http"
	"os"
	"path/filepath"
	"strings"
	"time"

	"oras.land/oras-go/v2/registry/remote"
	"oras.land/oras-go/v2/registry/remote/auth"
	"oras.land/oras-go/v2/registry/remote/errcode"

	"github.com/Justype/condatainer/internal/logging"
	"github.com/Justype/condatainer/internal/utils"
)

// What a push is checked against before it spends an hour finding out.
const (
	// requestsPerLayer is what one fresh layer costs once the presence probe has
	// been elided: POST to open the session, PUT to commit. The first layer costs
	// one more, for the probe that decides whether to elide the rest.
	requestsPerLayer = 2
	// fixedRequests is what surrounds the layers: the token exchange, the tag
	// probe, this preflight's own session and its cancel, the empty config blob,
	// the manifest, and the tags.
	fixedRequests = 12
	// expectedRequestCeiling is where a registry was observed to stop accepting
	// writes: GHCR refused a 41-layer push on chunk 34, about its 101st request.
	// Not exact — one measurement of a number nobody publishes.
	expectedRequestCeiling = 100

	// maxManifestBytes is the interoperability floor every registry accepts. At
	// descriptorBytes each it allows some 16,000 layers, so it will not fire in
	// practice; it is checked rather than assumed.
	maxManifestBytes = 4 << 20
	descriptorBytes  = 250

	// maxTagLength and maxReferenceLength are the OCI tag limit and the length
	// many clients assume a full registry/repository:tag reference stays within.
	maxTagLength       = 128
	maxReferenceLength = 255
)

// probeTimeout bounds the whole push preflight, including ORAS's own retries on
// a 5xx. Short because the probe's answer is optional: anything inconclusive is
// left to the push.
const probeTimeout = 5 * time.Second

// uploadPlan is what a push intends to do, settled before any bytes move.
type uploadPlan struct {
	Size      int64
	LayerSize int64
	Reason    layerSizeReason
	Layers    int
	// Requests is the estimated cost, reported rather than enforced: exceeding it
	// predicts a pause that retry carries through, not an error.
	Requests int
}

// preflightUpload settles the layer plan for artifactPath and checks everything
// knowable without touching the registry — each check being one that would
// otherwise surface after gigabytes have been digested and sent.
//
// The profile's per-layer maximum clamps the derived size here, so a registry's
// hard limit is respected before anything is hashed.
func preflightUpload(ctx context.Context, artifactPath, base string, tags []string, profile transferProfile) (uploadPlan, error) {
	info, err := os.Stat(artifactPath)
	if err != nil {
		return uploadPlan{}, err
	}
	if !info.Mode().IsRegular() {
		return uploadPlan{}, fmt.Errorf("%s is not a regular file, so its bytes are not a stable artifact", artifactPath)
	}
	if info.Size() == 0 {
		return uploadPlan{}, fmt.Errorf("%s is empty", artifactPath)
	}

	size := info.Size()
	layerSize, reason := planLayerSize(size, profile.MaxLayerSize)
	plan := uploadPlan{
		Size:      size,
		LayerSize: layerSize,
		Reason:    reason,
		Layers:    layerCount(size, layerSize),
	}
	plan.Requests = estimateRequests(plan.Layers)

	if err := checkManifestFits(plan.Layers); err != nil {
		return uploadPlan{}, err
	}
	if err := checkNames(base, artifactPath, tags, plan.Layers); err != nil {
		return uploadPlan{}, err
	}

	log := logging.FromContext(ctx)
	log.Info("upload plan",
		"size", utils.FormatBytes(plan.Size),
		"layer-size", fmt.Sprintf("%s (%s)", utils.FormatBytes(plan.LayerSize), plan.Reason),
		"layers", plan.Layers,
		"requests", "~"+fmt.Sprint(plan.Requests),
		"registry", registryHost(base))
	// An expectation, not a failure: the push is correct, it will simply pause,
	// and an operator told so in advance does not report it as a hang.
	if plan.Requests > expectedRequestCeiling {
		log.Warn("this push is larger than registries have been observed to accept in one run; expect it to pause and resume",
			"requests", plan.Requests, "observed-ceiling", expectedRequestCeiling)
	}
	return plan, nil
}

// probePushAccess opens a blob upload session and immediately abandons it.
//
// Two requests, answering what otherwise stays unanswered until the first layer
// has been digested and sent: whether the credential carries push scope, whether
// the repository exists — ECR does not create one on push — and whether its name
// is a shape this registry accepts. [checkTagIsFree] covers none of them; its
// Resolve is a pull-scope read.
//
// Only a refusal the registry was clear about stops the push. Anything else is
// left to the push, which has retry and better messages.
func probePushAccess(ctx context.Context, repository *remote.Repository) error {
	scheme := "https"
	if repository.PlainHTTP {
		scheme = "http"
	}
	ref := repository.Reference
	url := fmt.Sprintf("%s://%s/v2/%s/blobs/uploads/", scheme, ref.Registry, ref.Repository)

	// The same scope a real blob push asks for, or the token comes back without
	// push access and the probe fails for a reason of its own making.
	ctx = auth.AppendRepositoryScope(ctx, ref, auth.ActionPull, auth.ActionPush)
	ctx, stop := context.WithTimeout(ctx, probeTimeout)
	defer stop()
	req, err := http.NewRequestWithContext(ctx, http.MethodPost, url, nil)
	if err != nil {
		return err
	}
	resp, err := repository.Client.Do(req)
	if err != nil {
		logging.FromContext(ctx).Debug("push preflight could not reach the registry", "error", err)
		return nil
	}
	defer resp.Body.Close() //nolint:errcheck

	if resp.StatusCode == http.StatusAccepted {
		cancelProbeSession(ctx, repository, resp)
		return nil
	}

	err = classify(errutilParse(resp))
	switch {
	case errors.Is(err, ErrUnauthorized):
		return fmt.Errorf("%w: %s does not accept a push from this credential", err, FullRef(ref.Registry, ref.Repository, ""))
	case errors.Is(err, ErrNotFound):
		return fmt.Errorf("%w: %s/%s does not exist; some registries require the repository to be created before a push",
			err, ref.Registry, ref.Repository)
	case errors.Is(err, ErrIncompatibleRegistry):
		return err
	}
	logging.FromContext(ctx).Debug("push preflight was inconclusive",
		"status", resp.StatusCode, "repository", ref.Repository)
	return nil
}

// cancelProbeSession abandons the session the probe opened. Best effort:
// failing to tidy is not a reason to refuse a push that is otherwise fine.
func cancelProbeSession(ctx context.Context, repository *remote.Repository, opened *http.Response) {
	location, err := opened.Location()
	if err != nil {
		return
	}
	req, err := http.NewRequestWithContext(ctx, http.MethodDelete, location.String(), nil)
	if err != nil {
		return
	}
	resp, err := repository.Client.Do(req)
	if err != nil {
		return
	}
	resp.Body.Close() //nolint:errcheck
}

// errutilParse turns a failing response into the error ORAS would have built, so
// classify sees the same shape here as everywhere else.
func errutilParse(resp *http.Response) error {
	var body struct {
		Errors errcode.Errors `json:"errors"`
	}
	// Best effort: a registry answering with something other than an OCI error
	// document still has a status code, which is the part that decides.
	_ = json.NewDecoder(io.LimitReader(resp.Body, maxCapturedBody)).Decode(&body)
	return &errcode.ErrorResponse{
		Method:     resp.Request.Method,
		URL:        resp.Request.URL,
		StatusCode: resp.StatusCode,
		Errors:     body.Errors,
	}
}

// estimateRequests reports what a push of this many layers is expected to cost:
// one presence probe, then two requests per layer. A resumed push costs more —
// it probes every layer to find what landed — but those requests buy back
// gigabytes of upload.
func estimateRequests(layers int) int {
	return requestsPerLayer*layers + 1 + fixedRequests
}

// checkManifestFits refuses a layer count whose manifest could exceed the 4 MB
// every registry is expected to accept.
func checkManifestFits(layers int) error {
	if estimate := int64(layers) * descriptorBytes; estimate > maxManifestBytes {
		return fmt.Errorf("%d layers need about %s of manifest, past the %s registries are expected to accept",
			layers, utils.FormatBytes(estimate), utils.FormatBytes(maxManifestBytes))
	}
	return nil
}

// checkNames refuses a reference or layer title that would be rejected, or
// silently truncated, somewhere downstream.
//
// The longest title is checked rather than the first: a chunk suffix adds to
// every name equally, so the last layer is the one that overruns.
func checkNames(base, artifactPath string, tags []string, layers int) error {
	name := filepath.Base(artifactPath)
	if layers > 1 {
		name += fmt.Sprintf(chunkSuffix, layers-1)
	}
	if len(name) > maxReferenceLength {
		return fmt.Errorf("layer title %q is %d characters, past the %d clients assume", name, len(name), maxReferenceLength)
	}
	for _, tag := range tags {
		if len(tag) > maxTagLength {
			return fmt.Errorf("tag %q is %d characters, past the OCI limit of %d", tag, len(tag), maxTagLength)
		}
		if ref := TrimBaseScheme(base) + ":" + tag; len(ref) > maxReferenceLength {
			return fmt.Errorf("reference %q is %d characters, past the %d clients assume", ref, len(ref), maxReferenceLength)
		}
	}
	return nil
}

// registryHost returns just the host of a registry base, for the plan summary.
func registryHost(base string) string {
	host, _, _ := strings.Cut(TrimBaseScheme(base), "/")
	return host
}
