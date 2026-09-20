package freeze

import (
	"bytes"
	"context"
	"errors"
	"fmt"
	"os"
	"regexp"
	"runtime"
	"strings"

	"github.com/Justype/condatainer/internal/artifact/key"
	"github.com/Justype/condatainer/internal/artifact/meta"
	execpkg "github.com/Justype/condatainer/internal/runtime/exec"
	"github.com/Justype/condatainer/internal/utils"
)

// TreeIdentity hashes the payload of a packed artifact into its payload key.
//
// It reads the finished .sqf rather than the overlay it came from, so the key
// describes what the archive actually contains, and reading a .sqf through
// squashfuse costs a fraction of reading the .img through fuse2fs. The walk is
// this binary's own `_payload_key`, run inside the mount's namespace: the same
// code a build keys its payload with, and one that needs no host find or
// sha256sum.
func TreeIdentity(ctx context.Context, artifact string) (meta.KeyRef, error) {
	squashfuse, err := FindSquashfuse()
	if err != nil {
		return meta.KeyRef{}, err
	}
	self, err := executablePath()
	if err != nil {
		return meta.KeyRef{}, fmt.Errorf("locating condatainer binary: %w", err)
	}

	scratch := utils.GetTmpDir()
	if err := os.MkdirAll(scratch, 0o755); err != nil {
		return meta.KeyRef{}, fmt.Errorf("stage identity mount: %w", err)
	}
	mnt, err := os.MkdirTemp(scratch, "cnt-ident-")
	if err != nil {
		return meta.KeyRef{}, fmt.Errorf("stage identity mount: %w", err)
	}
	defer os.RemoveAll(mnt)

	work := fmt.Sprintf("%s _payload_key %s '' %d", shellQuote(self), shellQuote(mnt), runtime.GOMAXPROCS(0))
	var stdout, stderr bytes.Buffer
	err = MountedRun(ctx, squashfuse, []string{artifact}, mnt, work,
		execpkg.IO{Stdout: &stdout, Stderr: &stderr})
	if err != nil {
		return meta.KeyRef{}, fmt.Errorf("read %s to identify it: %w: %s", artifact, err, stderr.String())
	}

	// The mount may print notices around the command's own output; the key is
	// the one line that is a digest.
	for _, line := range strings.Fields(stdout.String()) {
		if sha256Hex.MatchString(line) {
			return meta.KeyRef{Scheme: string(key.PayloadTreeV1), SHA256: line}, nil
		}
	}
	return meta.KeyRef{}, fmt.Errorf("read %s to identify it: the payload walk printed no key", artifact)
}

var sha256Hex = regexp.MustCompile(`^[0-9a-f]{64}$`)

// ErrNoPayloadKey reports an artifact whose manifest records no payload key: a
// Conda or definition build, which its other keys already pin.
var ErrNoPayloadKey = errors.New("artifact records no payload key")

// VerifyPayload reads the packed artifact's payload and checks it against the key
// its manifest records: the payload key, or for a frozen environment its
// identity. An artifact with neither is ErrNoPayloadKey.
func VerifyPayload(ctx context.Context, artifact string) error {
	m, err := meta.ReadManifest(artifact)
	if err != nil {
		return err
	}
	want := m.Keys.Payload
	if m.BuildType == meta.BuildTypeSnapshot {
		want = m.Keys.Identity
	}
	if want.Empty() {
		return fmt.Errorf("%w: %s", ErrNoPayloadKey, artifact)
	}
	got, err := TreeIdentity(ctx, artifact)
	if err != nil {
		return err
	}
	if got.SHA256 != want.SHA256 {
		return fmt.Errorf("payload of %s is %s, its manifest records %s", artifact, got.Digest(), want.Digest())
	}
	return nil
}
