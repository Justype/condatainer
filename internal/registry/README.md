# internal/registry

`registry` publishes and fetches immutable CondaTainer artifacts through an OCI
registry. Registry is who the package talks to; OCI is the wire format it speaks.
The implementation uses oras-go directly, so compute nodes do not need an `oras`
binary.

The package has two related entry points:

- `Resolve(ctx, ref)` normalizes an upstream container image reference and
  returns the current machine's platform digest for build provenance.
- `Publish`, `ResolveArtifact`, `Check`, and `Pull` distribute completed `.sqf`
  overlays and base images. Writable `.img` overlays have no immutable identity
  and are never distributed.

## Artifact contract

Overlay and base images have distinct manifest and layer media types. Pull checks
both before downloading. An artifact larger than one layer is pushed as ordered
chunks without writing a second local copy; pull verifies every OCI descriptor,
checks destination free space, reassembles in manifest order, regenerates the
embedded identity/equivalence keys, and atomically renames into place.

Native artifacts are children of a multi-platform OCI image index. A `noarch`
artifact is a plain manifest. `ResolveArtifact` always returns the selected image
manifest descriptor and its annotations, including its artifact type.

## References

- Versioned artifacts map the final name segment to the tag:
  `grch38/genome/gencode49` becomes repository `grch38/genome`, tag
  `gencode49`.
- Version-less base and two-segment OS artifacts use the complete name as the
  repository and publish a `YYYYMMDD` date tag plus `latest`.
- A digest selector (`repository@sha256:...`) addresses exact manifest bytes.

The first tag for a versioned artifact is immutable by default. `Publish` only
allows replacement with `Force`; adding a previously absent architecture is an
addition rather than replacement.

### Project placement

`PublishRequest.Placement` overrides the derivation above with a repository and
an explicit tag list. It exists for one caller, `project push`, which keeps a
whole project in one repository and carries the artifact name in the tag —
`ProjectTags` renders those, `ParseProjectTag` reads one back. The two namespaces
cannot be confused, and each form ParseProjectTag accepts carries its own proof
of that: an identity-qualified tag carries `__`, which the segment grammar cannot
produce, so the name it qualifies may be a single component; a plain tag has only
its shape, so it must hold the `--` a name/version encodes to, or be the one name
with no version to encode — `env`, which every frozen environment is called.
Anything else reports false rather than inventing a name from a version segment.

The name is *refused* rather than truncated when it does not fit a tag. A
truncated tag is a different artifact's address, and the name is what a puller
reads back out of it.

Why a project needs a second, identity-qualified tag is a retention argument
rather than an addressing one, and it belongs where the decision was made:
[`docs/deployment/distribution.md`](../../docs/deployment/distribution.md).

## What a public endpoint takes

`Audience.Accepts` is the whole rule, and `docs/manuals/condatainer.md`'s
*Publishing rules* is the table it implements. One ordering in it is load
bearing.

**`#REDISTRIBUTE:` outranks every type default**, and the defaults are only what
the author already asserted by choosing a `#TYPE:` — never a guess about
licensing. Nothing overrides a refusal, no flag and no force, because the person
pressing it is rarely the person who agreed to the vendor's terms and a public
push cannot be taken back.

**A Conda build and a frozen environment are exempt from the app default** for
one mechanical reason: neither embeds a recipe, so there is nowhere for
`#REDISTRIBUTE:` to be written that travels with the artifact, and applying the
default would be a permanent refusal wearing a default's clothes. What the
channels cannot answer — a private or vendor channel — is reported from
`Build.Channels` at push time and decided by a person.

The argument for singling a frozen environment out was that nobody vetted what a
person installed by hand. It does not survive the sentence above it. The reason a
refusal cannot be overridden is that *the pusher is rarely the party who agreed
to the terms* — and for a snapshot that inverts: whoever freezes an environment
is exactly the person who installed everything in it, and the only party who can
answer the question at all. Refusing there withheld the decision from the one
person qualified to make it, while a Conda build assembled from the same channels
published freely. An environment is a conda environment — that is what
`meta.EnvPrefix` is — so treating the two differently was a distinction the
artifacts do not carry.

Three checks must never be added here, each for the same reason: a tool that
guesses in the permissive direction cannot take it back.

- **Adjudicating a Conda build's channels.** The only ways to try are a config
  allowlist nobody can maintain honestly, or interpreting a few hundred licence
  strings. Report them and let a person decide.
- **Deriving permission from `#LICENSE:`.** It is an SPDX expression published
  verbatim as `org.opencontainers.image.licenses`, never parsed. Reading
  redistribution rights out of one is a judgement a tool gets wrong.
- **Letting config enumerate types.** `types: [app]` beside a public endpoint
  would erase the rule with no error, where a wrong `audience` is a claim someone
  had to write down and defend. `audience` states a fact about the registry and
  the permitted set is derived from it — a declaration, not enforcement: nothing
  verifies the registry is actually restricted, and it is unrelated to a GitHub
  package's own visibility setting, which CondaTainer never reads or changes.

A writable `.img` is refused everywhere, which is structural rather than policy:
it has no identity, so there is nothing to publish it *as*.

## Verification and errors

Annotations are derived only from the embedded manifest. `Check` rejects an
unsupported metadata schema or a requested name/identity mismatch before payload
transfer. Pull then verifies that the payload regenerates the published keys
before installation.

Callers classify outcomes with `errors.Is`: `ErrNotFound`,
`ErrUnsupportedPlatform`, and `ErrUnavailable` may permit a local build fallback;
`ErrUnauthorized`, `ErrRateLimited`, `ErrIncompatible`, `ErrIncompatibleRegistry`,
`ErrInvalidArtifact`, `ErrMismatch`, and `ErrNoAnnotations` must be surfaced.

`ErrRateLimited` is deliberately not `ErrUnavailable`: a registry asking for a
pause is not a registry that is down, and treating it as one would turn a
sixty-second wait into a local rebuild. It is matched on positive evidence — a
`429`, the OCI `TOOMANYREQUESTS` code, or GitHub's secondary-limit wording, which
arrives as an ordinary `403 DENIED` and is distinguishable only by its message.

## Retry

Every registry *write* — blob, manifest, tag, index — runs under
`retryPolicy.run`, which waits out a rate limit rather than failing. Four retries
at roughly 1, 2, 4, and 8 minutes, jittered, honouring a server-directed delay
when one survives the transport.

Only `ErrRateLimited` is retried. An ordinary `403` fails immediately, because
spending fifteen minutes of backoff to arrive at "your token is wrong" would be
worse than the bug this fixes. A `401` is retried once without waiting: ORAS
fixes the `Authorization` header when it opens an upload session and reuses it to
finalize the blob, so a token that ages out mid-layer is only recoverable by
starting the layer over.

Two details are load-bearing. Each attempt builds a **new** `io.SectionReader`,
since ORAS streams from a body with no `GetBody` and a reused reader would send
nothing — that is why retry lives in `pushRange` rather than in the HTTP client.
And after a wait the blob's presence is **rechecked**, because a registry can
commit and then lose the response, and without the check the next attempt
retransmits gigabytes that already landed.

A pause reports itself and its resume, carrying the same attribute the progress
line does so the three read as one event:

```text
[CNT] upload progress done=1.25 GB total=2.00 GB layer=4/11
[CNT] upload rate limited layer=4/11 retry=1/4 wait=1m0s
[CNT] upload resumed layer=4/11
```

The pause is logged at `Warn`, which is also what closes the in-place progress
line: `clilog` ends an open line on the first record that is not progress, so no
stalled byte count survives the wait. A retried layer's counter restarts at zero
— the progress reader is rebuilt per attempt — because the display is bytes
accepted by the current request, not traffic spent across retries.

Tests inject a policy through the context (`withRetryPolicy`) so nothing sleeps.

## Pull

`download` states its concurrency rather than inheriting ORAS's, which leaves the
field zero and fills in its own number at copy time. A pull runs once per node
per artifact, so across a cluster its request rate is the larger of the two
directions and worth deciding on purpose. Three is kept: a pull sits on the path
of every job that needs the artifact, and `retryingSource` is the real protection
— a rate limit pauses the pull instead of failing it, so trading throughput away
to avoid one would be paying twice.

Progress is reported for the **artifact**, not per blob:

```text
[CNT] download progress done=8.00 GB total=26.00 GB layers=6/13
```

Layers arrive concurrently, so per-blob lines would put two or three readers on
one terminal line overwriting each other, and a finished download would print one
"complete" line per layer, none of them naming the whole. `layers=` counts how
many have arrived, where push's `layer=` names the one in flight — push reports
per blob for the opposite reason, since one upload is in flight at a time and a
retry names it.

Retrying a fetch is safe in a way retrying a push is not: it is addressed by
digest and ORAS verifies what arrives, so a repeat either produces the same bytes
or fails. What is covered is *opening* the stream, which is where a `429` or a
secondary limit is answered. A connection lost midway through a two-gigabyte body
is not retried — ORAS offers no way to resume one, and restarting would mean
beginning the whole copy again.

This matters more than it looks: `ErrRateLimited` is deliberately not
`ErrUnavailable`, so a throttled pull may not fall back to a local build. Waiting
is the only thing that turns it back into a working pull.

`Fetch` is `Pull` without the install: it verifies and writes the payload to an
exact path and stops there. `Pull` owns the managed images directory, its
producer lock and its naming; a caller that has already decided where the bytes
go — `project restore`, staging into a store transaction's `.part` — must not
have that decided again underneath it. `Kind` says which media types to accept,
because a destination path settled by the caller no longer implies the extension
the type would have been derived from.

`SplitCoordinate` cuts a recorded `remotes` coordinate into the base and
repository the transport takes separately. It splits at the first slash after the
scheme, and deliberately does not trim a leading one: `/lab/p` names no registry,
and accepting it would send a fetch to whatever host happened to be configured.

## Provider profiles

`profileFor(host)` returns what a registry is known to enforce — per-layer
maximum, upload timeout, token lifetime, minimum gap between writes. Guardrails,
not protocol branching: nothing changes about *how* CondaTainer talks to a
registry, only how much it asks for at once and how fast.

Every field is zero when nothing is documented, and zero means "make no claim",
never "no limit applies". The table is deliberately short — GHCR is in it because
it published its numbers and then enforced an unpublished one on a real push.
An undocumented guess written down here becomes stale policy that outlives
whoever could correct it.

Two effects follow. The profile's `MaxLayerSize` clamps the layer plan, rounding
*down* to whole GiB so the result lands below the limit rather than on it —
GHCR's 10 GB becomes 9 GiB. And `MinMutationGap` spaces completed writes through
a `pacer` on the context: the clock runs from when the last write *finished*, so
a four-minute layer pays nothing, and a presence hit pays nothing because it
creates no content and never reaches the retry loop where pacing is applied.

## Throughput guard

Where a registry documents a per-upload limit, `throughputGuard` refuses to start
a layer that measurement says cannot finish inside it. The alternative is not
success: an upload timeout is not a rate limit, so nothing retries it, and the
push fails anyway having spent the ten minutes and the bandwidth.

The budget is 80% of the tighter of the profile's `UploadTimeout` and
`TokenLifetime` — both bound a layer for the same reason, since ORAS fixes the
credential when it opens the session. The estimate is the **slowest** of the last
few completed layers: conservative, so a layer is refused on what the link has
actually done at its worst, and recent, so a recovered link is not held back by
one bad sample. Only transferred layers are measured; one the registry already
had returns at once and would read as infinite bandwidth.

The first layer is never refused — the only honest estimate is a measured one —
and a registry that documents nothing gets no guard at all.

The refusal names the measured rate, the rate needed, and the projection. There
is no remedy to suggest beyond that: with no layer-size setting, a link too slow
for the size this artifact requires is a link too slow for this artifact.

## Transport

`inspectTransport` wraps the HTTP round tripper inside `newAuthClient`. It makes
no policy decisions; it preserves what ORAS discards and sends what ORAS omits.

- **Keeps the response metadata.** `errcode.ErrorResponse` carries method, URL,
  status, and the parsed error document — and no headers. `Retry-After` and the
  rate-limit counters live outside the body and would reach no caller otherwise.
  The captured failure lands in a slot the retry loop puts in the context, one
  per attempt, and annotates that attempt's error so a server-directed wait beats
  the local backoff.
- **Restores every body it reads.** A bounded prefix is read and put back in
  front of the rest, because ORAS decodes the document itself and a consumed body
  decodes to nothing — turning a specific registry error into "Forbidden".
- **Asks permission before a large body** (`Expect: 100-continue` on a `PUT` or
  `PATCH` past 4 MiB). ORAS streams straight into the request, so without it an
  auth failure, quota rejection, proxy body limit, or rate limit is discovered
  only after a whole layer has crossed the wire.
- **Abandons a dead upload session.** A failed blob finalize leaves a session
  ORAS will never resume; the retry loop `DELETE`s it before waiting, using the
  credential from the failed request, which is the last place that has it.
- **Identifies CondaTainer**, rather than sending ORAS's default `oras-go`.

## Layer planning

The limit that decides how an artifact is cut is a count of *requests*, not of
bytes: a fresh layer costs three, and a registry may stop accepting them after
some number it does not publish. `planLayerSize` therefore holds the layer
*count* near `targetLayers` and derives the size from the artifact, in whole GiB
steps, with a 2 GiB floor and clamped to any known per-layer maximum. A fixed
size would make the count grow with the artifact and turn every choice into a
supported-size cliff.

The size is a per-push decision recorded nowhere. Pull reads each boundary from
the manifest, so artifacts published at any earlier size stay readable and
nothing needs migrating.

**There is no setting and no flag.** The size comes from the constants in
`layerplan.go` and the artifact in hand, and nothing else — a publisher who has
to work out a layer size by hand has been failed by `planLayerSize`.

A fresh push asks the registry **once** whether a layer is already present.
Layers upload in order from the first, so a first layer the registry does not
have means none of this artifact is committed, and every further probe is a
request spent to be told 404 — a third of the budget, and the cheapest third to
get back. When the first layer *is* present, every layer is probed: that is a
resumed push or a second architecture, exactly the case the probe pays for.

Eliding a probe is never a correctness risk. A registry accepts a blob it already
holds and deduplicates by digest, and `errdef.ErrAlreadyExists` is tolerated, so
guessing wrong costs bandwidth in a case the heuristic has established is
unlikely.

`preflightUpload` settles the plan and checks what can be known without touching
the registry: the source is a non-empty regular file, the manifest fits the 4 MB
every registry is expected to accept, and no layer title, tag, or reference
overruns the lengths clients assume. It ends with one summary —

```text
[CNT] upload plan size=26.00 GB layer-size=2.00 GB (floor) layers=13 requests=~39 registry=ghcr.io
```

— and, when the estimate exceeds the ceiling a registry was observed to enforce,
one warning that the push will pause. That is a prediction, not a refusal: the
push is still correct and retry carries it through.

`probePushAccess` then opens a blob upload session and immediately abandons it.
Two requests, answering what nothing else does until the first layer has been
hashed and sent: whether the credential carries **push** scope — `checkTagIsFree`
only does a pull-scope `Resolve` — whether the repository exists, since ECR does
not create one on push, and whether its name is a shape this registry accepts.

Only a refusal the registry was clear about stops the push. Unreachable,
throttled, or an unexpected status is left to the push itself, which has retry
and better messages; a preflight that invents failure modes is worse than none.
The whole probe is bounded by `probeTimeout`, because ORAS's transport would
otherwise retry a 5xx five times and a cheap check would stop being cheap.

The manifest ceiling is also the only lower bound on layer size, and the reason
no arbitrary one is imposed: it is derived from a spec limit and from the
artifact being pushed, rather than guessed.

## Authentication

Credential precedence is:

1. `GITHUB_TOKEN` for `ghcr.io` only;
2. Docker/OCI credential store, including configured credential helpers;
3. anonymous access.

`Login` and `Logout` manage the same store, and `StoredCredentials` lists its
hosts. Tokens are never logged or listed.

## Callers

The publisher/admin surface is `condatainer registry push|pull|tags|resolve|
login|logout|list`, in [`cmd/registry.go`](../../cmd/registry.go). Choosing an
endpoint from flags, config, or an artifact's recorded `build.source` happens
there and never here: this package is told where to go. See the
[manual](../../docs/manuals/condatainer.md) for that surface.

`build` reaches the same transport without the CLI, trying the selected recipe
source's ordered `oci.pull` endpoints before building locally, and requiring the
candidate to match the equivalence key derived from that recipe. It and
`registry pull` take the same destination producer lock.

`project` is the third caller and the only one that is told a placement rather
than deriving one. `project pin` resolves against a collection's
endpoints to record what it already serves, `project push` publishes with a
`Placement`, and `project restore` fetches recorded digests through `Fetch`. What
may be published is decided by `Audience.Accepts` here, from the embedded
manifest alone — never from a filename, a flag, or which command is running.
