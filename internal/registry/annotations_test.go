package registry

import (
	"errors"
	"strings"
	"testing"
	"time"

	"github.com/Justype/condatainer/catalog"
	"github.com/Justype/condatainer/internal/artifact/meta"
)

func fullManifest() meta.Manifest {
	return meta.Manifest{
		SchemaVersion: meta.SchemaVersion,
		Name:          "grch38/star/2.7.11b/gencode49-101",
		Type:          catalog.TypeData,
		BuildType:     "script",
		Description:   "STAR index for GENCODE 49",
		URL:           "https://github.com/alexdobin/STAR",
		Platform:      meta.Platform{OS: "linux", Arch: "x86_64"},
		Keys: meta.Keys{
			Identity: meta.KeyRef{Scheme: "script-identity-v1", SHA256: strings.Repeat("a", 64)},
			Equiv:    meta.KeyRef{Scheme: "script-equiv-v1", SHA256: strings.Repeat("b", 64)},
		},
		Build: meta.Build{
			Created: time.Date(2026, 7, 21, 13, 45, 0, 0, time.UTC),
			Source:  "https://github.com/Justype/cnt-scripts",
		},
	}
}

func TestAnnotations(t *testing.T) {
	ann := Annotations(fullManifest(), "zstd")

	want := map[string]string{
		AnnTitle:          "grch38/star/2.7.11b/gencode49-101",
		AnnVersion:        "gencode49-101",
		AnnDescription:    "STAR index for GENCODE 49",
		AnnURL:            "https://github.com/alexdobin/STAR",
		AnnSource:         "https://github.com/Justype/cnt-scripts",
		AnnCreated:        "2026-07-21T13:45:00Z",
		AnnSchema:         "1",
		AnnCompression:    "zstd",
		AnnIdentityScheme: "script-identity-v1",
		AnnIdentitySHA:    strings.Repeat("a", 64),
		AnnEquivScheme:    "script-equiv-v1",
		AnnEquivSHA:       strings.Repeat("b", 64),
	}
	for key, wantValue := range want {
		if got := ann[key]; got != wantValue {
			t.Errorf("%s = %q, want %q", key, got, wantValue)
		}
	}
	// A native artifact carries no arch annotation: the index child's platform
	// descriptor is the one place architecture lives.
	if _, ok := ann[AnnNoarch]; ok {
		t.Error("a native artifact was marked noarch")
	}
	if len(ann) != len(want) {
		t.Errorf("annotations = %v\nwant exactly the %d keys above", ann, len(want))
	}
}

// Every optional value is omitted rather than published empty: an empty
// annotation is a claim that the artifact knows something and it is blank.
func TestAnnotationsOmitsWhatIsNotRecorded(t *testing.T) {
	m := meta.Manifest{
		SchemaVersion: meta.SchemaVersion,
		Name:          "myenv",
		Type:          catalog.TypeApp,
		Platform:      meta.Platform{OS: "linux", Arch: "x86_64"},
	}
	ann := Annotations(m, "")

	for _, key := range []string{
		AnnDescription, AnnURL, AnnSource, AnnCreated, AnnCompression,
		AnnIdentityScheme, AnnIdentitySHA, AnnEquivScheme, AnnEquivSHA, AnnNoarch,
	} {
		if _, ok := ann[key]; ok {
			t.Errorf("%s was published as %q with nothing to say", key, ann[key])
		}
	}
	// Title, version, and schema are always known.
	if ann[AnnTitle] != "myenv" || ann[AnnVersion] != "myenv" || ann[AnnSchema] != "1" {
		t.Errorf("annotations = %v", ann)
	}
}

// #ARCH:noarch is the one architecture fact OCI cannot express structurally.
func TestAnnotationsMarksNoarch(t *testing.T) {
	m := fullManifest()
	m.Platform.Arch = meta.ArchNone

	if ann := Annotations(m, "zstd"); !IsNoarch(ann) {
		t.Errorf("noarch artifact not marked: %v", ann)
	}
	if ann := Annotations(fullManifest(), "zstd"); IsNoarch(ann) {
		t.Error("a native artifact was reported as noarch")
	}
}

// The keys round-trip: what push writes is what a puller reads back.
func TestKeysRoundTrip(t *testing.T) {
	m := fullManifest()
	ann := Annotations(m, "zstd")

	if got := Identity(ann); got != m.Keys.Identity {
		t.Errorf("Identity = %+v, want %+v", got, m.Keys.Identity)
	}
	if got := Equiv(ann); got != m.Keys.Equiv {
		t.Errorf("Equiv = %+v, want %+v", got, m.Keys.Equiv)
	}
}

// Half a key addresses nothing, so it reads back as absent rather than as a key
// with a blank half that could compare equal to another blank one.
func TestHalfAKeyIsAbsent(t *testing.T) {
	for _, ann := range []map[string]string{
		{AnnIdentityScheme: "script-identity-v1"},
		{AnnIdentitySHA: strings.Repeat("a", 64)},
		{},
	} {
		if got := Identity(ann); !got.Empty() {
			t.Errorf("Identity(%v) = %+v, want empty", ann, got)
		}
	}
}

func TestCheck(t *testing.T) {
	m := fullManifest()
	ann := Annotations(m, "zstd")

	t.Run("the right artifact passes", func(t *testing.T) {
		if err := Check(ann, Want{Name: m.Name, Identity: m.Keys.Identity}); err != nil {
			t.Errorf("Check: %v", err)
		}
	})

	t.Run("a pull by name insists on no identity", func(t *testing.T) {
		if err := Check(ann, Want{Name: m.Name}); err != nil {
			t.Errorf("Check: %v", err)
		}
	})

	t.Run("a different name is a mismatch", func(t *testing.T) {
		err := Check(ann, Want{Name: "cellranger/9.0.1"})
		if !errors.Is(err, ErrMismatch) {
			t.Errorf("err = %v, want ErrMismatch", err)
		}
	})

	t.Run("a different identity is a mismatch", func(t *testing.T) {
		other := meta.KeyRef{Scheme: "script-identity-v1", SHA256: strings.Repeat("c", 64)}
		err := Check(ann, Want{Name: m.Name, Identity: other})
		if !errors.Is(err, ErrMismatch) {
			t.Errorf("err = %v, want ErrMismatch", err)
		}
	})

	t.Run("one digest under two schemes is two artifacts", func(t *testing.T) {
		other := meta.KeyRef{Scheme: "conda-explicit-v1", SHA256: m.Keys.Identity.SHA256}
		err := Check(ann, Want{Name: m.Name, Identity: other})
		if !errors.Is(err, ErrMismatch) {
			t.Errorf("err = %v, want ErrMismatch — the scheme is half the key", err)
		}
	})

	t.Run("a manifest with no annotations is not ours", func(t *testing.T) {
		for _, empty := range []map[string]string{nil, {}, {"org.opencontainers.image.vendor": "x"}} {
			if err := Check(empty, Want{}); !errors.Is(err, ErrNoAnnotations) {
				t.Errorf("Check(%v) = %v, want ErrNoAnnotations", empty, err)
			}
		}
	})
}

// Both directions are named. Telling someone to upgrade when their artifact is
// too *old* sends them the wrong way.
func TestCheckSchemaNamesTheDirection(t *testing.T) {
	newer := Annotations(fullManifest(), "")
	newer[AnnSchema] = "99"
	err := Check(newer, Want{})
	if !errors.Is(err, ErrIncompatible) {
		t.Fatalf("err = %v, want ErrIncompatible", err)
	}
	if !strings.Contains(err.Error(), "upgrade") {
		t.Errorf("a future artifact should say to upgrade: %v", err)
	}

	older := Annotations(fullManifest(), "")
	older[AnnSchema] = "0"
	err = Check(older, Want{})
	if !errors.Is(err, ErrIncompatible) {
		t.Fatalf("err = %v, want ErrIncompatible", err)
	}
	if strings.Contains(err.Error(), "upgrade") {
		t.Errorf("an artifact predating this format must not say to upgrade: %v", err)
	}

	missing := Annotations(fullManifest(), "")
	delete(missing, AnnSchema)
	if err := Check(missing, Want{}); !errors.Is(err, ErrNoAnnotations) {
		t.Errorf("err = %v, want ErrNoAnnotations", err)
	}
}

// #LICENSE: reaches the standard slot verbatim. Nothing normalizes or parses it:
// the publication decision is #REDISTRIBUTE:, made before the push.
func TestAnnotationsCarryLicenseVerbatim(t *testing.T) {
	m := fullManifest()
	m.License = "GPL-3.0-or-later AND LicenseRef-Vendor"
	ann := Annotations(m, "")
	if got := ann[AnnLicenses]; got != m.License {
		t.Errorf("%s = %q, want %q", AnnLicenses, got, m.License)
	}
	if AnnLicenses != "org.opencontainers.image.licenses" {
		t.Errorf("licences belong in the standard slot, got %q", AnnLicenses)
	}
}

// Redistribution is a decision, not metadata a consumer acts on before fetching
// a blob, so it must never appear as an annotation — the embedded manifest is
// where it lives.
func TestAnnotationsCarryNoRedistributionClaim(t *testing.T) {
	no := false
	m := fullManifest()
	m.Redistribute = &no
	for key, value := range Annotations(m, "") {
		if strings.Contains(strings.ToLower(key+value), "redistribut") {
			t.Errorf("annotation %s=%s leaks a redistribution claim", key, value)
		}
	}
}

func TestAnnotationsOmitAnAbsentLicense(t *testing.T) {
	if _, ok := Annotations(fullManifest(), "")[AnnLicenses]; ok {
		t.Error("a recipe that declared no #LICENSE: must publish no licences annotation")
	}
}
