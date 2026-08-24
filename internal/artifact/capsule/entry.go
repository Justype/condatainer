package capsule

import (
	"encoding/json"
	"fmt"
	"os"
	"path/filepath"

	"github.com/Justype/condatainer/internal/artifact/key"
	"github.com/Justype/condatainer/internal/artifact/meta"
)

// An entry is one directory holding manifest.json plus exactly the rebuild
// sources that manifest names. The same directory appears in two places — inside
// an image at /.cnt/provenance, and in a checkout at cnt-lock/provenance — and
// it has to mean the same thing in both, because a lock is verified against
// records an image produced. Reading one therefore lives here rather than once
// per consumer: two readers that drift give a lock that verifies and an image
// that does not, or the reverse.

// Reader reads one file of an entry directory, given its full path. Callers
// supply their own so the policy for untrusted bytes stays theirs: a checkout's
// entries are tracked Git content and are read with a size bound and a symlink
// refusal, while an image's come straight off an extraction.
type Reader func(path string) ([]byte, error)

// Record is one entry read and verified end to end: the manifest, the source
// bytes beside it, and the keys those bytes regenerate.
type Record struct {
	Manifest meta.Manifest
	Sources  key.Sources
	Derived  key.Derived
}

// FileNames is everything an entry directory holds: the manifest, then the
// sources the manifest names. It is the single answer to "what is in an entry",
// used both to write one and to check one is whole.
func FileNames(m meta.Manifest) []string {
	return append([]string{meta.FileName}, m.Source.Files...)
}

// ReadRecord loads and fully verifies one entry directory. The manifest must
// parse and validate, every source it names must be present, and both keys must
// regenerate from those bytes to the digests recorded beside them — nothing is
// taken on trust, including the manifest's own keys.
func ReadRecord(dir string, read Reader) (Record, error) {
	manifest, err := ReadManifest(dir, read)
	if err != nil {
		return Record{}, err
	}
	sources, err := ReadSources(dir, manifest, read)
	if err != nil {
		return Record{}, err
	}
	derived, err := key.Verify(manifest, sources)
	if err != nil {
		return Record{}, fmt.Errorf("sources do not regenerate the recorded keys: %w", err)
	}
	return Record{Manifest: manifest, Sources: sources, Derived: derived}, nil
}

// ReadManifest parses and validates an entry's manifest.json.
func ReadManifest(dir string, read Reader) (meta.Manifest, error) {
	data, err := read(filepath.Join(dir, meta.FileName))
	if err != nil {
		return meta.Manifest{}, fmt.Errorf("cannot read %s: %w", meta.FileName, err)
	}
	var m meta.Manifest
	if err := json.Unmarshal(data, &m); err != nil {
		return meta.Manifest{}, fmt.Errorf("cannot parse %s: %w", meta.FileName, err)
	}
	m.Normalize()
	if err := meta.ValidateManifest(m); err != nil {
		return meta.Manifest{}, err
	}
	return m, nil
}

// ReadSources reads exactly the files the manifest names, keyed as it names
// them. A recipe build must carry its recipe byte for byte and a Conda build
// both of its exports; a URL or a catalog reference is never a substitute.
func ReadSources(dir string, m meta.Manifest, read Reader) (key.Sources, error) {
	sources := make(key.Sources, len(m.Source.Files))
	for _, name := range m.Source.Files {
		if name != filepath.Base(name) || name == "." || name == ".." {
			return nil, fmt.Errorf("source file %q is not a plain name", name)
		}
		data, err := read(filepath.Join(dir, name))
		if err != nil {
			return nil, fmt.Errorf("cannot read source %s: %w", name, err)
		}
		sources[name] = data
	}
	return sources, nil
}

// CheckDirName reports whether an entry directory is named for the record it
// holds. The name is addressing, never identity: the digest it carries is
// truncated, so the comparison is against the key regeneration produced.
func CheckDirName(dirName string, r Record) error {
	if want := EntryName(r.Manifest.Name, r.Derived.Identity.Ref.Digest()); want != dirName {
		return fmt.Errorf("directory name does not match its manifest; expected %s", want)
	}
	return nil
}

// PlainReader reads an entry file with no size bound, for entries already
// extracted from an image to local scratch. A checkout supplies its own reader
// instead, because tracked Git content is untrusted input.
var PlainReader Reader = os.ReadFile
