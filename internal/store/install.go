package store

import (
	"errors"
	"fmt"
	"io"
	"os"
	"path/filepath"
	"slices"

	"github.com/Justype/condatainer/catalog"
	artifactcache "github.com/Justype/condatainer/internal/artifact/cache"
	"github.com/Justype/condatainer/internal/artifact/compare"
	"github.com/Justype/condatainer/internal/artifact/meta"
	"github.com/Justype/condatainer/internal/config"
	"github.com/Justype/condatainer/internal/image"
	"github.com/Justype/condatainer/internal/image/producer"
	"github.com/Justype/condatainer/internal/utils"
)

var (
	ErrIdentityCollision = errors.New("artifact identity filename collision")
	ErrTransactionClosed = errors.New("store transaction is closed")
)

// BeginOptions selects the destination and lookup roots. Empty ImagesDir uses
// config.GetWritableImagesDir; nil SearchDirs uses configured image roots.
type BeginOptions struct {
	ImagesDir  string
	SearchDirs []string
	Equiv      meta.KeyRef
	// StoreOnly installs into store/ even when the bare name is free, because
	// something else is entitled to that name. An exact copy already at the flat
	// name is still adopted: this decides where a new artifact is written, not
	// what answers to the name.
	StoreOnly bool
}

// Transaction owns one specific target producer lock and its prepared sibling.
type Transaction struct {
	Name     string
	Identity meta.KeyRef
	Equiv    meta.KeyRef
	// Layout is where the target landed: flat when the bare name was free,
	// stored when a different identity already held it.
	Layout Layout
	// Root is the images directory containing the target, for either layout.
	Root       string
	TargetPath string
	Prepared   string
	Adopted    *Candidate
	guard      *producer.Guard
	closed     bool
	read       artifactReader
}

// Begin adopts an existing exact artifact or reserves one collision-safe store
// pathname. It never locks the whole store directory.
func Begin(name string, identity meta.KeyRef, opts BeginOptions) (*Transaction, error) {
	return begin(name, identity, opts, compare.Read)
}

func begin(name string, identity meta.KeyRef, opts BeginOptions, read artifactReader) (*Transaction, error) {
	name = catalog.Normalize(name)
	if name == "" {
		return nil, errors.New("store artifact name is empty")
	}
	if _, err := ParseKeyRef(FormatKeyRef(identity)); err != nil {
		return nil, err
	}
	imagesDir := opts.ImagesDir
	var err error
	if imagesDir == "" {
		imagesDir, err = config.GetWritableImagesDir()
		if err != nil {
			return nil, err
		}
	} else if !utils.CanWriteToDir(imagesDir) {
		return nil, fmt.Errorf("images directory is not writable: %s", imagesDir)
	}
	dirs := opts.SearchDirs
	if dirs == nil {
		dirs = config.GetImageSearchPaths()
	}
	if !slices.Contains(dirs, imagesDir) {
		dirs = append([]string{imagesDir}, dirs...)
	}
	tx := &Transaction{Name: name, Identity: identity, Equiv: opts.Equiv, read: read}
	report := scan(ScanOptions{Dirs: dirs, Name: name}, read, nil)
	if candidate, err := resolveReport(name, IdentityQuery{Scheme: identity.Scheme, SHA256: identity.SHA256}, report); err == nil {
		if !opts.Equiv.Empty() && candidate.Equiv != opts.Equiv {
			return nil, fmt.Errorf("existing exact identity has equivalence %s, expected %s", FormatKeyRef(candidate.Equiv), FormatKeyRef(opts.Equiv))
		}
		tx.Adopted = &candidate
		tx.closed = true
		return tx, nil
	} else if !errors.Is(err, ErrNotFound) {
		return nil, err
	}

	// No exact copy anywhere, so this identity has to be written. Take the bare
	// name when nothing holds it: a flat artifact answers to `exec -o`, to
	// `list`, and to every other project, so the same identity is not rebuilt
	// once per checkout. The store is the conflict destination, not the default —
	// and where a caller that has yielded the name sends this identity instead.
	if !opts.StoreOnly {
		reserved, err := reserveFlat(tx, name, identity, opts.Equiv, imagesDir, dirs, read)
		if err != nil {
			return nil, err
		}
		if reserved {
			return tx, nil
		}
	}

	storeDir := filepath.Join(imagesDir, DirName)
	if err := utils.MkdirAllShared(storeDir); err != nil {
		return nil, fmt.Errorf("cannot create store directory: %w", err)
	}
	for chars := DefaultPrefixChars; chars <= len(identity.SHA256); chars++ {
		filename, err := Filename(name, identity, chars)
		if err != nil {
			return nil, err
		}
		target := filepath.Join(storeDir, filename)
		if info, err := os.Lstat(target); err == nil {
			if !info.Mode().IsRegular() {
				continue
			}
			artifact, readErr := read(target)
			if readErr != nil {
				continue
			}
			got := keyRef(artifact.IdentityScheme, artifact.Identity)
			if got == identity && artifact.Name == name {
				candidate := candidateFromArtifact(target, imagesDir, LayoutStored, info.Size(), artifact)
				tx.Adopted, tx.closed = &candidate, true
				return tx, nil
			}
			if got.SHA256 == identity.SHA256 && got.Scheme != identity.Scheme {
				return nil, fmt.Errorf("%w: %s uses %s and %s", ErrIdentityCollision, identity.Digest(), got.Scheme, identity.Scheme)
			}
			continue
		} else if !os.IsNotExist(err) {
			return nil, err
		}
		guard, err := producer.AcquireLocal(target)
		if err != nil {
			return nil, err
		}
		if candidate, occupied, err := inspectTarget(target, imagesDir, LayoutStored, name, identity, opts.Equiv, read); err != nil {
			guard.Release() //nolint:errcheck
			return nil, err
		} else if candidate != nil {
			guard.Release() //nolint:errcheck
			tx.Adopted, tx.closed = candidate, true
			return tx, nil
		} else if occupied {
			guard.Release() //nolint:errcheck
			continue
		}
		tx.Layout, tx.Root = LayoutStored, imagesDir
		tx.TargetPath = target
		tx.Prepared = producer.PreparedPath(target, guard.Info())
		tx.guard = guard
		return tx, nil
	}
	return nil, fmt.Errorf("%w: no filename remains for %s", ErrIdentityCollision, name)
}

// reserveFlat claims the bare-name target when no readable root holds that name,
// reporting whether it did. A false return means the name is occupied by another
// identity and the caller must overflow to the store.
//
// Occupancy is judged across every readable root, not just the destination:
// reads are nearest-first while writes go furthest-first, so a flat copy in a
// nearer root would shadow a new flat install and leave two identities
// answering to one name. Any file found here is a different identity by
// construction, because the exact lookup over the same roots has already missed.
func reserveFlat(tx *Transaction, name string, identity, equiv meta.KeyRef, imagesDir string, dirs []string, read artifactReader) (bool, error) {
	flatName := image.EncodeArtifactName(name) + ".sqf"
	for _, root := range dirs {
		if _, err := os.Lstat(filepath.Join(root, flatName)); err == nil {
			return false, nil
		} else if !os.IsNotExist(err) {
			return false, nil
		}
	}

	target := filepath.Join(imagesDir, flatName)
	guard, err := producer.AcquireLocal(target)
	if err != nil {
		return false, err
	}
	// Re-check under the lock: another writer may have taken the name between
	// the scan above and now, and the loser overflows rather than replacing it.
	candidate, occupied, err := inspectTarget(target, imagesDir, LayoutFlat, name, identity, equiv, read)
	if err != nil {
		guard.Release() //nolint:errcheck
		return false, err
	}
	if candidate != nil {
		guard.Release() //nolint:errcheck
		tx.Adopted, tx.closed = candidate, true
		return true, nil
	}
	if occupied {
		guard.Release() //nolint:errcheck
		return false, nil
	}
	tx.Layout, tx.Root = LayoutFlat, imagesDir
	tx.TargetPath = target
	tx.Prepared = producer.PreparedPath(target, guard.Info())
	tx.guard = guard
	return true, nil
}

func inspectTarget(path, root string, layout Layout, name string, identity, equiv meta.KeyRef, read artifactReader) (*Candidate, bool, error) {
	info, err := os.Lstat(path)
	if os.IsNotExist(err) {
		return nil, false, nil
	}
	if err != nil {
		return nil, false, err
	}
	if !info.Mode().IsRegular() || info.Mode()&os.ModeSymlink != 0 {
		return nil, true, nil
	}
	artifact, err := read(path)
	if err != nil {
		return nil, true, nil
	}
	got := keyRef(artifact.IdentityScheme, artifact.Identity)
	if artifact.Name == name && got == identity {
		candidate := candidateFromArtifact(path, root, layout, info.Size(), artifact)
		if _, err := ParseKeyRef(FormatKeyRef(candidate.Equiv)); err != nil {
			return nil, true, nil
		}
		if !equiv.Empty() && candidate.Equiv != equiv {
			return nil, true, fmt.Errorf("existing exact identity has equivalence %s, expected %s", FormatKeyRef(candidate.Equiv), FormatKeyRef(equiv))
		}
		return &candidate, true, nil
	}
	if got.SHA256 == identity.SHA256 && got.Scheme != identity.Scheme {
		return nil, true, fmt.Errorf("%w: %s uses %s and %s", ErrIdentityCollision, identity.Digest(), got.Scheme, identity.Scheme)
	}
	return nil, true, nil
}

// Commit verifies and atomically publishes the prepared SQF, or adopts an exact
// target that appeared while the producer lock was being acquired.
func (tx *Transaction) Commit() (Candidate, error) {
	if tx == nil || tx.closed {
		if tx != nil && tx.Adopted != nil {
			return *tx.Adopted, nil
		}
		return Candidate{}, ErrTransactionClosed
	}
	defer tx.finish()
	info, err := os.Lstat(tx.Prepared)
	if err != nil || !info.Mode().IsRegular() || info.Mode()&os.ModeSymlink != 0 || filepath.Ext(tx.Prepared) != ".part" {
		return Candidate{}, fmt.Errorf("invalid prepared store artifact %s", tx.Prepared)
	}
	file, err := os.OpenFile(tx.Prepared, os.O_RDWR, 0)
	if err != nil {
		return Candidate{}, err
	}
	if err := file.Sync(); err != nil {
		file.Close()
		return Candidate{}, err
	}
	if err := file.Close(); err != nil {
		return Candidate{}, err
	}
	artifact, err := tx.read(tx.Prepared)
	if err != nil {
		return Candidate{}, err
	}
	identity, equiv := keyRef(artifact.IdentityScheme, artifact.Identity), keyRef(artifact.EquivScheme, artifact.Equiv)
	if artifact.Name != tx.Name || identity != tx.Identity || (!tx.Equiv.Empty() && equiv != tx.Equiv) {
		return Candidate{}, fmt.Errorf("prepared artifact does not match requested name and keys")
	}
	if _, err := ParseKeyRef(FormatKeyRef(equiv)); err != nil {
		return Candidate{}, fmt.Errorf("prepared artifact has invalid equivalence key: %w", err)
	}
	if candidate, occupied, err := inspectTarget(tx.TargetPath, tx.Root, tx.Layout, tx.Name, tx.Identity, tx.Equiv, tx.read); err != nil {
		return Candidate{}, err
	} else if candidate != nil {
		return *candidate, nil
	} else if occupied {
		return Candidate{}, fmt.Errorf("target appeared during publication: %s", tx.TargetPath)
	}
	if err := os.Rename(tx.Prepared, tx.TargetPath); err != nil {
		return Candidate{}, fmt.Errorf("failed to publish store artifact: %w", err)
	}
	utils.ShareWithParentGroup(tx.TargetPath)
	installed, err := os.Lstat(tx.TargetPath)
	if err != nil {
		return Candidate{}, err
	}
	cacheArtifact(tx.TargetPath, tx.Root, tx.Layout, artifact, installed)
	return candidateFromArtifact(tx.TargetPath, tx.Root, tx.Layout, installed.Size(), artifact), nil
}

// Reserved reports that this transaction holds a target nothing has written
// yet. False means an exact copy was adopted and there is nothing to produce.
func (tx *Transaction) Reserved() bool {
	return tx != nil && !tx.closed && tx.guard != nil
}

// Detach hands the reserved target to a producer that outlives this process —
// a scheduler job — and closes the transaction without releasing the lock.
//
// The lock is what holds the pathname across the queue wait. The job adopts it
// by job ID, publishes through its own transaction, and releases it there. A job
// that never runs leaves a lock that goes stale with it, which the next producer
// clears.
func (tx *Transaction) Detach(info producer.Info) error {
	if !tx.Reserved() {
		return ErrTransactionClosed
	}
	if err := producer.Overwrite(producer.Path(tx.TargetPath), info); err != nil {
		return err
	}
	tx.Prepared = producer.PreparedPath(tx.TargetPath, info)
	tx.guard, tx.closed = nil, true
	return nil
}

// Abort removes this producer's prepared output and releases its target lock.
func (tx *Transaction) Abort() {
	if tx != nil {
		tx.finish()
	}
}

func (tx *Transaction) finish() {
	if tx.closed {
		return
	}
	os.Remove(tx.Prepared) //nolint:errcheck
	if tx.guard != nil {
		tx.guard.Release()
	} //nolint:errcheck
	tx.closed = true
}

// InstallFile copies a completed external SQF through the same transaction.
func InstallFile(name string, identity meta.KeyRef, source string, opts BeginOptions) (Candidate, error) {
	tx, err := Begin(name, identity, opts)
	if err != nil {
		return Candidate{}, err
	}
	if tx.Adopted != nil {
		return *tx.Adopted, nil
	}
	defer tx.Abort()
	in, err := os.Open(source)
	if err != nil {
		return Candidate{}, err
	}
	defer in.Close()
	out, err := utils.CreateFileWritable(tx.Prepared)
	if err != nil {
		return Candidate{}, err
	}
	if _, err = io.Copy(out, in); err != nil {
		out.Close()
		return Candidate{}, err
	}
	if err = out.Close(); err != nil {
		return Candidate{}, err
	}
	return tx.Commit()
}

func candidateFromArtifact(path, root string, layout Layout, size int64, artifact compare.Artifact) Candidate {
	return Candidate{Name: artifact.Name, Path: path, Root: root, Layout: layout, Size: size,
		Identity: keyRef(artifact.IdentityScheme, artifact.Identity), Equiv: keyRef(artifact.EquivScheme, artifact.Equiv)}
}

func cacheArtifact(path, root string, layout Layout, artifact compare.Artifact, info os.FileInfo) {
	identity, equiv := keyRef(artifact.IdentityScheme, artifact.Identity), keyRef(artifact.EquivScheme, artifact.Equiv)
	artifactcache.Default().Merge(path, info, func(record *artifactcache.Record) {
		record.Name, record.Identity, record.Equiv, record.KeysVerified = artifact.Name,
			artifactcache.Key{Scheme: identity.Scheme, SHA256: identity.SHA256}, artifactcache.Key{Scheme: equiv.Scheme, SHA256: equiv.SHA256}, true
		record.Layout, record.Root, record.Size = string(layout), root, info.Size()
	})
}
