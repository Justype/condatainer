package meta

import (
	"encoding/json"
	"errors"
	"fmt"
	"log/slog"
	"os"
	"path/filepath"

	"github.com/Justype/condatainer/catalog"
	"github.com/Justype/condatainer/internal/image/sif"
	"github.com/Justype/condatainer/internal/image/squashfs"
	"github.com/Justype/condatainer/internal/image/tool"
	"github.com/Justype/condatainer/internal/utils"
)

// Stage writes the manifest into dir, ready to be packed into an image at
// meta.Path. The write goes through a temp file and an atomic rename, and the
// JSON is byte-identical for the same input.
func Stage(dir string, manifest Manifest) error {
	if err := utils.MkdirAllShared(dir); err != nil {
		return fmt.Errorf("failed to create metadata dir %s: %w", dir, err)
	}

	data, err := Marshal(manifest)
	if err != nil {
		return err
	}

	final := filepath.Join(dir, FileName)
	tmp := final + ".tmp"
	f, err := utils.CreateFileWritable(tmp)
	if err != nil {
		return fmt.Errorf("failed to create %s: %w", tmp, err)
	}
	defer os.Remove(tmp) // no-op after a successful rename

	if _, err := f.Write(data); err != nil {
		f.Close() //nolint:errcheck
		return fmt.Errorf("failed to write %s: %w", tmp, err)
	}
	if err := f.Close(); err != nil {
		return fmt.Errorf("failed to close %s: %w", tmp, err)
	}
	if err := os.Rename(tmp, final); err != nil {
		return fmt.Errorf("failed to install %s: %w", final, err)
	}
	utils.ShareWithParentGroup(final)
	return nil
}

// Marshal renders a manifest the way Stage writes it.
func Marshal(manifest Manifest) ([]byte, error) {
	data, err := json.MarshalIndent(manifest, "", "  ")
	if err != nil {
		return nil, fmt.Errorf("failed to encode manifest: %w", err)
	}
	return append(data, '\n'), nil
}

// Read returns the manifest embedded in an image. Only a genuinely absent
// manifest is ErrNoManifest; a host failure keeps its own cause. Cached by path,
// size and mtime, negative verdicts included. See the README's Manifests.
func Read(imagePath string) (Manifest, error) {
	abs, err := filepath.Abs(imagePath)
	if err != nil {
		abs = imagePath
	}
	fi, err := os.Stat(abs)
	if err != nil {
		globalCache.forget(abs)
		return Manifest{}, fmt.Errorf("%w: %s: %w", tool.ErrUnreadable, imagePath, err)
	}

	if manifest, cached, ok := globalCache.lookup(abs, fi); ok {
		if !cached {
			return Manifest{}, fmt.Errorf("%w: %s", ErrNoManifest, imagePath)
		}
		return manifest, nil
	}

	manifest, err := readUncached(abs)
	switch {
	case err == nil:
		globalCache.store(abs, fi, &manifest)
	case errors.Is(err, ErrNoManifest):
		globalCache.store(abs, fi, nil)
	}
	return manifest, err
}

// readUncached does the archive read and decode, bypassing the cache.
func readUncached(imagePath string) (Manifest, error) {
	data, err := readRaw(imagePath)
	if err != nil {
		if errors.Is(err, tool.ErrFileNotFound) {
			return Manifest{}, fmt.Errorf("%w: %s", ErrNoManifest, imagePath)
		}
		return Manifest{}, err
	}

	var manifest Manifest
	if err := json.Unmarshal(data, &manifest); err != nil {
		return Manifest{}, fmt.Errorf("%w: %s: %w", ErrInvalid, imagePath, err)
	}
	manifest.Normalize()
	if err := Validate(manifest); err != nil {
		return Manifest{}, fmt.Errorf("%s: %w", imagePath, err)
	}
	return manifest, nil
}

// CheckBase reports whether an image may serve as a container root. A manifest
// that reads has to declare type base; a missing one passes, an unreadable one
// warns and passes. See the README's Manifests.
func CheckBase(imagePath string) error {
	manifest, err := Read(imagePath)
	switch {
	case errors.Is(err, ErrNoManifest):
		return nil
	case err != nil:
		slog.Default().Warn("could not read base image metadata", "path", imagePath, "err", err)
		return nil
	case manifest.Type != catalog.TypeBase:
		return fmt.Errorf("%s is not a base image: its manifest says type %s", imagePath, manifest.Type)
	}
	return nil
}

// readRaw pulls the manifest bytes out of whichever image format this is.
// A writable .img is not handled: its environment comes from its .env sidecar.
func readRaw(imagePath string) ([]byte, error) {
	switch {
	case utils.IsSqf(imagePath):
		return squashfs.CatFile(imagePath, Path, 0)
	case utils.IsSif(imagePath):
		return sif.ReadFile(imagePath, Path)
	default:
		return nil, fmt.Errorf("%w: %s is not a CondaTainer image", tool.ErrCorrupt, imagePath)
	}
}
