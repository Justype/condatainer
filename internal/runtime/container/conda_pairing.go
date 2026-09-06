package container

import (
	"github.com/Justype/condatainer/internal/conda"
	"github.com/Justype/condatainer/internal/utils"
)

// PairedPackages returns every package installed at path — a .img or a bare
// .sqf — merged with an autoloaded snapshot's own packages (LookupSnapshot)
// when path is a .img that pairs with one, snapshot first so path's own
// entries win on conflict. This is the one entry point #IMG_PACKAGES: checks
// and similar callers need; they should not call
// conda.ListCondaPackages/ListCondaPackagesSqf directly and reimplement the
// pairing themselves. Returns (nil, nil) only when neither side has a conda
// environment at all.
func PairedPackages(path string) (map[string]string, error) {
	if utils.IsSqf(path) {
		return conda.ListCondaPackagesSqf(path)
	}
	imgPkgs, err := conda.ListCondaPackages(path)
	if err != nil {
		return nil, err
	}
	var snapshotPkgs map[string]string
	if utils.IsImg(path) {
		if lookup := LookupSnapshot(path); lookup.Path != "" {
			if snapshotPkgs, err = conda.ListCondaPackagesSqf(lookup.Path); err != nil {
				return nil, err
			}
		}
	}
	if snapshotPkgs == nil && imgPkgs == nil {
		return nil, nil
	}
	merged := make(map[string]string, len(snapshotPkgs)+len(imgPkgs))
	for name, ver := range snapshotPkgs {
		merged[name] = ver
	}
	for name, ver := range imgPkgs {
		merged[name] = ver
	}
	return merged, nil
}

// PairedInfo reads channels and explicitly-installed specs for path, merging
// a .img's own conda-meta/history with an autoloaded snapshot's (as one
// continuous log) when path is a .img that pairs with one (LookupSnapshot). A
// bare .sqf, or a .img with no pair, is read alone. This is the one entry
// point display code (`info`, the dashboard) needs — they should not call
// conda.ReadCondaInfo/ReadCondaInfoMerged directly and reimplement the
// pairing.
func PairedInfo(path, envPrefix string) *conda.CondaInfo {
	if utils.IsImg(path) {
		if lookup := LookupSnapshot(path); lookup.Path != "" {
			return conda.ReadCondaInfoMerged(lookup.Path, path, envPrefix)
		}
	}
	return conda.ReadCondaInfo(path, envPrefix)
}
