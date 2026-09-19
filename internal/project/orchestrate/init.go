package orchestrate

import (
	"context"
	"fmt"

	"github.com/Justype/condatainer/internal/config"
	"github.com/Justype/condatainer/internal/project/lock"
	"github.com/Justype/condatainer/internal/utils"
)

// Init creates root's cnt-lock/ and pins the configured default distro's base
// as the project's root, without scanning for declarations. A base that cannot
// be pinned leaves the empty project in place and is returned as the error.
func Init(ctx context.Context, root string) (*lock.Pinned, error) {
	if err := utils.MkdirAllShared(lock.Dir(root)); err != nil {
		return nil, err
	}
	if err := lock.EnsureIgnore(root); err != nil {
		return nil, err
	}
	current, err := lock.Load(root)
	if err != nil {
		return nil, err
	}
	pinned, err := lock.DeriveBase(root, current, config.ResolvedDefaultDistro(), lock.PinOptions{})
	if err != nil {
		return nil, fmt.Errorf("project created without a root: %w", err)
	}
	if pinned == nil {
		return nil, nil
	}
	RecordUpstream(ctx, root, current, pinned)
	if err := lock.Publish(root, current); err != nil {
		return nil, err
	}
	return pinned, nil
}
