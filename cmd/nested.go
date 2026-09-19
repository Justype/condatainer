package cmd

import (
	"context"
	"fmt"
	"slices"

	"github.com/Justype/condatainer/internal/build"
	"github.com/Justype/condatainer/internal/config"
	"github.com/Justype/condatainer/internal/libexec"
	"github.com/Justype/condatainer/internal/runtime/container"
	"github.com/Justype/condatainer/internal/utils"
)

// nestedPlan is how one container gets apptainer for nested running.
type nestedPlan struct {
	BindLibexec bool   // bind libexec/, which has an apptainer
	Overlay     string // installed apptainer version to mount, "" for none
	Build       bool   // no provider exists and the overlay must be built first
}

// planNested applies nested_run: libexec's apptainer wins and is bound; else the
// newest installed apptainer overlay is mounted; else "true" builds one and
// "auto" does nothing.
func planNested(mode string, insideContainer, libexecHasApptainer bool, installed []string) nestedPlan {
	if mode == config.NestedRunFalse || insideContainer {
		return nestedPlan{}
	}
	if libexecHasApptainer {
		return nestedPlan{BindLibexec: true}
	}
	if len(installed) > 0 {
		return nestedPlan{Overlay: utils.SortVersionsDescending(installed)[0]}
	}
	return nestedPlan{Build: mode == config.NestedRunTrue}
}

// insideContainer is config.IsInsideContainer, replaceable by tests.
var insideContainer = config.IsInsideContainer

// currentNestedPlan gathers the plan's inputs from the running system.
func currentNestedPlan() nestedPlan {
	return planNested(config.Global.NestedRun, insideContainer(),
		libexec.Installed("apptainer"), build.InstalledVersions(nil)("apptainer"))
}

// nestedRun adds what nested running needs to one e/exec/run launch: the
// overlays to mount, and whether to bind libexec. When nested_run is "true"
// and nothing provides apptainer, it builds the overlay first, like the base
// image, and fails if it cannot.
func nestedRun(ctx context.Context, overlays []string) ([]string, bool, error) {
	plan := currentNestedPlan()
	if plan.Build {
		utils.PrintNote("Building the apptainer overlay for nested running (nested_run: true)")
		err := buildNestedApptainer(ctx)
		build.InvalidateInstalledOverlays()
		if err != nil {
			return nil, false, fmt.Errorf("nested_run is true but apptainer cannot be provided for nested running: %w; install it with `condatainer install apptainer`, or set nested_run to auto", err)
		}
		if plan = currentNestedPlan(); plan.Overlay == "" && !plan.BindLibexec {
			return nil, false, fmt.Errorf("nested_run is true but the apptainer overlay build installed nothing")
		}
	}
	if plan.Overlay != "" {
		found, err := container.ResolveOverlayPaths([]string{"apptainer/" + plan.Overlay})
		if err != nil {
			if config.Global.NestedRun == config.NestedRunTrue {
				return nil, false, err
			}
			return overlays, false, nil
		}
		path := found[0]
		// First in the list, so its bin/ comes after every other overlay's on PATH
		// (BuildPathEnv prepends in list order) and shadows none of their tools.
		if !slices.Contains(overlays, path) {
			overlays = append([]string{path}, overlays...)
		}
	}
	return overlays, plan.BindLibexec, nil
}

// buildNestedApptainer builds the newest apptainer from the configured channels,
// locally, the way the default base image is built.
func buildNestedApptainer(ctx context.Context) error {
	name, _, err := solveCreateName(ctx, "apptainer")
	if err != nil {
		return err
	}
	imagesDir, err := config.GetWritableImagesDir()
	if err != nil {
		return err
	}
	obj, err := build.NewBuildObject(ctx, name, false, imagesDir, false)
	if err != nil {
		return err
	}
	graph, err := build.NewBuildGraph(ctx, []*build.BuildObject{obj}, imagesDir, false, false)
	if err != nil {
		return err
	}
	return graph.Run(ctx)
}
