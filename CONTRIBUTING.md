# Contributing

Thank you for helping improve this project! This document explains the two main contribution types:

- Editing `condatainer` (code/tooling changes)
- Adding `build-scripts` (new app/data/genome build definitions)

## Quick start

1. Fork the repository from the `dev` branch and create a feature branch off `dev`.
2. Make changes, run tests, and verify on HPC and/or locally if applicable.
3. Push your branch to your fork and open a Pull Request. Set the PR base/target to `dev` (do NOT open PRs against `main`).
4. In your PR description include testing steps, relevant logs, and any backward-compatibility notes.

> Branching model — `dev` = unstable/testing (for active development); `main` = stable/release-only (maintainer-managed).

## 1) Editing `condatainer`

What this covers
- Bug fixes, feature additions, refactors, tests, or docs related to the `condatainer` tool.

Development
- Run `make build` to build the binary
- Run `make test` or `go test ./...` to run unit tests

Testing
- Run existing tests before submitting changes
- Add tests for new functionality when applicable
- Test your changes on both HPC and local environments if applicable

Commit and PR notes
- Use clear commit messages describing the change.
- In your PR description, explain the change, why it's needed, and how to test it.
- Set the PR base/target to `dev` (unstable/testing). Do NOT open PRs against `main` unless a maintainer explicitly requests it.

## 2) Adding `build-scripts`

What this covers
- New `build-scripts` for building genome reference images, environment definitions, or packaging steps used by the project.

Where to add
- For system/text editor apps, add to `build-scripts/apps/<app-name>.def` (an apptainer definition file)
- For applications, add to `build-scripts/<app-name>/<version>` (version is a script not a folder)
- For references, add to `build-scripts/<assembly>/<data-type>/<version>` (sub-version is allowed)
  - Add the tool version if the index/data is tool-version dependent. (If you are not sure, add the tool version)
  - E.g., `build-scripts/grch38/star/2.7.11b/gencode47-101`
  - E.g., `build-scripts/grch38/bowtie2/ucsc_no_alt`

Guidelines
- Please check the [Build Script Manual](./docs/manuals/build_script.md) for naming conventions and available variables.

Testing
- Test the build script locally and on HPC if applicable.
- Ensure the resulting environment/image works as expected. (paths and environment variables)

PR notes
- Please include the sbatch running log or terminal output showing successful build and test.
- Open Pull Requests to the `dev` branch (base = `dev`). Do NOT target `main` unless requested by a maintainer.

## PR checklist (apply to both types)

- Title and description explaining change and how to test
- Small, focused commits with clear messages
- Updated or added tests where applicable
- PR base/target: `dev` (unstable/testing). `main` is reserved for stable releases and maintainer merges.

## Contacts / Maintainers

- If unsure about scope or breaking changes, open an issue first describing the planned work and tag maintainers.

Thanks again — contributions make this project better!
