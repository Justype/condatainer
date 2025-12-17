# Contributing

Thank you for helping improve this project! This document explains the two main contribution types:

- Editing `condatainer/modgen` (code/tooling changes)
- Adding `build-scripts` (new data/genome build definitions)

## Quick start

1. Fork the repository and create a feature branch.
2. Make changes, run a small test on both HPC and locally if applicable.
3. Push your branch and open a Pull Request with a clear description and testing notes

## 1) Editing `condatainer/modgen`

What this covers
- Bug fixes, feature additions, refactors, tests, or docs related to the `condatainer/modgen` tool.

Guidelines
- Do not use external dependencies. Make the script self-contained.
- Make sure condatainer and modgen have the same parameters and behavior where applicable.
  - Don't forget to update both scripts if you change shared logic.

Testing
- There is no formal test suite yet, but please test your changes on both HPC and local environments if applicable.

Commit and PR notes
- Use clear commit messages describing the change.
- In your PR description, explain the change, why it's needed, and how to test it.

## 2) Adding `build-scripts`

What this covers
- New `build-scripts` for building genome reference images, environment definitions, or packaging steps used by the project.

Where to add
- For applications, add to `build-scripts/<app-name>/<version>` (version is a script not a folder)
- For references, add to `build-scripts/<assembly>/<data-type>/<version>` (sub-version is allowed)
  - Add the tool version if the index/data is tool-version dependent. (If you are not sure, add the tool version)
  - E.g., `build-scripts/grch38/star/2.7.11b/gencode47-101`
  - E.g., `build-scripts/grch38/bowtie2/ucsc_no_alt`

Guidelines
- Please check the [Build Script Manual](assets/MANUAL_BUILD_SCRIPT.md) for naming conventions and available variables.

Testing
- Test the build script locally and on HPC if applicable.
- Ensure the resulting environment/image works as expected. (paths and environment variables)

PR notes
- Please include the sbatch running log or terminal output showing successful build and test.

## PR checklist (apply to both types)

- Title and description explaining change and how to test
- Small, focused commits with clear messages
- Updated or added tests where applicable

## Contacts / Maintainers

- If unsure about scope or breaking changes, open an issue first describing the planned work and tag maintainers.

Thanks again — contributions make this project better!
