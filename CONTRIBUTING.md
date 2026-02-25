# Contributing

Thank you for helping improve this project! This document covers contributions to `condatainer` (CLI code and docs).

> For build scripts and helper scripts, contribute to the [`cnt-scripts`](https://github.com/Justype/cnt-scripts) repo instead.

## Repositories

| Repo | Contents |
|------|----------|
| [`condatainer`](https://github.com/Justype/condatainer) | CLI binary (Go code, docs) |
| [`cnt-scripts`](https://github.com/Justype/cnt-scripts) | Build scripts and helper scripts |

## Quick start

1. Fork this repository from the `dev` branch and create a feature branch off `dev`.
2. Make changes, run tests, and verify on HPC and/or locally if applicable.
3. Push your branch to your fork and open a Pull Request. Set the PR base/target to `dev` (do NOT open PRs against `main`).
4. In your PR description include testing steps, relevant logs, and any backward-compatibility notes.

> Branching model — `dev` = unstable/testing (for active development); `main` = stable/release-only (maintainer-managed).

## Editing `condatainer`

What this covers
- Bug fixes, feature additions, refactors, tests, or docs related to the `condatainer` tool.

Development
- Run `make build` to build the binary
- Run `make test` or `go test ./...` to run unit tests

Testing
- Run existing tests before submitting changes
- Add tests for new functionality when applicable
- Test your changes on both HPC and local environments if applicable

PR notes
- Use clear commit messages describing the change.
- In your PR description, explain the change, why it's needed, and how to test it.
- Set the PR base/target to `dev` (unstable/testing). Do NOT open PRs against `main` unless a maintainer explicitly requests it.

## PR checklist

- Title and description explaining change and how to test
- Small, focused commits with clear messages
- Updated or added tests where applicable
- PR base/target: `dev` (unstable/testing). `main` is reserved for stable releases and maintainer merges.

## Contacts / Maintainers

- If unsure about scope or breaking changes, open an issue first describing the planned work and tag maintainers.

Thanks again — contributions make this project better!
