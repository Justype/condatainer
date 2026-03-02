---
name: "\U0001F4A1 Feature Request"
about: Suggest a new feature or improvement for CondaTainer
title: "[Feature] <summary>"
labels: enhancement
assignees: ''

---

## Problem / Motivation

<!-- What limitation or need is driving this request? -->

## Proposed Solution

<!-- Describe your proposed solution or API change. Be as specific as possible. -->

## Affected Component

<!-- Check all that apply -->

**Module Management**
- [ ] `build` — build system, dependency graph (`#DEP` resolution)
- [ ] `avail` — remote script fetching
- [ ] `list` — listing installed modules / overlays
- [ ] `remove` — removing installed modules / overlays

**Overlay**
- [ ] `overlay info` — overlay image inspection
- [ ] `overlay create` / `o` — ext3 overlay creation (sparse, fakeroot, conda init)
- [ ] `overlay resize` — resizing existing overlay images
- [ ] `overlay check` — filesystem integrity verification

**Execution**
- [ ] `exec` — explicit overlay execution (read-only by default)
- [ ] `e` — shortcut exec with writable overlay and auto-load

**Script Runner**
- [ ] `run` (local) — local script execution with `#DEP` auto-resolution
- [ ] `run` (scheduler submit) — job submission to SLURM / PBS / LSF / HTCondor
- [ ] `check` — pre-flight dependency check for scripts

**Scheduler**
- [ ] Script parsing — reading `#SBATCH` / `#PBS` / `#BSUB` directives or `.sub` script
- [ ] Script generation — creating scheduler job scripts
- [ ] Info gathering — partition / queue / resource discovery

**Instance**
- [ ] `instance start` — starting named persistent Apptainer instances
- [ ] `instance stop` — stopping instances
- [ ] `instance exec` — executing commands inside running instances
- [ ] `instance list` — listing active instances

**Configuration**
- [ ] `config get` / `config set` — configuration key management
- [ ] `config init` — config file initialization
- [ ] `config env` — `CNT_*` environment variable handling

**Helper / Scripts**
- [ ] `helper` — fetching / updating helper scripts from `cnt-scripts`
- [ ] Apptainer wrapper — low-level Apptainer binary abstraction

**Other**
- [ ] CLI (new commands or parameters)
- [ ] Update (binary or `--base` base image)
- [ ] Shell completion (bash, zsh, fish)
- [ ] Other: <!-- describe -->

## Alternatives Considered

<!-- Any alternative approaches you've thought about? -->

## Additional Notes

<!-- Links to related issues, upstream tools, or examples from other projects -->
