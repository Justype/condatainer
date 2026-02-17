# utils

Shared utilities for console output, file operations, downloads, and script parsing.

## Architecture

```
console.go   Styled console output, logging, user prompts
download.go  HTTP downloads with progress bars
files.go     File and directory utilities
parser.go    Build script metadata parsing
```

## Console Output

**Print Functions:** `PrintMessage`, `PrintSuccess`, `PrintWarning`, `PrintError`, `PrintNote`, `PrintDebug`

**Style Functions:** `StyleError` (red), `StyleSuccess` (green), `StyleWarning` (yellow), `StyleHint` (cyan), `StyleInfo` (magenta), `StylePath`, `StyleName`, `StyleNumber`, `StyleAction`, `StyleCommand`

**Modes:** `DebugMode`, `QuietMode`, `YesMode`

```go
utils.PrintMessage("Installing %s", utils.StyleName("package"))
utils.PrintSuccess("Build completed")
utils.PrintDebug("Path: %s", utils.StylePath(path))
yes, err := utils.AskYesNo("Continue?", defaultYes=true)
```

## File Operations

**Checks:** `FileExists`, `DirExists`, `IsImg`, `IsSqf`, `IsOverlay`  
**Operations:** `EnsureDir`, `CopyFile`, `MoveFile`, `DeleteFile`  
**Permissions:** `PermFile` (0664), `PermDir` (0775), `PermExec` (0775)

```go
utils.FileExists("/path/file")
utils.IsImg("env.img")
utils.EnsureDir("/path/dir")
```

## Downloads

```go
utils.DownloadFile(url, destPath)           // With progress bar
utils.DownloadExecutable(url, destPath)     // Sets exec permissions
```

## Script Parsing

**Dependencies:**
```go
deps, err := utils.GetDependenciesFromScript(scriptPath)
// Extracts #DEP:name/version and module load commands
```

**Interactive Prompts:**
```go
prompts, err := utils.GetInteractivePromptsFromScript(scriptPath)
// Extracts #INTERACTIVE:prompt directives
```

**Scheduler Specs** (from `internal/scheduler`):
```go
specs, err := scheduler.ReadScriptSpecsFromPath(scriptPath)
// Extracts #SBATCH/#PBS/#BSUB directives and resource requirements (HTCondor uses native .sub files)
```
