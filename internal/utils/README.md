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

**Print:** `PrintMessage`, `PrintSuccess`, `PrintWarning`, `PrintError`, `PrintHint`, `PrintNote`, `PrintDebug`

**Style:** `StyleName` (yellow), `StylePath`, `StyleAction`, `StyleCommand`, `StyleHint` (cyan), `StyleError` (red), `StyleSuccess` (green), `StyleWarning` (yellow)

**Modes:** `DebugMode`, `QuietMode`, `YesMode`

```go
utils.PrintMessage("Installing %s", utils.StyleName("package"))
utils.PrintSuccess("Build completed")
utils.ShouldAnswerYes()                      // true when YesMode or non-interactive
utils.ReadLineContext(ctx context.Context)    // read line with cancellation
```

## File Operations

**Checks:** `FileExists`, `DirExists`, `IsImg`, `IsSqf`, `IsSif`, `IsOverlay`
**Operations:** `EnsureDir`
**Permissions:** `PermFile` (0664), `PermDir` (0775), `PermExec` (0775)
**Permission Fixers:** `FixPermissionsDefault(path)`, `ShareToUGORecursive(path)`

## Downloads

```go
utils.DownloadFile(url, destPath)       // with progress bar
utils.DownloadExecutable(url, destPath) // sets exec permissions
```

## Script Parsing

```go
// parseModuleLoad: also extract "module load" / "ml" lines as deps
deps, err := utils.GetDependenciesFromScript(scriptPath, parseModuleLoad)

prompts, err := utils.GetInteractivePromptsFromScript(scriptPath)

// extract scheduler directives and apply defaults (scheduler package)
specs, err := scheduler.ReadScriptSpecsFromPath(scriptPath)
```
