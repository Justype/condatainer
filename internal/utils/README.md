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

## Downloads

```go
utils.DownloadFile(url, destPath)       // with progress bar
utils.DownloadExecutable(url, destPath) // sets exec permissions
```

## Script Parsing

These read *user* scripts and helper scripts, not recipes — a recipe is parsed by
`catalog.ParseRecipe`.

```go
// parseModuleLoad: also extract "module load" / "ml" lines as deps
deps, err := utils.GetDependenciesFromScript(scriptPath, parseModuleLoad)

description := utils.GetDescriptionFromScript(scriptPath)

// #TYPE: from an external build script, defaulting to "app"
kind, err := utils.GetExternalBuildTypeFromScript(scriptPath)

// extract scheduler directives and apply defaults (scheduler package)
specs, err := scheduler.ReadScriptSpecsFromPath(scriptPath)
```

### Version Helpers

```go
// Sort version strings descending (newest first). Returns a new slice.
sorted := utils.SortVersionsDescending(versions)

// Render a #TARGET: pattern with {placeholder} tokens in bold-yellow.
styled := utils.HighlightTemplatePlaceholders(pattern)
```

Names, versions and dependency constraints live in `catalog`, so one string
resolves one way everywhere: `catalog.Normalize`, `catalog.CompareVersions`,
`catalog.ParseDep` and `Dep.Satisfies`. Template matching and `#PH:`/`#VALUE:`
value lists likewise: `catalog.NewTemplate` and `catalog.ParseValues`.
