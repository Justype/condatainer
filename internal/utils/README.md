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
**Permission Fixers:** `FixPermissionsDefault(path)`

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

### Template Helpers

```go
// Extract {var} → value from a concrete name matched against a #TARGET: pattern.
// E.g. MatchTemplateTarget("salmon/{ver}", "salmon/1.10.2") → {"ver":"1.10.2"}, true
vars, ok := utils.MatchTemplateTarget(pattern, concrete)

// Interpolate {key} tokens in a string from a vars map.
result := utils.InterpolateVars(template, vars)

// Sort version strings descending (newest first). Returns a new slice.
sorted := utils.SortVersionsDescending(versions)

// Render a #TARGET: pattern with {placeholder} tokens in bold-yellow.
styled := utils.HighlightTemplatePlaceholders(pattern)
```

### Version Constraint Helpers

```go
// Split "samtools/1.22.1>=1.10" → ("samtools/1.22.1", ">=", "1.10")
nv, op, minVer := utils.SplitDepConstraint(raw)

// Compare partial version strings ("1.10" == "1.10.0"). Returns -1/0/1.
cmp := utils.CompareVersions(a, b)

// True if installedVersion satisfies op+minVersion and does not exceed preferredVersion.
// preferredVersion="" skips the upper bound check.
ok := utils.DepSatisfiedByVersion(installed, op, minVersion, preferredVersion)
```

`#DEP:name/version>=min` semantics: the preferred version is the implicit upper bound. An installed version is accepted if `min <= installed <= preferred`.
