# Go validation

Run all Go checks from the repository root.

- Use `go vet ./...` for the standard Go analyzers.
- Use `staticcheck ./...` for the broader correctness and style checks.
- Use `rg --files -g '*.go' -0 | xargs -0 gopls check -severity=warning`
  to reproduce the warning-level diagnostics shown by the VS Code Go extension,
  including the `golang.org/x/tools/go/analysis/passes/nilness` analyzer.
- Run `go test ./...` after changing Go code.
- Run `gofmt` on every modified Go file before committing.

If the default Go build cache is not writable, set `GOCACHE` to a writable
temporary directory when running these commands.
