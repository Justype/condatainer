package web

import "embed"

// Files contains the dashboard web assets embedded into the binary.
//
//go:embed *
var Files embed.FS
