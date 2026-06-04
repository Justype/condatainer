package web

import "embed"

// Files contains the dashboard web assets embedded into the binary.
//
//go:embed *.html *.js *.css *.svg
var Files embed.FS
