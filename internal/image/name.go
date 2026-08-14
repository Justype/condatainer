package image

import "strings"

// EncodeArtifactName maps a normalized artifact name to the single-component
// spelling used by image files, provenance capsules, and store entries.
func EncodeArtifactName(name string) string {
	return strings.ReplaceAll(strings.Trim(name, "/"), "/", "--")
}

// DecodeArtifactName reverses EncodeArtifactName for artifact names.
func DecodeArtifactName(name string) string {
	return strings.ReplaceAll(name, "--", "/")
}
