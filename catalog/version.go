package catalog

import "strings"

// CompareVersions orders two versions naturally: digit runs compare as numbers,
// the runs between them as text. Returns -1, 0 or 1.
//
// A prefix sorts below what extends it (1.16 < 1.16.0), except a pre-release,
// which sorts below its own release:
//
//	1.1.0-alpha < 1.1.0-beta < 1.1.0-rc1 < 1.1.0-rc2 < 1.1.0
//
// Labels order alphabetically, so -dev sorts by spelling rather than by meaning.
func CompareVersions(a, b string) int {
	i, j := 0, 0
	for i < len(a) || j < len(b) {
		// Whichever ran out is the shorter version, and wins only if what
		// remains of the other is a pre-release.
		if i == len(a) {
			return endsBefore(b[j:])
		}
		if j == len(b) {
			return -endsBefore(a[i:])
		}

		si, sj := i, j
		for i < len(a) && !isDigit(a[i]) {
			i++
		}
		for j < len(b) && !isDigit(b[j]) {
			j++
		}
		if c := strings.Compare(a[si:i], b[sj:j]); c != 0 {
			return c
		}

		si, sj = i, j
		for i < len(a) && isDigit(a[i]) {
			i++
		}
		for j < len(b) && isDigit(b[j]) {
			j++
		}
		if c := compareDigits(a[si:i], b[sj:j]); c != 0 {
			return c
		}
	}
	return 0
}

// endsBefore reports how a version that has ended compares to one continuing
// with rest: greater when rest is a pre-release suffix, lesser otherwise.
//
// A pre-release is a dash followed by a letter. The letter separates -rc1 from
// an ordinary trailing component like RStudio's 2026.06.0-242.
func endsBefore(rest string) int {
	if len(rest) > 1 && rest[0] == '-' && isLetter(rest[1]) {
		return 1
	}
	return -1
}

// compareDigits compares two digit runs as numbers of any length: once leading
// zeros are gone the longer run is the larger number.
func compareDigits(a, b string) int {
	a, b = strings.TrimLeft(a, "0"), strings.TrimLeft(b, "0")
	if len(a) != len(b) {
		if len(a) < len(b) {
			return -1
		}
		return 1
	}
	return strings.Compare(a, b)
}

func isDigit(c byte) bool  { return c >= '0' && c <= '9' }
func isLetter(c byte) bool { return c >= 'a' && c <= 'z' || c >= 'A' && c <= 'Z' }
