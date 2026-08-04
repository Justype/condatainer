package conda

import (
	"bufio"
	"bytes"
	"context"
	"fmt"
	"os"
	"regexp"
	"strings"
)

var matchSpecName = regexp.MustCompile(`^(?:[^:[:space:]]+::)?([A-Za-z0-9_.-]+)`)

func matchSpecPackage(spec string) (string, error) {
	m := matchSpecName.FindStringSubmatch(strings.TrimSpace(spec))
	if len(m) != 2 {
		return "", fmt.Errorf("invalid package MatchSpec %q", spec)
	}
	return m[1], nil
}

func ReadPins(path string) ([]string, error) {
	data, err := os.ReadFile(path)
	if os.IsNotExist(err) {
		return nil, nil
	}
	if err != nil {
		return nil, err
	}
	var pins []string
	scanner := bufio.NewScanner(bytes.NewReader(data))
	for scanner.Scan() {
		line := strings.TrimSpace(scanner.Text())
		if line != "" && !strings.HasPrefix(line, "#") {
			pins = append(pins, line)
		}
	}
	return pins, scanner.Err()
}

func AddPin(path, spec string) error {
	name, err := matchSpecPackage(spec)
	if err != nil {
		return err
	}
	pins, err := ReadPins(path)
	if err != nil {
		return err
	}
	updated := make([]string, 0, len(pins)+1)
	for _, pin := range pins {
		pinName, parseErr := matchSpecPackage(pin)
		if parseErr != nil || pinName != name {
			updated = append(updated, pin)
		}
	}
	updated = append(updated, strings.TrimSpace(spec))
	return writePins(path, updated)
}

// ResolvePin converts a bare package name to an exact installed-version pin.
// An explicit MatchSpec is returned unchanged.
func (e *Environment) ResolvePin(ctx context.Context, spec string) (string, error) {
	name, err := matchSpecPackage(spec)
	if err != nil {
		return "", err
	}
	if strings.TrimSpace(spec) != name {
		return strings.TrimSpace(spec), nil
	}
	version, err := e.InstalledVersion(ctx, name)
	if err != nil {
		return "", err
	}
	return fmt.Sprintf("%s ==%s", name, version), nil
}

func (e *Environment) Pin(ctx context.Context, spec string) (string, error) {
	resolved, err := e.ResolvePin(ctx, spec)
	if err != nil {
		return "", err
	}
	return resolved, AddPin(e.PinnedPath, resolved)
}

func (e *Environment) Unpin(packageName string) error {
	return RemovePin(e.PinnedPath, packageName)
}

func RemovePin(path, packageName string) error {
	name, err := matchSpecPackage(packageName)
	if err != nil {
		return err
	}
	pins, err := ReadPins(path)
	if err != nil {
		return err
	}
	found := false
	updated := make([]string, 0, len(pins))
	for _, pin := range pins {
		pinName, parseErr := matchSpecPackage(pin)
		if parseErr == nil && pinName == name {
			found = true
			continue
		}
		updated = append(updated, pin)
	}
	if !found {
		return fmt.Errorf("package %q is not pinned", name)
	}
	return writePins(path, updated)
}

func writePins(path string, pins []string) error {
	var out strings.Builder
	for _, pin := range pins {
		out.WriteString(pin)
		out.WriteByte('\n')
	}
	return writeFileAtomic(path, []byte(out.String()), 0o644)
}
