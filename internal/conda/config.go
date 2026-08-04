package conda

import (
	"bytes"
	"fmt"
	"os"
	"path/filepath"
	"slices"

	"go.yaml.in/yaml/v3"
)

// ParseChannels returns the top-level channels list in priority order.
func ParseChannels(data []byte) ([]string, error) {
	if len(bytes.TrimSpace(data)) == 0 {
		return nil, nil
	}
	doc, err := parseYAMLDocument(data)
	if err != nil {
		return nil, err
	}
	value := mappingValue(doc.Content[0], "channels")
	if value == nil {
		return nil, nil
	}
	if value.Kind != yaml.SequenceNode {
		return nil, fmt.Errorf(".condarc channels must be a YAML list")
	}
	channels := make([]string, 0, len(value.Content))
	for _, item := range value.Content {
		if item.Kind != yaml.ScalarNode || item.Value == "" {
			return nil, fmt.Errorf(".condarc contains an invalid channel entry")
		}
		channels = append(channels, item.Value)
	}
	return channels, nil
}

// SetChannels replaces only the top-level channels value, retaining all other
// configuration nodes and comments.
func SetChannels(data []byte, channels []string) ([]byte, error) {
	doc, err := parseYAMLDocument(data)
	if err != nil {
		return nil, err
	}
	root := doc.Content[0]
	value := &yaml.Node{Kind: yaml.SequenceNode, Tag: "!!seq"}
	for _, channel := range channels {
		if channel == "" {
			return nil, fmt.Errorf("channel cannot be empty")
		}
		value.Content = append(value.Content, &yaml.Node{Kind: yaml.ScalarNode, Tag: "!!str", Value: channel})
	}

	for i := 0; i+1 < len(root.Content); i += 2 {
		if root.Content[i].Value == "channels" {
			root.Content[i+1] = value
			return encodeYAML(doc)
		}
	}
	root.Content = append(root.Content,
		&yaml.Node{Kind: yaml.ScalarNode, Tag: "!!str", Value: "channels"}, value)
	return encodeYAML(doc)
}

func ReadChannels(path string) ([]string, error) {
	data, err := os.ReadFile(path)
	if os.IsNotExist(err) {
		return nil, nil
	}
	if err != nil {
		return nil, err
	}
	return ParseChannels(data)
}

func WriteChannels(path string, channels []string) error {
	data, err := os.ReadFile(path)
	if os.IsNotExist(err) {
		data = nil
	} else if err != nil {
		return err
	}
	updated, err := SetChannels(data, channels)
	if err != nil {
		return err
	}
	return writeFileAtomic(path, updated, 0o644)
}

type ChannelChange int

const (
	ChannelAppend ChannelChange = iota
	ChannelPrepend
	ChannelRemove
)

// ChangeChannel updates the saved channel order and reports whether the
// channel was already present.
func (e *Environment) ChangeChannel(channel string, change ChannelChange) (bool, error) {
	if channel == "" {
		return false, fmt.Errorf("channel cannot be empty")
	}
	channels, err := ReadChannels(e.CondarcPath)
	if err != nil {
		return false, err
	}
	found := slices.Contains(channels, channel)
	updated := make([]string, 0, len(channels)+1)
	for _, existing := range channels {
		if existing != channel {
			updated = append(updated, existing)
		}
	}
	switch change {
	case ChannelRemove:
		if !found {
			return false, fmt.Errorf("channel %q is not configured", channel)
		}
	case ChannelPrepend:
		updated = append([]string{channel}, updated...)
	case ChannelAppend:
		updated = append(updated, channel)
	default:
		return false, fmt.Errorf("unknown channel change")
	}
	return found, WriteChannels(e.CondarcPath, updated)
}

func parseYAMLDocument(data []byte) (*yaml.Node, error) {
	if len(bytes.TrimSpace(data)) == 0 {
		return &yaml.Node{Kind: yaml.DocumentNode, Content: []*yaml.Node{{Kind: yaml.MappingNode, Tag: "!!map"}}}, nil
	}
	var doc yaml.Node
	if err := yaml.Unmarshal(data, &doc); err != nil {
		return nil, fmt.Errorf("invalid YAML: %w", err)
	}
	if len(doc.Content) != 1 || doc.Content[0].Kind != yaml.MappingNode {
		return nil, fmt.Errorf("configuration must be a YAML mapping")
	}
	return &doc, nil
}

func mappingValue(mapping *yaml.Node, key string) *yaml.Node {
	for i := 0; i+1 < len(mapping.Content); i += 2 {
		if mapping.Content[i].Value == key {
			return mapping.Content[i+1]
		}
	}
	return nil
}

func encodeYAML(doc *yaml.Node) ([]byte, error) {
	var out bytes.Buffer
	enc := yaml.NewEncoder(&out)
	enc.SetIndent(2)
	if err := enc.Encode(doc); err != nil {
		return nil, err
	}
	if err := enc.Close(); err != nil {
		return nil, err
	}
	return out.Bytes(), nil
}

func writeFileAtomic(path string, data []byte, defaultMode os.FileMode) error {
	dir := filepath.Dir(path)
	if err := os.MkdirAll(dir, 0o755); err != nil {
		return err
	}
	mode := defaultMode
	if info, err := os.Stat(path); err == nil {
		mode = info.Mode().Perm()
	}
	tmp, err := os.CreateTemp(dir, ".condatainer-*")
	if err != nil {
		return err
	}
	tmpPath := tmp.Name()
	defer os.Remove(tmpPath) //nolint:errcheck
	if err := tmp.Chmod(mode); err != nil {
		tmp.Close()
		return err
	}
	if _, err := tmp.Write(data); err != nil {
		tmp.Close()
		return err
	}
	if err := tmp.Sync(); err != nil {
		tmp.Close()
		return err
	}
	if err := tmp.Close(); err != nil {
		return err
	}
	return os.Rename(tmpPath, path)
}
