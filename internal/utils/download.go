package utils

import (
	"encoding/json"
	"fmt"
	"io"
	"net/http"
	"os"
	"runtime"
	"sort"
	"strings"
	"time"
)

// DownloadFile downloads a file from url to destPath.
func DownloadFile(url, destPath string) error {
	client := &http.Client{
		Timeout: 5 * time.Minute,
	}

	resp, err := client.Get(url)
	if err != nil {
		return fmt.Errorf("failed to fetch %s: %w", url, err)
	}
	defer resp.Body.Close()

	if resp.StatusCode != http.StatusOK {
		return fmt.Errorf("failed to download: HTTP %d", resp.StatusCode)
	}

	tmpPath := destPath + ".tmp"
	file, err := CreateFileWritable(tmpPath)
	if err != nil {
		return fmt.Errorf("failed to create file: %w", err)
	}

	_, err = io.Copy(file, resp.Body)
	file.Close()
	if err != nil {
		os.Remove(tmpPath)
		return fmt.Errorf("failed to write file: %w", err)
	}

	if err := os.Rename(tmpPath, destPath); err != nil {
		os.Remove(tmpPath)
		return fmt.Errorf("failed to rename file: %w", err)
	}

	// Set permissions for downloaded file (664 for group-writable)
	if err := os.Chmod(destPath, PermFile); err != nil {
		return fmt.Errorf("failed to set file permissions: %w", err)
	}

	return nil
}

// URLExists checks if a URL exists (returns HTTP 200) using a HEAD request.
func URLExists(url string) bool {
	client := &http.Client{
		Timeout: 30 * time.Second,
	}

	resp, err := client.Head(url)
	if err != nil {
		return false
	}
	defer resp.Body.Close()

	return resp.StatusCode == http.StatusOK
}

// FetchCondaSummary queries the anaconda.org API for a package summary.
// Channels are tried in the order provided; returns "" if not found or on error.
// The summary is truncated to at most 100 words.
func FetchCondaSummary(packageName string, channels []string) string {
	client := &http.Client{Timeout: 10 * time.Second}
	for _, channel := range channels {
		url := "https://api.anaconda.org/package/" + channel + "/" + packageName
		resp, err := client.Get(url)
		if err != nil || resp.StatusCode != http.StatusOK {
			if resp != nil {
				resp.Body.Close()
			}
			continue
		}
		var result struct {
			Summary string `json:"summary"`
		}
		err = json.NewDecoder(resp.Body).Decode(&result)
		resp.Body.Close()
		if err != nil || result.Summary == "" {
			continue
		}
		return TruncateWords(result.Summary, 100)
	}
	return ""
}

// CondaSearchResult holds a single result from the anaconda.org search API.
type CondaSearchResult struct {
	Name      string   `json:"name"`
	Channel   string   `json:"owner"` // anaconda.org uses "owner" for the channel name
	Summary   string   `json:"summary"`
	Versions  []string `json:"versions"`
	Platforms []string `json:"conda_platforms"` // supported conda platforms, e.g. ["linux-64", "osx-arm64"]
}

// SearchCondaPackages looks up conda packages on anaconda.org.
//
// Exact match (default): queries /package/<channel>/<name> for each channel in
// order and returns the first hit — matching install behaviour.
//
// Fuzzy match (fuzzy=true): issues two requests to /search?name=<query> —
// one for the current platform (e.g. linux-64) and one for noarch — then
// merges and deduplicates the results. Two requests are required because the
// API only accepts a single platform value and silently drops noarch-only
// packages when a platform filter is set. Returns capped=true if either
// request hit the limit (results may be incomplete).
//
// Versions within each result are sorted descending.
func SearchCondaPackages(name string, channels []string, fuzzy bool, limit int) ([]CondaSearchResult, bool, error) {
	client := &http.Client{Timeout: 10 * time.Second}

	if fuzzy {
		return searchCondaFuzzy(client, name, channels, limit)
	}
	results, err := searchCondaExact(client, name, channels)
	return results, false, err
}

// CurrentCondaPlatform returns the conda platform string for the current OS and
// architecture, e.g. "linux-64", "linux-aarch64", "osx-arm64".
// Note: Linux ARM64 is "linux-aarch64" while macOS ARM64 is "osx-arm64".
func CurrentCondaPlatform() string {
	goos := runtime.GOOS
	goarch := runtime.GOARCH
	condaOS := map[string]string{"linux": "linux", "darwin": "osx", "windows": "win"}[goos]
	if condaOS == "" {
		return ""
	}
	var condaArch string
	switch goarch {
	case "amd64":
		condaArch = "64"
	case "386":
		condaArch = "32"
	case "arm64":
		if goos == "linux" {
			condaArch = "aarch64"
		} else {
			condaArch = "arm64"
		}
	case "ppc64le":
		condaArch = "ppc64le"
	default:
		return ""
	}
	return condaOS + "-" + condaArch
}

// searchCondaExact queries each channel's package API in order and returns the
// first channel that has the package — matching install behaviour.
// Versions are filtered to those available for the current platform (or noarch).
// The result includes the list of supported conda platforms.
func searchCondaExact(client *http.Client, name string, channels []string) ([]CondaSearchResult, error) {
	currentPlatform := CurrentCondaPlatform()
	for _, ch := range channels {
		pkgURL := "https://api.anaconda.org/package/" + ch + "/" + name
		resp, err := client.Get(pkgURL)
		if err != nil || resp.StatusCode != http.StatusOK {
			if resp != nil {
				resp.Body.Close()
			}
			continue
		}
		var pkg struct {
			Name      string   `json:"name"`
			Summary   string   `json:"summary"`
			Versions  []string `json:"versions"`
			Platforms []string `json:"conda_platforms"`
			Files     []struct {
				Version string `json:"version"`
				Attrs   struct {
					Subdir string `json:"subdir"`
				} `json:"attrs"`
			} `json:"files"`
		}
		decodeErr := json.NewDecoder(resp.Body).Decode(&pkg)
		resp.Body.Close()
		if decodeErr != nil || pkg.Name != name {
			continue
		}
		// Filter versions to those available for the current platform or noarch.
		versions := pkg.Versions
		if currentPlatform != "" && len(pkg.Files) > 0 {
			seen := make(map[string]bool)
			for _, f := range pkg.Files {
				subdir := f.Attrs.Subdir
				if subdir == currentPlatform || subdir == "noarch" {
					seen[f.Version] = true
				}
			}
			versions = make([]string, 0, len(seen))
			for v := range seen {
				versions = append(versions, v)
			}
		}
		return []CondaSearchResult{{
			Name:      pkg.Name,
			Channel:   ch,
			Summary:   TruncateWords(pkg.Summary, 100),
			Versions:  SortVersionsDescending(versions),
			Platforms: pkg.Platforms,
		}}, nil
	}
	return nil, nil
}

// searchCondaFuzzy queries the anaconda.org /search API twice — once for the
// current platform and once for noarch — then merges and deduplicates results.
// Results are channel-filtered client-side. Returns capped=true if either
// request hit the limit.
func searchCondaFuzzy(client *http.Client, name string, channels []string, limit int) ([]CondaSearchResult, bool, error) {
	if limit <= 0 {
		limit = 100
	}
	// The API only accepts a single platform value and drops noarch-only packages
	// when platform is set. Query platform and noarch separately, then merge.
	platform := CurrentCondaPlatform()
	platforms := []string{"noarch"}
	if platform != "" {
		platforms = append([]string{platform}, platforms...)
	}

	// index by channel/name so versions from both requests can be merged
	type entry struct {
		result   CondaSearchResult
		versions map[string]bool
	}
	index := make(map[string]*entry)
	var order []string // preserve insertion order for stable sort base
	var capped bool
	for _, plat := range platforms {
		apiURL := fmt.Sprintf("https://api.anaconda.org/search?name=%s&limit=%d&platform=%s", name, limit, plat)
		resp, err := client.Get(apiURL)
		if err != nil || resp.StatusCode != http.StatusOK {
			if resp != nil {
				resp.Body.Close()
			}
			continue
		}
		var batch []CondaSearchResult
		json.NewDecoder(resp.Body).Decode(&batch) //nolint:errcheck
		resp.Body.Close()
		if len(batch) >= limit {
			capped = true
		}
		for _, r := range batch {
			key := r.Channel + "/" + r.Name
			if e, ok := index[key]; ok {
				// Merge versions from the second request
				for _, v := range r.Versions {
					e.versions[v] = true
				}
			} else {
				vset := make(map[string]bool, len(r.Versions))
				for _, v := range r.Versions {
					vset[v] = true
				}
				index[key] = &entry{result: r, versions: vset}
				order = append(order, key)
			}
		}
	}

	allowed := make(map[string]bool, len(channels))
	for _, ch := range channels {
		allowed[ch] = true
	}

	var filtered []CondaSearchResult
	for _, key := range order {
		e := index[key]
		if len(channels) > 0 && !allowed[e.result.Channel] {
			continue
		}
		vs := make([]string, 0, len(e.versions))
		for v := range e.versions {
			vs = append(vs, v)
		}
		e.result.Versions = SortVersionsDescending(vs)
		filtered = append(filtered, e.result)
	}

	// Sort results alphabetically by name
	sort.Slice(filtered, func(i, j int) bool {
		return filtered[i].Name < filtered[j].Name
	})

	return filtered, capped, nil
}

// TruncateWords returns s truncated to at most n words, appending "..." if truncated.
func TruncateWords(s string, n int) string {
	words := strings.Fields(s)
	if len(words) <= n {
		return s
	}
	return strings.Join(words[:n], " ") + "..."
}

// DownloadExecutable downloads a file and sets it as executable (PermExec).
func DownloadExecutable(url, destPath string) error {
	if err := DownloadFile(url, destPath); err != nil {
		return err
	}

	if err := os.Chmod(destPath, PermExec); err != nil {
		return fmt.Errorf("failed to set executable permissions: %w", err)
	}

	return nil
}
