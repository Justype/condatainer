package utils

import (
	"encoding/json"
	"fmt"
	"io"
	"net/http"
	"os"
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
	Name     string   `json:"name"`
	Channel  string   `json:"owner"` // anaconda.org uses "owner" for the channel name
	Summary  string   `json:"summary"`
	Versions []string `json:"versions"`
}

// SearchCondaPackages looks up conda packages on anaconda.org.
//
// Exact match (default): queries the per-channel package API
// (/package/<channel>/<name>) for each configured channel in order.
// This is fast and reliable even for packages missing from the search index.
//
// Fuzzy match (fuzzy=true): queries the search API (/search?name=<query>)
// which does server-side substring matching and returns up to 100 results,
// filtered to the given channels. The returned bool is true when the API
// result was capped at 100 entries (results may be incomplete).
//
// Versions within each result are sorted descending.
func SearchCondaPackages(name string, channels []string, fuzzy ...bool) ([]CondaSearchResult, bool, error) {
	client := &http.Client{Timeout: 10 * time.Second}

	if len(fuzzy) > 0 && fuzzy[0] {
		results, capped, err := searchCondaFuzzy(client, name, channels)
		return results, capped, err
	}
	results, err := searchCondaExact(client, name, channels)
	return results, false, err
}

// searchCondaExact queries each channel's package API in order and returns the
// first channel that has the package — matching install behaviour.
func searchCondaExact(client *http.Client, name string, channels []string) ([]CondaSearchResult, error) {
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
			Name     string   `json:"name"`
			Summary  string   `json:"summary"`
			Versions []string `json:"versions"`
		}
		decodeErr := json.NewDecoder(resp.Body).Decode(&pkg)
		resp.Body.Close()
		if decodeErr != nil || pkg.Name != name {
			continue
		}
		return []CondaSearchResult{{
			Name:     pkg.Name,
			Channel:  ch,
			Summary:  TruncateWords(pkg.Summary, 100),
			Versions: SortVersionsDescending(pkg.Versions),
		}}, nil
	}
	return nil, nil
}

// searchCondaFuzzy queries the anaconda.org search API for a substring match.
// The API does server-side substring matching on package names and caps results at 100.
// Results are filtered to the given channels (all channels if empty).
// Returns capped=true when the raw API response hit the 100-result limit.
func searchCondaFuzzy(client *http.Client, name string, channels []string) ([]CondaSearchResult, bool, error) {
	apiURL := "https://api.anaconda.org/search?name=" + name
	resp, err := client.Get(apiURL)
	if err != nil {
		return nil, false, fmt.Errorf("search request failed: %w", err)
	}
	defer resp.Body.Close()
	if resp.StatusCode != http.StatusOK {
		return nil, false, fmt.Errorf("search API returned HTTP %d", resp.StatusCode)
	}

	var all []CondaSearchResult
	if err := json.NewDecoder(resp.Body).Decode(&all); err != nil {
		return nil, false, fmt.Errorf("failed to parse search response: %w", err)
	}

	capped := len(all) >= 100

	allowed := make(map[string]bool, len(channels))
	for _, ch := range channels {
		allowed[ch] = true
	}

	var filtered []CondaSearchResult
	for _, r := range all {
		if len(channels) > 0 && !allowed[r.Channel] {
			continue
		}
		r.Versions = SortVersionsDescending(r.Versions)
		filtered = append(filtered, r)
	}
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
