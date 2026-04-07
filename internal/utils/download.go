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
// Tries bioconda first, then conda-forge. Returns "" if not found or on error.
// The summary is truncated to at most 50 words.
func FetchCondaSummary(packageName string) string {
	client := &http.Client{Timeout: 10 * time.Second}
	for _, channel := range []string{"bioconda", "conda-forge"} {
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
		return truncateWords(result.Summary, 50)
	}
	return ""
}

// truncateWords returns s truncated to at most n words, appending "..." if truncated.
func truncateWords(s string, n int) string {
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
