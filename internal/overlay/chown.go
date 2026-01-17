package overlay

import (
	"bufio"
	"bytes"
	"fmt"
	"os/exec"
	"path/filepath"
	"regexp"
	"strconv"
	"strings"

	"github.com/Justype/condatainer/internal/utils"
)

// ChownRecursively changes the UID/GID of files inside an unmounted overlay image.
// It maps the provided internalPath (e.g., "/") to the overlay's "/upper" structure.
func ChownRecursively(imagePath string, uid, gid int, internalPath string) error {
	// 1. Dependency Check
	if err := checkDependencies([]string{"debugfs"}); err != nil {
		return err
	}

	absPath, err := filepath.Abs(imagePath)
	if err != nil {
		return fmt.Errorf("failed to resolve path %s: %w", imagePath, err)
	}

	// 2. Map logical path to internal /upper structure
	// e.g. "/" -> "/upper", "/etc" -> "/upper/etc"
	targetPath := filepath.Join("/upper", strings.TrimPrefix(internalPath, "/"))

	utils.PrintMessage("Scanning %s inside %s...",
		utils.StylePath(targetPath),
		utils.StylePath(filepath.Base(absPath)),
	)

	// 3. Recursive Discovery
	inodes, err := scanInodes(absPath, targetPath)
	if err != nil {
		return &Error{
			Op:      "scan inodes",
			Path:    absPath,
			Tool:    "debugfs",
			Output:  err.Error(),
			BaseErr: err,
		}
	}

	if len(inodes) == 0 {
		utils.PrintWarning("No inodes found to modify at %s.", utils.StylePath(targetPath))
		return nil
	}

	// 4. Batch Update
	utils.PrintMessage("Modifying %s inodes: UID=%s GID=%s...",
		utils.StyleNumber(len(inodes)),
		utils.StyleNumber(uid), // Using Number style for consistency
		utils.StyleNumber(gid),
	)

	var cmds []string

	// A. Ensure skeleton paths are accessible (if we are chowning root)
	if internalPath == "/" || internalPath == "" {
		skeletonPaths := []string{"/upper", "/work", "/work/work"}
		for _, p := range skeletonPaths {
			cmds = append(cmds, fmt.Sprintf("sif %s uid %d", p, uid))
			cmds = append(cmds, fmt.Sprintf("sif %s gid %d", p, gid))
		}
	}

	// B. Fix discovered inodes
	for _, inode := range inodes {
		cmds = append(cmds, fmt.Sprintf("sif <%s> uid %d", inode, uid))
		cmds = append(cmds, fmt.Sprintf("sif <%s> gid %d", inode, gid))
	}
	cmds = append(cmds, "quit")

	// C. Execute Batch
	script := strings.Join(cmds, "\n")
	cmd := exec.Command("debugfs", "-w", absPath)
	cmd.Stdin = strings.NewReader(script)

	if out, err := cmd.CombinedOutput(); err != nil {
		return &Error{
			Op:      "chown inodes",
			Path:    absPath,
			Tool:    "debugfs",
			Output:  string(out),
			BaseErr: err,
		}
	}

	utils.PrintSuccess("Permissions updated successfully.")
	return nil
}

// scanInodes implements a BFS walker to find all inodes under startPath.
func scanInodes(imgPath, startPath string) ([]string, error) {
	// A. Get Start Inode via 'stat'
	cmd := exec.Command("debugfs", "-R", fmt.Sprintf("stat %s", startPath), imgPath)
	out, err := cmd.Output()
	if err != nil {
		// This usually means the path doesn't exist in the overlay layer yet
		// (e.g. user hasn't written to it), so we return empty without error.
		return nil, fmt.Errorf("path %s not found in overlay layer (it may only exist in base image)", startPath)
	}

	reInode := regexp.MustCompile(`Inode:\s+(\d+)`)
	match := reInode.FindStringSubmatch(string(out))
	if match == nil {
		return nil, fmt.Errorf("could not parse inode for %s", startPath)
	}
	startInode := match[1]

	// B. Walk the Tree
	visited := make(map[string]bool)
	toVisit := []string{startInode}
	results := []string{startInode}

	count := 0

	for len(toVisit) > 0 {
		current := toVisit[0]
		toVisit = toVisit[1:]

		if visited[current] {
			continue
		}
		visited[current] = true

		// List contents: ls -l <INODE>
		cmdLs := exec.Command("debugfs", "-R", fmt.Sprintf("ls -l <%s>", current), imgPath)
		outLs, _ := cmdLs.Output()

		scanner := bufio.NewScanner(bytes.NewReader(outLs))
		for scanner.Scan() {
			line := strings.TrimSpace(scanner.Text())
			parts := strings.Fields(line)

			// Expected debugfs output format: "  2   40755 (2)      0      0    ..."
			// We need at least inode, mode, and name
			if len(parts) < 6 {
				continue
			}

			inode := parts[0]
			mode := parts[1] // e.g. "40755" (dir) or "100644" (file)
			name := parts[len(parts)-1]

			if name == "." || name == ".." {
				continue
			}

			// Validate inode is a number
			if _, err := strconv.Atoi(inode); err != nil {
				continue
			}

			results = append(results, inode)

			// Check if directory (Mode starts with 4 or 04)
			// 40755 -> Directory
			// 100644 -> File
			if strings.HasPrefix(mode, "4") || strings.HasPrefix(mode, "04") {
				toVisit = append(toVisit, inode)
			}

			count++
			if count%1000 == 0 {
				fmt.Printf("\r    Scanning inodes... %s found", utils.StyleNumber(count))
			}
		}
	}
	fmt.Printf("\r%s\r", strings.Repeat(" ", 50)) // Clear line

	return results, nil
}
