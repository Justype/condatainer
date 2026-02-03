package overlay

import (
	"bufio"
	"bytes"
	"context"
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
func ChownRecursively(ctx context.Context, imagePath string, uid, gid int, internalPath string) error {
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

	utils.PrintMessage("Scanning %s inside %s (uid=%d gid=%d)",
		utils.StylePath(targetPath), utils.StylePath(absPath), uid, gid)

	// 3. Recursive Discovery
	inodes, err := scanInodes(ctx, absPath, targetPath)
	if err != nil {
		return &Error{
			Op:      "scan inodes",
			Path:    absPath,
			Tool:    "debugfs",
			Output:  err.Error(),
			BaseErr: err,
		}
	}

	// 4. Batch Update
	// We always attempt to chown the skeleton structure (/upper, /work, /work/work) if they exist.
	// This ensures the overlay is usable by the target user even if the target path is empty.
	skeletonInodes := []string{"/upper", "/work", "/work/work"}

	uniqueInodes := make(map[string]bool)
	for _, inode := range inodes {
		uniqueInodes[inode] = true
	}

	if len(uniqueInodes) == 0 {
		utils.PrintDebug("No inodes found to modify at %s", targetPath)
	}

	utils.PrintMessage("Updating %d inodes (uid=%d gid=%d) in %s",
		len(uniqueInodes), uid, gid, utils.StylePath(absPath))

	var cmds []string

	// A. Ensure skeleton paths are accessible
	for _, p := range skeletonInodes {
		cmds = append(cmds, fmt.Sprintf("sif %s uid %d", p, uid))
		cmds = append(cmds, fmt.Sprintf("sif %s gid %d", p, gid))
	}

	// B. Fix discovered inodes
	for inode := range uniqueInodes {
		cmds = append(cmds, fmt.Sprintf("sif <%s> uid %d", inode, uid))
		cmds = append(cmds, fmt.Sprintf("sif <%s> gid %d", inode, gid))
	}
	cmds = append(cmds, "quit")

	// C. Execute Batch
	script := strings.Join(cmds, "\n")
	cmd := exec.CommandContext(ctx, "debugfs", "-w", absPath)
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

	utils.PrintSuccess("Permissions updated for %s", utils.StylePath(absPath))
	return nil
}

// scanInodes implements a BFS walker to find all inodes under startPath.
func scanInodes(ctx context.Context, imgPath, startPath string) ([]string, error) {
	// A. Get Start Inode via 'stat'
	cmd := exec.CommandContext(ctx, "debugfs", "-R", fmt.Sprintf("stat %s", startPath), imgPath)
	utils.PrintDebug("[chown] debugfs -R stat %s %s", startPath, imgPath)
	out, err := cmd.Output()
	if err != nil {
		// This usually means the path doesn't exist in the overlay layer yet
		// (e.g. user hasn't written to it).
		utils.PrintDebug("Path %s not found in overlay layer: %v", startPath, err)
		return nil, nil // Return empty without error
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
	results := make(map[string]bool)
	results[startInode] = true

	count := 0

	for len(toVisit) > 0 {
		select {
		case <-ctx.Done():
			return nil, ctx.Err()
		default:
		}
		current := toVisit[0]
		toVisit = toVisit[1:]

		if visited[current] {
			continue
		}
		visited[current] = true

		// List contents: ls -l <INODE>
		cmdLs := exec.CommandContext(ctx, "debugfs", "-R", fmt.Sprintf("ls -l <%s>", current), imgPath)
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

			results[inode] = true

			// Check if directory (Mode starts with 4 or 04)
			// 40755 -> Directory (octal)
			if strings.HasPrefix(mode, "4") || strings.HasPrefix(mode, "04") {
				if !visited[inode] {
					toVisit = append(toVisit, inode)
				}
			}

			count++
			if count%1000 == 0 {
				fmt.Printf("\r    Scanning inodes... %s found", utils.StyleNumber(count))
			}
		}
	}
	fmt.Printf("\r%s\r", strings.Repeat(" ", 50)) // Clear line

	inodeList := make([]string, 0, len(results))
	for inode := range results {
		inodeList = append(inodeList, inode)
	}

	utils.PrintMessage("Finished scanning %s inodes", utils.StyleNumber(len(inodeList)))

	return inodeList, nil
}
