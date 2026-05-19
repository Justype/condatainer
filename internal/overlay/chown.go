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

	"github.com/Justype/condatainer/internal/logging"
)

// ChownRecursively changes the UID/GID of files inside an unmounted overlay image.
// It maps the provided internalPath (e.g., "/") to the overlay's "/upper" structure.
func ChownRecursively(ctx context.Context, imagePath string, uid, gid int, internalPath string) error {
	if err := checkDependencies([]string{"debugfs"}); err != nil {
		return err
	}

	absPath, err := filepath.Abs(imagePath)
	if err != nil {
		return fmt.Errorf("failed to resolve path %s: %w", imagePath, err)
	}

	targetPath := filepath.Join("/upper", strings.TrimPrefix(internalPath, "/"))
	log := logging.FromContext(ctx)

	log.Info(fmt.Sprintf("scanning %s inside %s (uid=%d gid=%d)", targetPath, absPath, uid, gid))

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

	skeletonInodes := []string{"/upper", "/work", "/work/work"}

	uniqueInodes := make(map[string]bool)
	for _, inode := range inodes {
		uniqueInodes[inode] = true
	}

	if len(uniqueInodes) == 0 {
		log.Debug(fmt.Sprintf("no inodes found to modify at %s", targetPath))
	}

	log.Info(fmt.Sprintf("updating %d inodes (uid=%d gid=%d) in %s", len(uniqueInodes), uid, gid, absPath))

	var cmds []string

	for _, p := range skeletonInodes {
		cmds = append(cmds, fmt.Sprintf("sif %s uid %d", p, uid))
		cmds = append(cmds, fmt.Sprintf("sif %s gid %d", p, gid))
	}

	for inode := range uniqueInodes {
		cmds = append(cmds, fmt.Sprintf("sif <%s> uid %d", inode, uid))
		cmds = append(cmds, fmt.Sprintf("sif <%s> gid %d", inode, gid))
	}
	cmds = append(cmds, "quit")

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

	log.Info(fmt.Sprintf("permissions updated for %s", absPath), "kind", "success")
	return nil
}

// scanInodes implements a BFS walker to find all inodes under startPath.
func scanInodes(ctx context.Context, imgPath, startPath string) ([]string, error) {
	log := logging.FromContext(ctx)

	cmd := exec.CommandContext(ctx, "debugfs", "-R", fmt.Sprintf("stat %s", startPath), imgPath)
	log.Debug(fmt.Sprintf("debugfs -R stat %s %s", startPath, imgPath))
	out, err := cmd.Output()
	if err != nil {
		log.Debug(fmt.Sprintf("path %s not found in overlay layer: %v", startPath, err))
		return nil, nil
	}

	reInode := regexp.MustCompile(`Inode:\s+(\d+)`)
	match := reInode.FindStringSubmatch(string(out))
	if match == nil {
		return nil, fmt.Errorf("could not parse inode for %s", startPath)
	}
	startInode := match[1]

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

		cmdLs := exec.CommandContext(ctx, "debugfs", "-R", fmt.Sprintf("ls -l <%s>", current), imgPath)
		outLs, _ := cmdLs.Output()

		scanner := bufio.NewScanner(bytes.NewReader(outLs))
		for scanner.Scan() {
			line := strings.TrimSpace(scanner.Text())
			parts := strings.Fields(line)

			if len(parts) < 6 {
				continue
			}

			inode := parts[0]
			mode := parts[1]
			name := parts[len(parts)-1]

			if name == "." || name == ".." {
				continue
			}

			if _, err := strconv.Atoi(inode); err != nil {
				continue
			}

			results[inode] = true

			if strings.HasPrefix(mode, "4") || strings.HasPrefix(mode, "04") {
				if !visited[inode] {
					toVisit = append(toVisit, inode)
				}
			}

			count++
			if count%1000 == 0 {
				log.Info(fmt.Sprintf("scanning inodes... %d found", count))
			}
		}
	}

	inodeList := make([]string, 0, len(results))
	for inode := range results {
		inodeList = append(inodeList, inode)
	}

	log.Info(fmt.Sprintf("finished scanning %d inodes", len(inodeList)))

	return inodeList, nil
}
