package registry

import (
	"fmt"
	"os"
	"os/user"
	"path/filepath"
	"strconv"
	"syscall"
)

// Reach is a class of people, from the narrowest to the widest.
type Reach int

const (
	// ReachOwner is the file's owner alone.
	ReachOwner Reach = iota
	// ReachGroup is the owner and the file's group.
	ReachGroup
	// ReachEveryone is every user on the system.
	ReachEveryone
)

// Access is who can read a credential file and who can replace it, from the
// file's and its directory's mode bits. ACLs are not consulted.
type Access struct {
	Readers Reach
	Writers Reach
	Group   string // the file's group, named in a message about the group
}

// reader names who can read the file, for a message.
func (a Access) reader() string { return a.describe(a.Readers) }

func (a Access) describe(aud Reach) string {
	switch aud {
	case ReachGroup:
		return "group " + a.Group
	case ReachEveryone:
		return "everyone"
	}
	return "you"
}

// accessOf reads who can reach path. A class reads when the directory lets it in
// and the file lets it read. It replaces when the directory lets it in and
// either the file is writable by it or the directory is, unless the directory is
// sticky, where deleting another user's file is refused.
func accessOf(path string) (Access, error) {
	info, err := os.Stat(path)
	if err != nil {
		return Access{}, err
	}
	dir, err := os.Stat(filepath.Dir(path))
	if err != nil {
		return Access{}, err
	}
	fm, dm := info.Mode().Perm(), dir.Mode().Perm()
	sticky := dir.Mode()&os.ModeSticky != 0

	access := Access{}
	switch {
	case dm&0o001 != 0 && fm&0o004 != 0:
		access.Readers = ReachEveryone
	case dm&0o010 != 0 && fm&0o040 != 0:
		access.Readers = ReachGroup
	}
	switch {
	case dm&0o001 != 0 && (fm&0o002 != 0 || dm&0o002 != 0 && !sticky):
		access.Writers = ReachEveryone
	case dm&0o010 != 0 && (fm&0o020 != 0 || dm&0o020 != 0 && !sticky):
		access.Writers = ReachGroup
	}
	if stat, ok := info.Sys().(*syscall.Stat_t); ok {
		access.Group = strconv.FormatUint(uint64(stat.Gid), 10)
		if g, err := user.LookupGroupId(access.Group); err == nil {
			access.Group = g.Name
		}
	}
	return access, nil
}

// Finding is one thing to tell the person who saved a credential.
type Finding struct {
	Warn bool // false is a note
	Text string
}

// Findings reports what is wrong with who can reach a layer's credential file. The
// user layer is personal, so anyone else reading it is a warning; a shared layer
// is meant to be read by others, so its readers are only shown. Anyone other than
// the owner being able to replace the file is always reported.
func Findings(layer, path string) []Finding {
	access, err := accessOf(path)
	if err != nil {
		return nil
	}
	var out []Finding
	if layer == "user" && access.Readers != ReachOwner {
		out = append(out, Finding{Warn: true, Text: fmt.Sprintf(
			"%s is readable by %s. Run: chmod 600 %s", path, access.reader(), path)})
	}
	switch access.Writers {
	case ReachEveryone:
		out = append(out, Finding{Warn: true, Text: fmt.Sprintf(
			"anyone can replace the credential in %s. Run: chmod go-w %s, and check that %s is not writable by others",
			path, path, filepath.Dir(path))})
	case ReachGroup:
		out = append(out, Finding{Warn: layer == "user", Text: fmt.Sprintf(
			"members of %s can replace the credential in %s", access.describe(ReachGroup), path)})
	}
	return out
}

// ReadableBy names who can read the file saved for a layer: "you",
// "group <name>" or "everyone". It is empty when the file cannot be inspected.
func ReadableBy(path string) string {
	access, err := accessOf(path)
	if err != nil {
		return ""
	}
	return access.reader()
}

// LayerFile is the credential file of a config layer, wherever it may or may not
// exist yet.
func LayerFile(layer string) (string, error) { return credentialFilePath(layer) }
