package helper

import "github.com/Justype/condatainer/internal/helperhistory"

// UsedCombination is one project-shared overlay combination, recorded and
// read through internal/helperhistory. Re-exported here so a helper-package
// caller (the CLI, the dashboard) never has to import that package
// directly; internal/project imports it on its own, since importing
// internal/helper from there would cycle back through
// internal/helper's own dependency on internal/project.
type UsedCombination = helperhistory.UsedCombination

// RecordUsed records that name was run at location (project-root-relative,
// slash-separated) with overlays, in root's shared history. See
// internal/helperhistory.RecordUsed.
func RecordUsed(root, name, location string, overlays []string) error {
	return helperhistory.RecordUsed(root, name, location, overlays)
}

// ListUsed returns every combination recorded for name at location, newest
// first. See internal/helperhistory.ListUsed.
func ListUsed(root, name, location string) ([]*UsedCombination, error) {
	return helperhistory.ListUsed(root, name, location)
}

// ListAll groups every combination in root's shared history by helper, then
// location. See internal/helperhistory.ListAll.
func ListAll(root string) (map[string]map[string][]*UsedCombination, error) {
	return helperhistory.ListAll(root)
}
