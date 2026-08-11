package key

import (
	"fmt"

	"github.com/Justype/condatainer/catalog"
	"github.com/Justype/condatainer/internal/artifact/meta"
)

// CondaExplicitV1 answers: "Which exact Conda package set is installed?"
//
// Its preimage is explicit.txt exactly as stored, byte for byte. It performs no
// parsing, normalization, sorting, or model encoding.
func deriveCondaExplicitV1(m meta.Manifest, explicit []byte) (Value, error) {
	if m.Type != catalog.TypeApp {
		return Value{}, fmt.Errorf("%s cannot derive type %s", CondaExplicitV1, m.Type)
	}
	if len(m.Dependencies) > 0 {
		return Value{}, fmt.Errorf("%s cannot have dependencies", CondaExplicitV1)
	}
	return value(CondaExplicitV1, explicit), nil
}
