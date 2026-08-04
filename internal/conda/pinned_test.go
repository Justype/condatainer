package conda

import (
	"reflect"
	"testing"
)

func TestPinAddReplaceAndRemove(t *testing.T) {
	path := t.TempDir() + "/conda-meta/pinned"
	if err := AddPin(path, "numpy ==2.1.0"); err != nil {
		t.Fatal(err)
	}
	if err := AddPin(path, "scipy 1.14.*"); err != nil {
		t.Fatal(err)
	}
	if err := AddPin(path, "numpy ==2.2.0"); err != nil {
		t.Fatal(err)
	}
	pins, err := ReadPins(path)
	if err != nil {
		t.Fatal(err)
	}
	if want := []string{"scipy 1.14.*", "numpy ==2.2.0"}; !reflect.DeepEqual(pins, want) {
		t.Fatalf("pins = %v, want %v", pins, want)
	}
	if err := RemovePin(path, "numpy"); err != nil {
		t.Fatal(err)
	}
	pins, err = ReadPins(path)
	if err != nil {
		t.Fatal(err)
	}
	if want := []string{"scipy 1.14.*"}; !reflect.DeepEqual(pins, want) {
		t.Fatalf("pins = %v, want %v", pins, want)
	}
}

func TestRemovePinUsesLiteralPackageName(t *testing.T) {
	path := t.TempDir() + "/conda-meta/pinned"
	if err := AddPin(path, "r-base ==4.4.0"); err != nil {
		t.Fatal(err)
	}
	if err := RemovePin(path, "r.base"); err == nil {
		t.Fatal("regex-like package name unexpectedly removed a different pin")
	}
}
