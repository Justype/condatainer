package conda

import (
	"reflect"
	"strings"
	"testing"
)

func TestSetChannelsPreservesOtherConfig(t *testing.T) {
	in := []byte("# environment policy\nchannel_priority: strict\nchannels:\n  - defaults\nssl_verify: /etc/ssl/cert.pem\n")
	out, err := SetChannels(in, []string{"conda-forge", "bioconda"})
	if err != nil {
		t.Fatal(err)
	}
	channels, err := ParseChannels(out)
	if err != nil {
		t.Fatal(err)
	}
	if want := []string{"conda-forge", "bioconda"}; !reflect.DeepEqual(channels, want) {
		t.Fatalf("channels = %v, want %v", channels, want)
	}
	text := string(out)
	for _, retained := range []string{"# environment policy", "channel_priority: strict", "ssl_verify: /etc/ssl/cert.pem"} {
		if !strings.Contains(text, retained) {
			t.Errorf("updated config lost %q:\n%s", retained, text)
		}
	}
}

func TestParseChannelsRejectsNonList(t *testing.T) {
	if _, err := ParseChannels([]byte("channels: conda-forge\n")); err == nil {
		t.Fatal("expected scalar channels to fail")
	}
}
