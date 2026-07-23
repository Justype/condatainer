package overlay

import (
	"reflect"
	"testing"
)

func TestParseCondarcChannels(t *testing.T) {
	cases := []struct {
		name string
		in   string
		want []string
	}{
		{
			name: "channels then other key",
			in:   "channels:\n  - conda-forge\n  - bioconda\nchannel_priority: strict\n",
			want: []string{"conda-forge", "bioconda"},
		},
		{
			name: "quoted and blank lines",
			in:   "channel_priority: strict\nchannels:\n  - \"conda-forge\"\n\n  - bioconda\n",
			want: []string{"conda-forge", "bioconda"},
		},
		{
			name: "no channels key",
			in:   "channel_priority: strict\n",
			want: nil,
		},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			if got := parseCondarcChannels(c.in); !reflect.DeepEqual(got, c.want) {
				t.Errorf("got %v, want %v", got, c.want)
			}
		})
	}
}
