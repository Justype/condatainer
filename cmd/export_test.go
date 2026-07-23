package cmd

import "testing"

func TestReorderChannelBlock(t *testing.T) {
	priority := []string{"conda-forge", "bioconda"}

	t.Run("reorders to priority order", func(t *testing.T) {
		in := "name: base\nchannels:\n  - bioconda\n  - conda-forge\ndependencies:\n  - fastp=1.3.6\n"
		want := "name: base\nchannels:\n  - conda-forge\n  - bioconda\ndependencies:\n  - fastp=1.3.6\n"
		if got := string(reorderChannelBlock([]byte(in), priority)); got != want {
			t.Errorf("got:\n%s\nwant:\n%s", got, want)
		}
	})

	t.Run("channel not in condarc is prepended (matches conda insert(0))", func(t *testing.T) {
		in := "channels:\n  - bioconda\n  - pytorch\n  - conda-forge\ndependencies:\n"
		want := "channels:\n  - pytorch\n  - conda-forge\n  - bioconda\ndependencies:\n"
		if got := string(reorderChannelBlock([]byte(in), priority)); got != want {
			t.Errorf("got:\n%s\nwant:\n%s", got, want)
		}
	})

	t.Run("no channels block is unchanged", func(t *testing.T) {
		in := "@EXPLICIT\nhttps://conda.anaconda.org/conda-forge/linux-64/x.conda\n"
		if got := string(reorderChannelBlock([]byte(in), priority)); got != in {
			t.Errorf("expected unchanged, got:\n%s", got)
		}
	})

	t.Run("single channel is unchanged", func(t *testing.T) {
		in := "channels:\n  - conda-forge\ndependencies:\n"
		if got := string(reorderChannelBlock([]byte(in), priority)); got != in {
			t.Errorf("expected unchanged, got:\n%s", got)
		}
	})
}
