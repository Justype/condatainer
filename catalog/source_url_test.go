package catalog

import (
	"errors"
	"strings"
	"testing"
)

func TestParseSourceTakesNameThenURLs(t *testing.T) {
	rec := parse(t, "recipes/grch38/idx", strings.Join([]string{
		"#TYPE:data",
		"#SOURCE:gtf https://example.invalid/a.gtf.gz",
		"#SOURCE:genome https://example.invalid/g.fa.gz",
		"echo build",
	}, "\n"))

	if len(rec.Sources) != 2 {
		t.Fatalf("parsed %d sources, want 2: %+v", len(rec.Sources), rec.Sources)
	}
	if got := rec.Sources[0]; got.Name != "gtf" || got.URL != "https://example.invalid/a.gtf.gz" {
		t.Errorf("first source = %+v", got)
	}
	if got := rec.Sources[1]; got.Name != "genome" || got.URL != "https://example.invalid/g.fa.gz" {
		t.Errorf("second source = %+v", got)
	}
}

func TestMalformedSourceIsSkippedThenReported(t *testing.T) {
	for _, tc := range []struct{ name, line string }{
		{"no url", "#SOURCE:gtf"},
		{"two urls", "#SOURCE:gtf https://example.invalid/a https://mirror.invalid/a"},
		{"name with a dash", "#SOURCE:my-gtf https://example.invalid/a"},
		{"name starting with a digit", "#SOURCE:1gtf https://example.invalid/a"},
		{"empty", "#SOURCE:"},
	} {
		t.Run(tc.name, func(t *testing.T) {
			rec := parse(t, "recipes/grch38/idx", "#TYPE:data\n"+tc.line+"\necho build\n")
			if len(rec.Sources) != 0 {
				t.Fatalf("parsing kept %+v; a malformed line must be skipped", rec.Sources)
			}
			err := rec.Validate()
			if err == nil || !strings.Contains(err.Error(), "#SOURCE:") {
				t.Errorf("Validate() = %v, want it to name the bad declaration", err)
			}
		})
	}
}

func TestDuplicateSourceNameIsRejected(t *testing.T) {
	rec := parse(t, "recipes/grch38/idx", strings.Join([]string{
		"#TYPE:data",
		"#SOURCE:gtf https://example.invalid/a",
		"#SOURCE:gtf https://example.invalid/b",
		"echo build",
	}, "\n"))

	err := rec.Validate()
	if err == nil || !strings.Contains(err.Error(), "twice") {
		t.Errorf("Validate() = %v, want a duplicate-name rejection", err)
	}
}

func TestExpandSubstitutesSourceURLs(t *testing.T) {
	rec := parse(t, "recipes/grch38/gtf-gencode", strings.Join([]string{
		"#TARGET:grch38/gtf-gencode/{gencode_version}",
		"#PH:gencode_version:49,48",
		"#SOURCE:gtf https://example.invalid/release_{gencode_version}/v{gencode_version}.gtf.gz",
		"echo build",
	}, "\n"))

	out, err := Expand(rec, map[string]string{"gencode_version": "49"})
	if err != nil {
		t.Fatal(err)
	}
	want := "https://example.invalid/release_49/v49.gtf.gz"
	if got := out.Sources[0].URL; got != want {
		t.Errorf("expanded URL = %q, want %q", got, want)
	}
	if out.Sources[0].Name != "gtf" {
		t.Errorf("name = %q, want it unchanged by expansion", out.Sources[0].Name)
	}
	// The template itself must keep its tokens: it is what the artifact embeds.
	if got := rec.Sources[0].URL; !strings.Contains(got, "{gencode_version}") {
		t.Errorf("template URL = %q, want its tokens intact", got)
	}
}

// #SOURCE: is a comment, so it is stripped from the recipe preimage like every
// other annotation. Changing a URL must not move a key by itself — only the
// bytes it served can, through the src= lines.
func TestSourceLineLeavesTheRecipePreimage(t *testing.T) {
	const body = "#TYPE:data\n#SOURCE:gtf %s\necho build\n"
	a := StripComments([]byte(strings.ReplaceAll(body, "%s", "https://example.invalid/a")))
	b := StripComments([]byte(strings.ReplaceAll(body, "%s", "https://mirror.invalid/a")))
	if string(a) != string(b) {
		t.Errorf("a changed #SOURCE: URL moved the recipe preimage:\n%s\n%s", a, b)
	}
}

func TestParseSourceAskTakesAPromptNotURLs(t *testing.T) {
	rec := parse(t, "recipes/cellranger/9.0.1", strings.Join([]string{
		"#TYPE:app",
		"#INPUT:Accept the licence, then enter y",
		"#SOURCE:crx ask:10x links expire after one day. Paste the tar.gz link",
		"#SOURCE:gtf https://example.invalid/a.gtf",
		"echo build",
	}, "\n"))

	if len(rec.Sources) != 2 {
		t.Fatalf("parsed %d sources, want 2: %+v", len(rec.Sources), rec.Sources)
	}
	if got := rec.Sources[0]; got.Name != "crx" || got.URL != "" ||
		got.Prompt != "10x links expire after one day. Paste the tar.gz link" {
		t.Errorf("ask source = %+v, want the whole prompt and no URLs", got)
	}
	if got := rec.Sources[1]; got.Prompt != "" || got.URL == "" {
		t.Errorf("literal source = %+v, want a URL and no prompt", got)
	}
}

// #INPUT: answers reach the recipe and ask: answers do not, so the two must stay
// distinguishable by position: recipe prompts first, then asks in order.
func TestPromptsListInputThenAsks(t *testing.T) {
	rec := parse(t, "recipes/cellranger/9.0.1", strings.Join([]string{
		"#TYPE:app",
		"#SOURCE:crx ask:first ask",
		"#INPUT:the recipe's own question",
		"#SOURCE:idx ask:second ask",
		"echo build",
	}, "\n"))

	want := []string{"the recipe's own question", "first ask", "second ask"}
	if got := rec.Prompts(); strings.Join(got, "|") != strings.Join(want, "|") {
		t.Errorf("Prompts() = %v, want %v", got, want)
	}
	if len(rec.Inputs) != 1 {
		t.Errorf("Inputs = %v, want only the recipe's own", rec.Inputs)
	}
}

func TestAskWithNoPromptIsMalformed(t *testing.T) {
	rec := parse(t, "recipes/cellranger/9.0.1", "#TYPE:app\n#SOURCE:crx ask:\necho build\n")
	if len(rec.Sources) != 0 {
		t.Fatalf("parsing kept %+v; an ask with no prompt must be skipped", rec.Sources)
	}
	if err := rec.Validate(); err == nil || !strings.Contains(err.Error(), "#SOURCE:") {
		t.Errorf("Validate() = %v, want it to name the bad declaration", err)
	}
}

func TestExpandSubstitutesTheAskPrompt(t *testing.T) {
	rec := parse(t, "recipes/grch38/x", strings.Join([]string{
		"#TARGET:grch38/x/{v}",
		"#PH:v:1,2",
		"#SOURCE:crx ask:paste the link for release {v}",
		"echo build",
	}, "\n"))
	out, err := Expand(rec, map[string]string{"v": "2"})
	if err != nil {
		t.Fatal(err)
	}
	if got := out.Sources[0].Prompt; got != "paste the link for release 2" {
		t.Errorf("prompt = %q, want the placeholder substituted", got)
	}
}

// A definition is built by apptainer, not run as a script, so nothing would
// fetch a #SOURCE: or read an #INPUT: answer. Declaring either is refused
// rather than silently ignored.
func TestDefinitionRefusesSourceAndInput(t *testing.T) {
	for _, tc := range []struct{ name, header, key string }{
		{"source", "#SOURCE:gtf https://example.invalid/a", "#SOURCE"},
		{"ask source", "#SOURCE:crx ask:paste the link", "#SOURCE"},
		{"malformed source", "#SOURCE:gtf", "#SOURCE"},
		{"input", "#INPUT:accept the licence", "#INPUT"},
	} {
		t.Run(tc.name, func(t *testing.T) {
			rec := parse(t, "recipes/ubuntu24/thing.def",
				"Bootstrap: docker\nFrom: ubuntu:24.04\n"+tc.header+"\n")
			err := rec.Validate()
			if err == nil || !strings.Contains(err.Error(), tc.key+":") ||
				!strings.Contains(err.Error(), "definition") {
				t.Errorf("Validate() = %v, want a refusal naming %s on a definition", err, tc.key)
			}
		})
	}

	// The same headers stay legal in a script recipe.
	script := parse(t, "recipes/grch38/idx",
		"#TYPE:data\n#SOURCE:gtf https://example.invalid/a\n#INPUT:go on\necho build\n")
	if err := script.Validate(); err != nil {
		t.Errorf("a script recipe was refused: %v", err)
	}
}

func TestSourcePerArchitecture(t *testing.T) {
	rec := parse(t, "recipes/tool/1", strings.Join([]string{
		"#SOURCE:bin amd64 https://example.invalid/tool-x64",
		"#SOURCE:bin arm64 https://example.invalid/tool-a64",
		"#SOURCE:pkg arm64 ask:paste the arm64 link",
	}, "\n")+"\n")
	got, err := rec.SourcesFor("arm64")
	if err != nil || len(got) != 2 || got[0].URL != "https://example.invalid/tool-a64" || got[1].Prompt != "paste the arm64 link" {
		t.Fatalf("SourcesFor(arm64) = %+v, %v", got, err)
	}
	if _, err := rec.SourcesFor("amd64"); !errors.Is(err, ErrNoSourceForArch) {
		t.Errorf("SourcesFor(amd64) error = %v, want ErrNoSourceForArch", err)
	}
	for _, text := range []string{
		"#SOURCE:bin amd64 https://a\n#SOURCE:bin amd64 https://b\n",
		"#SOURCE:bin https://a\n#SOURCE:bin arm64 https://b\n",
		"#SOURCE:bin x86_64 https://a\n",
	} {
		if err := parse(t, "recipes/tool/1", text).Validate(); !errors.Is(err, ErrInvalidRecipe) {
			t.Errorf("Validate(%q) = %v, want ErrInvalidRecipe", text, err)
		}
	}
}
