module github.com/Justype/condatainer

go 1.26.4

// spf13/cobra         for Arg Parse
// spf13/viper         for Config Management
// fatih/color         for ANSI Colors
// golang.org/x/term   for Terminal Detection

replace golang.org/x/crypto => github.com/Justype/x-crypto v0.50.0-hostbased.1

require (
	github.com/chzyer/readline v1.5.1
	github.com/fatih/color v1.19.0
	github.com/spf13/cobra v1.10.2
	github.com/spf13/pflag v1.0.10
	github.com/spf13/viper v1.21.0
	golang.org/x/crypto v0.52.0-hostbased.1
	golang.org/x/mod v0.36.0
	golang.org/x/sys v0.45.0
	golang.org/x/term v0.43.0
)

require (
	github.com/fsnotify/fsnotify v1.10.1 // indirect
	github.com/go-viper/mapstructure/v2 v2.5.0 // indirect
	github.com/inconshreveable/mousetrap v1.1.0 // indirect
	github.com/mattn/go-colorable v0.1.15 // indirect
	github.com/mattn/go-isatty v0.0.22 // indirect
	github.com/pelletier/go-toml/v2 v2.3.1 // indirect
	github.com/sagikazarmark/locafero v0.12.0 // indirect
	github.com/spf13/afero v1.15.0 // indirect
	github.com/spf13/cast v1.10.0 // indirect
	github.com/subosito/gotenv v1.6.0 // indirect
	go.yaml.in/yaml/v3 v3.0.4 // indirect
	golang.org/x/text v0.37.0 // indirect
	gopkg.in/check.v1 v1.0.0-20190902080502-41f04d3bba15 // indirect
)
