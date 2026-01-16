package cmd

import (
	"github.com/Justype/condatainer/internal/config"
	"github.com/Justype/condatainer/internal/utils"
	"github.com/spf13/cobra"
)

// Variables to hold flag values
var (
	createName      string
	createPrefix    string
	createFile      string
	createBaseImage string
	createSource    string
	createTempSize  string
	
	// Compression flags
	compZstd       bool
	compZstdFast   bool
	compZstdMedium bool
	compZstdHigh   bool
	compGzip       bool
	compLz4        bool
)

var createCmd = &cobra.Command{
	Use:     "create [packages...]",
	Aliases: []string{"install", "i"},
	Short:   "Create a new SquashFS overlay",
	Long:    `Create a new SquashFS overlay using available build scripts or Conda packages.`,
	Run: func(cmd *cobra.Command, args []string) {
		// 1. Validation Logic (Ported from Python)
		if len(args) == 0 && createFile == "" && createSource == "" {
			utils.PrintError("At least one of [packages], --file, or --source must be provided.")
			return
		}

		if len(args) > 0 && createPrefix != "" {
			utils.PrintError("Packages cannot be used with --prefix.")
			return
		}

		// 2. Handle Compression Config
		// (We update the global config based on flags)
		if compLz4 {
			config.Global.CompressArgs = "-comp lz4"
		} else if compZstd {
			config.Global.CompressArgs = "-comp zstd -Xcompression-level 14"
		} else if compZstdFast {
			config.Global.CompressArgs = "-comp zstd -Xcompression-level 3"
		} else if compZstdMedium {
			config.Global.CompressArgs = "-comp zstd -Xcompression-level 8"
		} else if compZstdHigh {
			config.Global.CompressArgs = "-comp zstd -Xcompression-level 19"
		} else if compGzip {
			config.Global.CompressArgs = "-comp gzip"
		}

		// 3. Temporary Debug Output (Until we implement the builder)
		utils.PrintMessage("Starting creation process...")
		utils.PrintDebug("Compression: %s", utils.StyleName(config.Global.CompressArgs))
		if len(args) > 0 {
			utils.PrintDebug("Target Packages: %v", args)
		}
		
		// TODO: Initialize internal/builder and run the build logic
	},
}

func init() {
	rootCmd.AddCommand(createCmd)

	// Register Flags
	f := createCmd.Flags()
	f.StringVarP(&createName, "name", "n", "", "Custom name for the resulting overlay file")
	f.StringVarP(&createPrefix, "prefix", "p", "", "Custom prefix path for the overlay file")
	f.StringVarP(&createFile, "file", "f", "", "Path to definition file (.yaml, .sh, .def)")
	f.StringVarP(&createBaseImage, "base-image", "b", "", "Base image to use instead of default")
	f.StringVarP(&createSource, "source", "s", "", "Remote source URI (e.g., docker://ubuntu:22.04)")
	f.StringVar(&createTempSize, "temp-size", "20G", "Size of temporary overlay")

	// Compression Flags
	f.BoolVar(&compZstdFast, "zstd-fast", false, "Use zstd compression level 3")
	f.BoolVar(&compZstdMedium, "zstd-medium", false, "Use zstd compression level 8")
	f.BoolVar(&compZstd, "zstd", false, "Use zstd compression level 14")
	f.BoolVar(&compZstdHigh, "zstd-high", false, "Use zstd compression level 19")
	f.BoolVar(&compGzip, "gzip", false, "Use gzip compression")
	f.BoolVar(&compLz4, "lz4", false, "Use LZ4 compression (default)")
}