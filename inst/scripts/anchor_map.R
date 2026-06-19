#!/usr/bin/env Rscript
# Thin CLI wrapper for the AnchorMap engine. Install the package, then:
#   Rscript -e 'anchormap:::cli_anchor_map()' --config <yaml> ...
# or run this file:  Rscript <path>/anchor_map.R --config <yaml> ...
# or use the `anchor_map` shim on PATH (in the Docker image).
suppressPackageStartupMessages(library(anchormap))
anchormap:::cli_anchor_map(commandArgs(TRUE))
