#!/usr/bin/env Rscript
# Thin CLI wrapper for the AnchorMap figures. Install the package, then:
#   Rscript -e 'anchormap:::cli_plot_anchors()' --config <plots.yaml> ...
# or run this file:  Rscript <path>/plot_anchors.R --config <plots.yaml> ...
# or use the `plot_anchors` shim on PATH (in the Docker image).
suppressPackageStartupMessages(library(anchormap))
anchormap:::cli_plot_anchors(commandArgs(TRUE))
