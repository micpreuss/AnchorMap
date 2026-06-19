# cli.R - argument parsing + usage text for the command-line entry points.
# Hand-rolled (no argparse/optparse dependency). parse_*() are pure + unit-testable; cli_*() are the
# Rscript-facing wrappers that print help / errors and exit.

engine_usage <- function() {
  paste(
    "AnchorMap engine - score which ontology domain each cluster anchors to.",
    "",
    "Usage:",
    "  anchor_map --config <config.yaml|name> [options]",
    "",
    "Required:",
    "  --config PATH        YAML config, or a bare shipped-config name (e.g. synthetic_rds).",
    "",
    "Input overrides (else taken from the config; paths are relative to the working dir):",
    "  --rds PATH           GenomicSEM ldsc() .rds artifact (Input C route).",
    "  --rg-long PATH       cluster x trait rg long-table TSV (Input A).",
    "  --trait-rg PATH      trait x trait LDSC --rg summary TSV (redundancy source).",
    "  --ontology PATH      ontology TSV.",
    "",
    "Output control:",
    "  --out-dir PATH       output directory (default: out_dir from the config).",
    "  --run-label NAME     run label (logged; fallback out_dir = results/<NAME>).",
    "",
    "Engine:",
    "  --threads N          worker/thread count (default: 1).",
    "  --z-vector a,b,c     reliability-threshold sweep (default: config z_vector, e.g. 3,4,5).",
    "  -h, --help           show this help and exit.",
    "",
    "Example:",
    "  anchor_map --config example_anthro.yaml --out-dir results/run1 --threads 4",
    "",
    sep = "\n")
}

# Parse the engine CLI. Returns a named list; errors on unknown flags. `--help` short-circuits.
parse_engine_args <- function(a) {
  out <- list(config = NULL, threads = 1L, rds = NULL, z_vector = NULL, out_dir = NULL,
              run_label = NULL, rg_long = NULL, trait_rg = NULL, ontology = NULL, help = FALSE)
  i <- 1L
  need <- function(i) { if (i > length(a)) stop(sprintf("option %s needs a value", a[i - 1L])); a[i] }
  while (i <= length(a)) {
    switch(a[i],
      "--config"    = { out$config    <- need(i + 1L); i <- i + 2L },
      "--threads"   = { out$threads   <- as.integer(need(i + 1L)); i <- i + 2L },
      "--rds"       = { out$rds       <- need(i + 1L); i <- i + 2L },
      "--rg-long"   = { out$rg_long   <- need(i + 1L); i <- i + 2L },
      "--trait-rg"  = { out$trait_rg  <- need(i + 1L); i <- i + 2L },
      "--ontology"  = { out$ontology  <- need(i + 1L); i <- i + 2L },
      "--out-dir"   = { out$out_dir   <- need(i + 1L); i <- i + 2L },
      "--run-label" = { out$run_label <- need(i + 1L); i <- i + 2L },
      "--z-vector"  = { out$z_vector  <- as.numeric(strsplit(need(i + 1L), "[, ]+")[[1]]); i <- i + 2L },
      "--help"      = { out$help <- TRUE; i <- i + 1L },
      "-h"          = { out$help <- TRUE; i <- i + 1L },
      stop(sprintf("unknown option: %s", a[i])))
  }
  out
}

# Rscript-facing engine entry. Prints help / usage and exits non-zero on bad input.
cli_anchor_map <- function(a = commandArgs(TRUE)) {
  args <- tryCatch(parse_engine_args(a),
                   error = function(e) { message("error: ", conditionMessage(e)); cat(engine_usage()); quit(status = 2L) })
  if (isTRUE(args$help)) { cat(engine_usage()); quit(status = 0L) }
  if (is.null(args$config)) { message("error: --config is required\n"); cat(engine_usage()); quit(status = 2L) }
  run_anchormap(args$config, args$threads, args$rds, args$z_vector, args$out_dir,
                args$run_label, args$rg_long, args$trait_rg, args$ontology)
}

plots_usage <- function() {
  paste(
    "AnchorMap figures - render anchoring + specificity figures from scored TSVs.",
    "",
    "Usage:",
    "  plot_anchors --config <plots.yaml|name> [options]",
    "",
    "Required:",
    "  --config PATH        plot YAML config, or a bare shipped-config name.",
    "",
    "Options:",
    "  --out-dir PATH       output directory for figures (default: out_dir from the config).",
    "  --q-sig N            significance threshold for rings/masks (default: 0.05).",
    "  --rg-floor N         min |pooled_rg| for a specificity cell (default: 0.10).",
    "  --min-clusters N     min clusters scoring a category for a stable z (default: 5).",
    "  -h, --help           show this help and exit.",
    "",
    "Example:",
    "  plot_anchors --config example_plots.yaml --out-dir results/run1/figures",
    "",
    sep = "\n")
}

# Parse the figures CLI. Returns a named list; errors on unknown flags. `--help` short-circuits.
parse_plots_args <- function(a) {
  out <- list(config = NULL, q_sig = NULL, rg_floor = NULL, min_clusters = NULL,
              out_dir = NULL, help = FALSE)
  i <- 1L
  need <- function(i) { if (i > length(a)) stop(sprintf("option %s needs a value", a[i - 1L])); a[i] }
  while (i <= length(a)) {
    switch(a[i],
      "--config"       = { out$config       <- need(i + 1L); i <- i + 2L },
      "--q-sig"        = { out$q_sig        <- as.numeric(need(i + 1L)); i <- i + 2L },
      "--rg-floor"     = { out$rg_floor     <- as.numeric(need(i + 1L)); i <- i + 2L },
      "--min-clusters" = { out$min_clusters <- as.integer(need(i + 1L)); i <- i + 2L },
      "--out-dir"      = { out$out_dir      <- need(i + 1L); i <- i + 2L },
      "--help"         = { out$help <- TRUE; i <- i + 1L },
      "-h"             = { out$help <- TRUE; i <- i + 1L },
      stop(sprintf("unknown option: %s", a[i])))
  }
  out
}

# Rscript-facing figures entry.
cli_plot_anchors <- function(a = commandArgs(TRUE)) {
  args <- tryCatch(parse_plots_args(a),
                   error = function(e) { message("error: ", conditionMessage(e)); cat(plots_usage()); quit(status = 2L) })
  if (isTRUE(args$help)) { cat(plots_usage()); quit(status = 0L) }
  if (is.null(args$config)) { message("error: --config is required\n"); cat(plots_usage()); quit(status = 2L) }
  run_plots(args$config, args$q_sig, args$rg_floor, args$min_clusters, args$out_dir)
}
