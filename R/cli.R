# cli.R - argument parsing + usage text for the command-line entry points.
# Built on optparse: the option lists are declarative, parse_*() return a normalized named list
# (pure + unit-testable), and cli_*() are the Rscript-facing wrappers that print help / errors and
# exit. All optparse calls are fully qualified (optparse is in Imports), so no NAMESPACE import.

# ---- engine -----------------------------------------------------------------

engine_parser <- function() {
  optparse::OptionParser(
    usage = "anchor_map --config <config.yaml|name> [options]",
    description = paste(
      "",
      "AnchorMap engine - score which ontology domain each cluster anchors to.",
      "--config takes a YAML path or a bare shipped-config name (e.g. synthetic_rds).",
      sep = "\n"),
    epilogue = paste(
      "Input overrides default to the config; paths are relative to the working dir.",
      "Example:",
      "  anchor_map --config example_anthro.yaml --out-dir results/run1 --threads 4",
      "", sep = "\n"),
    option_list = list(
      optparse::make_option("--config", type = "character", default = NULL,
        help = "YAML config, or a bare shipped-config name (required)."),
      optparse::make_option("--rds", type = "character", default = NULL,
        help = "GenomicSEM ldsc() .rds artifact (Input C route)."),
      optparse::make_option("--rg-long", type = "character", default = NULL,
        help = "cluster x trait rg long-table TSV (Input A)."),
      optparse::make_option("--trait-rg", type = "character", default = NULL,
        help = "trait x trait LDSC --rg summary TSV (redundancy source)."),
      optparse::make_option("--ontology", type = "character", default = NULL,
        help = "ontology TSV."),
      optparse::make_option("--out-dir", type = "character", default = NULL,
        help = "output directory (default: out_dir from the config)."),
      optparse::make_option("--run-label", type = "character", default = NULL,
        help = "run label (logged; fallback out_dir = results/<NAME>)."),
      optparse::make_option("--threads", type = "integer", default = 1L,
        help = "worker/thread count [default %default]."),
      optparse::make_option("--z-vector", type = "character", default = NULL,
        help = "reliability-threshold sweep, comma-separated (default: config z_vector, e.g. 3,4,5).")))
}

engine_usage <- function() {
  paste(utils::capture.output(optparse::print_help(engine_parser())), collapse = "\n")
}

# Parse the engine CLI. Returns a normalized named list (z_vector split to numeric); errors on
# unknown flags or missing values. `--help` / `-h` set $help without exiting.
parse_engine_args <- function(a) {
  o <- optparse::parse_args(engine_parser(), args = a,
                            convert_hyphens_to_underscores = TRUE,
                            print_help_and_exit = FALSE)
  if (!is.null(o$z_vector)) o$z_vector <- as.numeric(strsplit(o$z_vector, "[, ]+")[[1]])
  o
}

# Rscript-facing engine entry. Prints help / usage and exits non-zero on bad input.
cli_anchor_map <- function(a = commandArgs(TRUE)) {
  args <- tryCatch(parse_engine_args(a),
                   error = function(e) { message("error: ", conditionMessage(e)); quit(status = 2L) })
  if (isTRUE(args$help)) { cat(engine_usage()); quit(status = 0L) }
  if (is.null(args$config)) { message("error: --config is required"); cat(engine_usage()); quit(status = 2L) }
  run_anchormap(args$config, args$threads, args$rds, args$z_vector, args$out_dir,
                args$run_label, args$rg_long, args$trait_rg, args$ontology)
}

# ---- figures ----------------------------------------------------------------

plots_parser <- function() {
  optparse::OptionParser(
    usage = "plot_anchors --config <plots.yaml|name> [options]",
    description = paste(
      "",
      "AnchorMap figures - render anchoring + specificity figures from scored TSVs.",
      "--config takes a plot YAML path or a bare shipped-config name.",
      sep = "\n"),
    epilogue = paste(
      "Example:",
      "  plot_anchors --config example_plots.yaml --out-dir results/run1/figures",
      "", sep = "\n"),
    option_list = list(
      optparse::make_option("--config", type = "character", default = NULL,
        help = "plot YAML config, or a bare shipped-config name (required)."),
      optparse::make_option("--out-dir", type = "character", default = NULL,
        help = "output directory for figures (default: out_dir from the config)."),
      optparse::make_option("--in-dir", type = "character", default = NULL,
        help = "read scored TSVs from this dir (by basename) instead of the config's track paths."),
      optparse::make_option("--q-sig", type = "double", default = NULL,
        help = "significance threshold for rings/masks [default 0.05]."),
      optparse::make_option("--rg-floor", type = "double", default = NULL,
        help = "min |pooled_rg| for a specificity cell [default 0.10]."),
      optparse::make_option("--min-clusters", type = "integer", default = NULL,
        help = "min clusters scoring a category for a stable z [default 5].")))
}

plots_usage <- function() {
  paste(utils::capture.output(optparse::print_help(plots_parser())), collapse = "\n")
}

# Parse the figures CLI. Returns a named list; errors on unknown flags. `--help` sets $help.
parse_plots_args <- function(a) {
  optparse::parse_args(plots_parser(), args = a,
                       convert_hyphens_to_underscores = TRUE,
                       print_help_and_exit = FALSE)
}

# Rscript-facing figures entry.
cli_plot_anchors <- function(a = commandArgs(TRUE)) {
  args <- tryCatch(parse_plots_args(a),
                   error = function(e) { message("error: ", conditionMessage(e)); quit(status = 2L) })
  if (isTRUE(args$help)) { cat(plots_usage()); quit(status = 0L) }
  if (is.null(args$config)) { message("error: --config is required"); cat(plots_usage()); quit(status = 2L) }
  run_plots(args$config, args$q_sig, args$rg_floor, args$min_clusters, args$out_dir, args$in_dir)
}
