#!/usr/bin/env Rscript
# plot_anchors.R — AnchorMap Phase 4 figures CLI.
#
# Renders the cluster-anchoring + cross-cluster specificity figures from the scored TSVs. Reads a
# plot YAML config (out_dir + tracks:[{name, level, scores, labels}]) and writes, into out_dir:
#   anchor_lollipops_<track>.{png,pdf}          one per track
#   anchor_dotheatmap.{png,pdf}                 all tracks side by side
#   anchor_auc_coherence.{png,pdf}              diagnostic scatter (if cfg$scatter)
#   anchor_specificity_<track>.{png,pdf}        one per track
#   anchor_specificity_diagonal_<track>.{png,pdf}   one per track (skipped if no sig cell)
#   cluster_distinctive_categories.tsv          per-cluster most-distinctive category
#
#   Usage: Rscript R/plot_anchors.R --config configs/carey_rint15_plots.yaml
#          [--q-sig 0.05] [--rg-floor 0.10] [--min-clusters 5]

suppressPackageStartupMessages({ library(data.table); library(yaml) })
options(bitmapType = "cairo")

get_script_dir <- function() {
  a <- commandArgs(FALSE)
  f <- grep("^--file=", a, value = TRUE)
  if (length(f)) dirname(normalizePath(sub("^--file=", "", f[1]))) else getwd()
}
# This entry lives inside R/, so the engine modules are siblings in the same dir.
.SDIR <- get_script_dir()
for (m in c("io.R", "plot.R")) source(file.path(.SDIR, m))

parse_args <- function(a) {
  out <- list(config = NULL, q_sig = NULL, rg_floor = NULL, min_clusters = NULL)
  i <- 1L
  while (i <= length(a)) {
    if (a[i] == "--config")            { out$config <- a[i + 1L]; i <- i + 2L }
    else if (a[i] == "--q-sig")        { out$q_sig <- as.numeric(a[i + 1L]); i <- i + 2L }
    else if (a[i] == "--rg-floor")     { out$rg_floor <- as.numeric(a[i + 1L]); i <- i + 2L }
    else if (a[i] == "--min-clusters") { out$min_clusters <- as.integer(a[i + 1L]); i <- i + 2L }
    else i <- i + 1L
  }
  if (is.null(out$config))
    stop("usage: plot_anchors.R --config <plots.yaml> [--q-sig N] [--rg-floor N] [--min-clusters N]")
  out
}

run_plots <- function(config_path, q_sig = NULL, rg_floor = NULL, min_clusters = NULL) {
  cfg   <- yaml::read_yaml(config_path)
  sroot <- stage_root_of(config_path)
  # plot-config defaults (mirror plot_anchors.py / plot_specificity.py)
  defaults <- list(top_k = 8, lollipop_ncols = 3, scatter = TRUE, rg_cap = 0.55, q_sig = 0.05,
                   spec_rg_floor = 0.10, spec_min_clusters = 5)
  for (k in names(defaults)) if (is.null(cfg[[k]])) cfg[[k]] <- defaults[[k]]
  # CLI overrides
  if (!is.null(q_sig))        cfg[["q_sig"]] <- q_sig
  if (!is.null(rg_floor))     cfg[["spec_rg_floor"]] <- rg_floor
  if (!is.null(min_clusters)) cfg[["spec_min_clusters"]] <- min_clusters

  out_dir <- resolve_path(sroot, cfg[["out_dir"]])
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  written <- character(0)
  emit <- function(p) { message(sprintf("[write] %s", p)); written[[length(written) + 1L]] <<- p }

  tracks <- lapply(cfg[["tracks"]], load_track, stage_root = sroot)
  row_order <- natural_order(unlist(lapply(tracks, function(tr) tr[["s"]][["cluster_label"]])))

  # 1. lollipops (per track)
  for (tr in tracks) {
    nc <- min(cfg[["lollipop_ncols"]], length(unique(tr[["s"]][["cluster_label"]])))
    n_cl <- length(intersect(row_order, unique(tr[["s"]][["cluster_label"]])))
    nrw  <- ceiling(n_cl / max(nc, 1))
    p <- fig_lollipops(tr, row_order, cfg)
    png <- file.path(out_dir, sprintf("anchor_lollipops_%s.png", tr[["name"]]))
    save_fig(p, png, width = nc * 3.6 + 1.5, height = nrw * 2.3 + 0.6); emit(png)
  }

  # 2. dot-heatmap (combined)
  totcats <- sum(vapply(tracks, function(tr) length(unique(tr[["s"]][["category"]])), integer(1)))
  png <- file.path(out_dir, "anchor_dotheatmap.png")
  save_fig(fig_dotheatmap(tracks, row_order, cfg), png,
           width = totcats * 0.34 + 3.5, height = length(row_order) * 0.34 + 2.2); emit(png)

  # 3. AUC-vs-coherence scatter (combined)
  if (isTRUE(cfg[["scatter"]])) {
    png <- file.path(out_dir, "anchor_auc_coherence.png")
    save_fig(fig_scatter(tracks, cfg), png, width = length(tracks) * 4.6, height = 4.4); emit(png)
  }

  # 4-5. specificity heatmap + diagonal (per track) + distinctive table
  dist_all <- list()
  for (tr in tracks) {
    spec <- specificity(tr[["s"]], cfg[["q_sig"]], cfg[["spec_rg_floor"]], cfg[["spec_min_clusters"]])
    ro   <- natural_order(rownames(spec[["M"]]))
    ncat <- ncol(spec[["Z"]])
    png  <- file.path(out_dir, sprintf("anchor_specificity_%s.png", tr[["name"]]))
    save_fig(fig_specificity(spec, ro, tr[["name"]]), png,
             width = ncat * 0.42 + 3.5, height = length(ro) * 0.34 + 2.0); emit(png)
    pd <- fig_diagonal(spec, ro, tr[["name"]])
    if (is.null(pd)) {
      message(sprintf("[skip] anchor_specificity_diagonal_%s — no significant distinctive cells",
                      tr[["name"]]))
    } else {
      ncols <- length(diagonal_column_order(.boxed(spec, ro), ro))
      png <- file.path(out_dir, sprintf("anchor_specificity_diagonal_%s.png", tr[["name"]]))
      save_fig(pd, png, width = ncols * 0.6 + 3.5, height = length(ro) * 0.34 + 2.0); emit(png)
    }
    dist_all[[length(dist_all) + 1L]] <- distinctive_table(spec, tr[["name"]])
  }
  dist <- data.table::rbindlist(dist_all)
  dist_path <- file.path(out_dir, "cluster_distinctive_categories.tsv")
  data.table::fwrite(dist, dist_path, sep = "\t", na = "", quote = FALSE)
  emit(dist_path)

  message(sprintf("FINISHED ok | %d figure files + 1 table | out_dir: %s",
                  length(written) - 1L, out_dir))
  invisible(written)
}

if (sys.nframe() == 0L) {
  a <- parse_args(commandArgs(TRUE))
  run_plots(a$config, a$q_sig, a$rg_floor, a$min_clusters)
}
