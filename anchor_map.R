#!/usr/bin/env Rscript
# anchor_map.R — AnchorMap engine CLI (Phase 2).
#
# Faithful R port of cluster_anchoring/anchor_categories.py. Reads the identical YAML configs.
#   Usage: Rscript anchor_map.R --config <config.yaml> [--threads N] [--rds <ldsc_output.rds>]
#
# Writes (flat, mirroring the Python reference, into the config's out_dir):
#   category_anchor_scores.tsv   one row per (cluster, level, category)
#   cluster_anchor_labels.tsv    one row per cluster (auto-label + anchor shape)
#   anchormap.log                timestamped steps, ending in a FINISHED statement
# Still single-z (the parallel z-sweep is Phase 3). Phase 2 adds: the GenomicSEM .rds input route
# (--rds / cfg$rds, via R/ingest_rds.R) and the `vif_correlation: auto` redundancy auto-fallback.
# Explicit `trait_rg`/`cluster_profile` modes keep Phase-1 behaviour byte-for-byte.

suppressPackageStartupMessages({ library(data.table); library(yaml) })

# ---- locate + source the engine modules ------------------------------------
get_script_dir <- function() {
  a <- commandArgs(FALSE)
  f <- grep("^--file=", a, value = TRUE)
  if (length(f)) dirname(normalizePath(sub("^--file=", "", f[1]))) else getwd()
}
.SDIR <- get_script_dir()
for (m in c("io.R","gate.R","redundancy.R","score.R","label.R","ingest_rds.R"))
  source(file.path(.SDIR, "R", m))

# ---- tiny arg parser (no argparse dependency) ------------------------------
parse_args <- function(a) {
  out <- list(config = NULL, threads = 1L, rds = NULL)
  i <- 1L
  while (i <= length(a)) {
    if (a[i] == "--config")  { out$config  <- a[i + 1L]; i <- i + 2L }
    else if (a[i] == "--threads") { out$threads <- as.integer(a[i + 1L]); i <- i + 2L }
    else if (a[i] == "--rds")     { out$rds     <- a[i + 1L]; i <- i + 2L }
    else i <- i + 1L
  }
  if (is.null(out$config)) stop("usage: anchor_map.R --config <config.yaml> [--threads N] [--rds <file>]")
  out
}

# ---- output column contracts (must match the Python reference) -------------
.SCORE_COLS <- c("cluster_label","level","category","eligible","n","n_eff","n_hit","rho_bar","vif",
                 "auc_abs","auc_signed","perm_p","vif_z","vif_p","pooled_rg","pooled_rg_ci_lo",
                 "pooled_rg_ci_hi","coherence","mean_abs_rg","mean_signed_rg","odds_ratio",
                 "fisher_p","q","rank")
.LABEL_COLS <- c("cluster_label","auto_label","anchor_shape","anchor_margin","anchor_focus",
                 "n_sig_domains","top_auc","top_q","top_pooled_rg","top_coherence","profile")

run_anchormap <- function(config_path, threads = 1L, rds = NULL) {
  t0 <- Sys.time()
  log <- character(0)
  emit <- function(...) { line <- sprintf(...); message(line); log[[length(log) + 1L]] <<- line }
  emit("[start] %s  AnchorMap engine (Phase 2)  config=%s",
       format(t0, "%Y-%m-%d %H:%M:%S"), config_path)

  data.table::setDTthreads(max(1L, threads))
  cfg   <- load_config(config_path)
  sroot <- stage_root_of(config_path)
  set.seed(as.integer(cfg[["random_seed"]]))
  emit("[config] random_seed=%d permutation_K=%d vif_correlation=%s",
       as.integer(cfg[["random_seed"]]), as.integer(cfg[["permutation_K"]]), cfg[["vif_correlation"]])

  ontology <- resolve_path(sroot, cfg[["ontology"]])
  out_dir  <- resolve_path(sroot, cfg[["out_dir"]])
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  # ---- input route: GenomicSEM .rds (Input C) OR rg long-TSV (Input A) ----
  rds_path <- if (!is.null(rds)) rds else cfg[["rds"]]
  trait_rg_override <- NULL
  if (!is.null(rds_path)) {
    rds_path <- resolve_path(sroot, rds_path)
    emit("[ingest] GenomicSEM .rds route: %s", rds_path)
    route <- read_rds_route(rds_path, cfg, sroot, emit)
    df <- route[["df"]]; trait_rg_override <- route[["trait_rg"]]
    emit("[ingest] %d cluster factors x %d panel traits -> %d long rows",
         route[["n_factors"]], route[["n_panel"]], nrow(df))
  } else {
    rg_long <- resolve_path(sroot, cfg[["rg_long"]])
    emit("[load] %s", rg_long)
    df <- read_long(rg_long)
  }

  g   <- apply_universe_gate(df, cfg)
  ont <- read_ontology(ontology, cfg[["ontology_key"]])
  g   <- attach_ontology(g, ont, cfg[["ontology_key"]], cfg[["levels"]])
  med <- as.integer(stats::median(table(g[["cluster_label"]])))
  emit("[gate] %s track: %d gated rows, %d clusters, %d traits/cluster (median); gate h2_z>%g",
       cfg[["trait_group"]], nrow(g), length(unique(g[["cluster_label"]])), med,
       as.numeric(cfg[["h2_z_threshold"]]))

  # ---- redundancy source (Phase-1 explicit modes unchanged; `auto` self-selects) ----
  sel  <- select_corr_source(g, cfg, sroot, trait_rg_override, emit)
  corr <- sel[["corr"]]
  emit("[vif] reason: %s", sel[["reason"]])

  rows <- list()
  for (cl in unique(g[["cluster_label"]])) {
    gc <- g[g[["cluster_label"]] == cl, , drop = FALSE]; rownames(gc) <- NULL
    for (lvl in cfg[["levels"]]) rows <- c(rows, score_cluster_level(gc, lvl, corr, cfg))
  }
  if (!length(rows)) stop("No category scores produced - check gate thresholds / input.")
  scores <- do.call(rbind, rows)
  rl <- rank_and_label(scores, cfg)
  ranked <- rl[["ranked"]]; labels <- rl[["labels"]]

  # order + write the two TSVs in the exact contract schema
  ranked <- ranked[order(ranked[["level"]], ranked[["cluster_label"]], ranked[["rank"]]), .SCORE_COLS]
  ranked[["eligible"]] <- ifelse(ranked[["eligible"]], "True", "False")   # match pandas bool repr
  labels <- labels[order(labels[["cluster_label"]]), .LABEL_COLS]

  scores_path <- file.path(out_dir, "category_anchor_scores.tsv")
  labels_path <- file.path(out_dir, "cluster_anchor_labels.tsv")
  data.table::fwrite(ranked, scores_path, sep = "\t", na = "", quote = FALSE)
  data.table::fwrite(labels, labels_path, sep = "\t", na = "", quote = FALSE)
  emit("[write] %s", scores_path)
  emit("[write] %s", labels_path)
  for (i in seq_len(nrow(labels)))
    emit("  %-10s -> %-18s [%s]", labels[["cluster_label"]][i],
         labels[["auto_label"]][i], labels[["anchor_shape"]][i])

  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  emit("FINISHED ok | %d score rows, %d clusters | %.1fs | outputs: %s, %s",
       nrow(ranked), nrow(labels), elapsed, basename(scores_path), basename(labels_path))
  writeLines(log, file.path(out_dir, "anchormap.log"))
  invisible(list(ranked = ranked, labels = labels))
}

if (sys.nframe() == 0L) {
  a <- parse_args(commandArgs(TRUE))
  run_anchormap(a$config, a$threads, a$rds)
}
