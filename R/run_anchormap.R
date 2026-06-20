# run_anchormap.R - the engine orchestrator (gate -> redundancy -> score -> label -> z-sweep -> TSVs).

# ---- output column contracts -----------------------------------------------
.SCORE_COLS <- c("cluster_label","level","category","eligible","n","n_eff","n_hit","rho_bar","vif",
                 "auc_abs","auc_signed","perm_p","vif_z","vif_p","pooled_rg","pooled_rg_ci_lo",
                 "pooled_rg_ci_hi","coherence","mean_abs_rg","mean_signed_rg","odds_ratio",
                 "fisher_p","q","rank")
.LABEL_COLS <- c("cluster_label","auto_label","anchor_shape","anchor_margin","anchor_focus",
                 "n_sig_domains","top_auc","top_q","top_pooled_rg","top_coherence","profile")
# Sensitivity contracts: primary cols + z_threshold (+ label_stable for the labels table).
.SENS_SCORE_COLS <- c(.SCORE_COLS, "z_threshold")
.SENS_LABEL_COLS <- c(.LABEL_COLS, "z_threshold", "label_stable")

#' Run the AnchorMap engine
#'
#' Drives the full engine from a YAML config: ingest (rg long-table TSV **or** a GenomicSEM `.rds`),
#' reliability gate, within-category redundancy, competitive scoring, auto-labelling, and a parallel
#' reliability-threshold sensitivity sweep. Writes `category_anchor_scores.tsv`,
#' `cluster_anchor_labels.tsv`, the two `sensitivity_z_*` TSVs, and `anchormap.log` into the output
#' directory.
#'
#' @param config_path Path to a YAML config, or a bare shipped-config name (e.g. `"synthetic_rds"`).
#' @param threads Worker/thread count (`setDTthreads` + the z-sweep workers).
#' @param rds Optional GenomicSEM `.rds` input (overrides `cfg$rds`).
#' @param z_vector Optional numeric vector overriding `cfg$z_vector` for the sweep.
#' @param out_dir Optional output-directory override (else `cfg$out_dir`, else `results/<run_label>`).
#' @param run_label Optional run label (logging; fallback output dir when `out_dir` is unset).
#' @param rg_long,trait_rg,ontology Optional input-path overrides for the rg long-table, the
#'   trait x trait rg matrix, and the ontology TSV (else taken from the config).
#' @return Invisibly, a list with `ranked`, `labels`, `sens_scores`, `sens_labels`.
#' @export
run_anchormap <- function(config_path, threads = 1L, rds = NULL, z_vector = NULL,
                          out_dir = NULL, run_label = NULL,
                          rg_long = NULL, trait_rg = NULL, ontology = NULL) {
  t0 <- Sys.time()
  log <- character(0)
  emit <- function(...) { line <- sprintf(...); message(line); log[[length(log) + 1L]] <<- line }

  config_path <- resolve_config_path(config_path)
  emit("[start] %s  AnchorMap engine  config=%s", format(t0, "%Y-%m-%d %H:%M:%S"), config_path)

  data.table::setDTthreads(max(1L, threads))
  cfg   <- load_config(config_path)
  sroot <- stage_root_of(config_path)
  set.seed(as.integer(cfg[["random_seed"]]))

  # ---- CLI input overrides (resolved relative to the working dir) ----
  if (!is.null(rg_long))  cfg[["rg_long"]]         <- .abs_cwd(rg_long)
  if (!is.null(trait_rg)) cfg[["trait_rg_matrix"]] <- .abs_cwd(trait_rg)
  if (!is.null(ontology)) cfg[["ontology"]]        <- .abs_cwd(ontology)
  if (!is.null(run_label)) cfg[["run_label"]]      <- run_label

  emit("[config] random_seed=%d permutation_K=%d vif_correlation=%s",
       as.integer(cfg[["random_seed"]]), as.integer(cfg[["permutation_K"]]), cfg[["vif_correlation"]])

  # ---- output directory: --out-dir (cwd-relative) > cfg$out_dir > results/<run_label> ----
  rl <- run_label %||% cfg[["run_label"]]
  out_dir <-
    if (!is.null(out_dir))          .abs_cwd(out_dir)
    else if (!is.null(cfg[["out_dir"]])) resolve_path(sroot, cfg[["out_dir"]])
    else if (!is.null(rl))          resolve_path(sroot, file.path("results", rl))
    else stop("No output directory: pass out_dir, or set out_dir/run_label in the config.")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  ontology_path <- resolve_path(sroot, cfg[["ontology"]])

  # ---- input route: GenomicSEM .rds (Input C) OR rg long-TSV (Input A) ----
  rds_path <- if (!is.null(rds)) .abs_cwd(rds) else cfg[["rds"]]
  validate_config_sources(cfg, rds_active = !is.null(rds_path))
  trait_rg_override <- NULL
  if (!is.null(rds_path)) {
    rds_path <- resolve_path(sroot, rds_path)
    emit("[ingest] GenomicSEM .rds route: %s", rds_path)
    route <- read_rds_route(rds_path, cfg, sroot, emit)
    df <- route[["df"]]; trait_rg_override <- route[["trait_rg"]]
    emit("[ingest] %d cluster factors x %d panel traits -> %d long rows",
         route[["n_factors"]], route[["n_panel"]], nrow(df))
  } else {
    rg_long_path <- resolve_path(sroot, cfg[["rg_long"]])
    emit("[load] %s", rg_long_path)
    df <- read_long(rg_long_path)
  }

  g   <- apply_universe_gate(df, cfg)
  ont <- read_ontology(ontology_path, cfg[["ontology_key"]])
  g   <- attach_ontology(g, ont, cfg[["ontology_key"]], cfg[["levels"]])
  med <- as.integer(stats::median(table(g[["cluster_label"]])))
  emit("[gate] %s track: %d gated rows, %d clusters, %d traits/cluster (median); gate h2_z>%g",
       cfg[["trait_group"]], nrow(g), length(unique(g[["cluster_label"]])), med,
       as.numeric(cfg[["h2_z_threshold"]]))

  # ---- parallel z-threshold sensitivity sweep ----
  # Each z is a full independent re-run; the primary z (h2_z_threshold) is always folded in, so the
  # primary TSVs are its slice (primary == sweep[z==primary]).
  z_vec <- if (!is.null(z_vector)) as.numeric(z_vector) else as.numeric(cfg[["z_vector"]])
  zs    <- sort(unique(c(z_vec, as.numeric(cfg[["h2_z_threshold"]]))))

  # Build the trait x trait matrix ONCE over the loosest-z (superset) gated traits for the TSV route,
  # so the LDSC --rg summary isn't re-read per z (the .rds route already supplies the override).
  if (is.null(trait_rg_override) && cfg[["vif_correlation"]] %in% c("trait_rg","auto") &&
      !is.null(cfg[["trait_rg_matrix"]])) {
    loose_tids <- unique(apply_universe_gate(df, cfg, min(zs))[["trait_id"]])
    trait_rg_override <- build_trait_rg_matrix(resolve_path(sroot, cfg[["trait_rg_matrix"]]),
                                               loose_tids, isTRUE(cfg[["trait_rg_require_converged"]]))
  }

  workers <- max(1L, min(as.integer(threads), length(zs)))
  emit("[sweep] z in {%s} on %d worker(s)", paste(zs, collapse = ","), workers)
  sw <- run_sensitivity(df, ont, cfg, sroot, z_vec, threads, trait_rg_override, emit)
  for (m in sw[["meta"]])
    emit("[sweep z=%g] gated=%d clusters=%d source=%s coverage=%s", m[["z"]], m[["n_gated"]],
         m[["n_clusters"]], ifelse(is.na(m[["source"]]), "-", m[["source"]]),
         ifelse(is.na(m[["coverage"]]), "-", sprintf("%.0f%%", 100 * m[["coverage"]])))

  prim <- sw[["primary"]]
  if (is.null(prim[["ranked"]]) || !nrow(prim[["ranked"]]))
    stop("No category scores produced at the primary z - check gate thresholds / input.")

  # ---- primary TSVs (exact .SCORE_COLS/.LABEL_COLS contract) ----
  ranked <- prim[["ranked"]]; labels <- prim[["labels"]]
  ranked <- ranked[order(ranked[["level"]], ranked[["cluster_label"]], ranked[["rank"]]), .SCORE_COLS]
  ranked[["eligible"]] <- ifelse(ranked[["eligible"]], "True", "False")   # match pandas bool repr
  labels <- labels[order(labels[["cluster_label"]]), .LABEL_COLS]
  scores_path <- file.path(out_dir, "category_anchor_scores.tsv")
  labels_path <- file.path(out_dir, "cluster_anchor_labels.tsv")
  data.table::fwrite(ranked, scores_path, sep = "\t", na = "", quote = FALSE)
  data.table::fwrite(labels, labels_path, sep = "\t", na = "", quote = FALSE)
  emit("[write] %s", scores_path)
  emit("[write] %s", labels_path)

  # ---- sensitivity TSVs (.SENS_* contracts: primary cols + z_threshold [+ label_stable]) ----
  ss <- sw[["scores"]]
  ss <- ss[order(ss[["z_threshold"]], ss[["level"]], ss[["cluster_label"]], ss[["rank"]]), .SENS_SCORE_COLS]
  ss[["eligible"]] <- ifelse(ss[["eligible"]], "True", "False")
  sl <- sw[["labels"]]
  sl <- sl[order(sl[["z_threshold"]], sl[["cluster_label"]]), .SENS_LABEL_COLS]
  sl[["label_stable"]] <- ifelse(sl[["label_stable"]], "True", "False")
  sens_scores_path <- file.path(out_dir, "sensitivity_z_scores.tsv")
  sens_labels_path <- file.path(out_dir, "sensitivity_z_labels.tsv")
  data.table::fwrite(ss, sens_scores_path, sep = "\t", na = "", quote = FALSE)
  data.table::fwrite(sl, sens_labels_path, sep = "\t", na = "", quote = FALSE)
  emit("[write] %s", sens_scores_path)
  emit("[write] %s", sens_labels_path)

  n_stable <- sum(tapply(sl[["label_stable"]], sl[["cluster_label"]], function(v) v[1] == "True"))
  emit("[stable] %d/%d clusters label-stable across z in {%s}",
       n_stable, length(sw[["all_clusters"]]), paste(sw[["zs"]], collapse = ","))
  for (i in seq_len(nrow(labels)))
    emit("  %-10s -> %-18s [%s]", labels[["cluster_label"]][i],
         labels[["auto_label"]][i], labels[["anchor_shape"]][i])

  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  emit("FINISHED ok | %d score rows, %d clusters | %d z swept | %.1fs | outputs: %s, %s, %s, %s",
       nrow(ranked), nrow(labels), length(sw[["zs"]]), elapsed, basename(scores_path),
       basename(labels_path), basename(sens_scores_path), basename(sens_labels_path))
  writeLines(log, file.path(out_dir, "anchormap.log"))
  invisible(list(ranked = ranked, labels = labels, sens_scores = ss, sens_labels = sl))
}
