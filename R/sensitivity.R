# sensitivity.R - Phase-3 parallel z-threshold sensitivity sweep.
#
# Re-runs the whole single-z engine (gate -> redundancy -> score -> label) at each h2-reliability
# threshold in the z-vector, in parallel over z via future.apply, then stacks the two output tables
# and flags per-cluster auto-label stability.
#
# Determinism is engineered, not hoped-for: each z-task re-seeds with cfg$random_seed and preserves
# the exact cluster x level loop order of anchor_map.R. Hence (a) the z == h2_z_threshold task
# reproduces the Phase-1/2 single-z output bit-for-bit (incl. perm_p), and (b) the whole sweep is
# invariant to worker count AND backend. future.seed=TRUE only silences the parallel-RNG warning;
# the inner set.seed dominates. score.R / label.R are deliberately untouched - perm_p stays serial
# so its RNG stream matches the single-z run; the sweep parallelizes the OUTER z axis only.

# Map FUN over X with `threads` workers via future.apply. workers==1 -> sequential plan (plan
# selection, NOT an availability fallback; future.apply is a hard dependency). setDTthreads(1)
# avoids data.table x workers oversubscription.
#
# future.globals = FALSE: we use only `sequential` and `multicore` (fork) plans, so the worker sees
# the calling process's globals directly (fork-inherited memory) - no export needed, and it sidesteps
# future's globals scanner. (It would be unsafe only under multisession/cluster plans, which
# parallel_lapply never selects.)
parallel_lapply <- function(X, FUN, threads = 1L) {
  workers <- max(1L, min(as.integer(threads), length(X)))
  data.table::setDTthreads(1L)
  if (workers == 1L) future::plan("sequential")           # sequential takes no `workers` arg
  else               future::plan("multicore", workers = workers)
  on.exit(future::plan("sequential"), add = TRUE)
  future.apply::future_lapply(X, FUN, future.seed = TRUE, future.globals = FALSE)
}

#' Score the engine at a single reliability threshold
#'
#' Runs the full single-z engine (gate -> redundancy -> score -> label) at one h2-reliability
#' threshold `z`. Deterministic: re-seeds with `cfg$random_seed` and preserves the exact
#' cluster x level loop order, so `z == cfg$h2_z_threshold` reproduces the single-z primary output
#' bit-for-bit.
#'
#' @param df Gated long-table data frame (Input A schema).
#' @param ont Ontology data frame (read by `read_ontology()`).
#' @param cfg Config list (from `load_config()`).
#' @param sroot Stage root for resolving relative paths.
#' @param z Numeric h2-reliability threshold.
#' @param trait_rg_override Optional precomputed trait x trait rg matrix.
#' @param emit Logging function (default [message()]).
#' @return A list with `ranked`, `labels`, and per-z metadata.
#' @export
score_at_z <- function(df, ont, cfg, sroot, z, trait_rg_override = NULL, emit = message) {
  # Pin the RNG to R's defaults BEFORE seeding: future.seed=TRUE switches the kind to L'Ecuyer-CMRG,
  # and set.seed(seed) with no `kind` keeps the current kind - so without this the sweep would draw
  # from a different generator than the legacy serial run and perm_p would not match.
  set.seed(as.integer(cfg[["random_seed"]]),
           kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
  g <- apply_universe_gate(df, cfg, z)
  g <- attach_ontology(g, ont, cfg[["ontology_key"]], cfg[["levels"]])
  out <- function(ranked, labels, source = NA_character_, coverage = NA_real_)
    list(ranked = ranked, labels = labels, z = z, n_gated = nrow(g),
         n_clusters = length(unique(g[["cluster_label"]])), source = source, coverage = coverage)
  if (!nrow(g)) return(out(NULL, NULL))

  sel <- select_corr_source(g, cfg, sroot, trait_rg_override, emit)
  rows <- list()
  for (cl in unique(g[["cluster_label"]])) {              # SAME order as anchor_map.R (parity)
    gc <- g[g[["cluster_label"]] == cl, , drop = FALSE]; rownames(gc) <- NULL
    for (lvl in cfg[["levels"]]) rows <- c(rows, score_cluster_level(gc, lvl, sel[["corr"]], cfg))
  }
  if (!length(rows)) return(out(NULL, NULL, sel[["source"]], sel[["coverage"]]))
  rl <- rank_and_label(do.call(rbind, rows), cfg)
  out(rl[["ranked"]], rl[["labels"]], sel[["source"]], sel[["coverage"]])
}

#' Parallel reliability-threshold sensitivity sweep
#'
#' Re-runs [score_at_z()] over `sort(unique(c(z_vector, cfg$h2_z_threshold)))` in parallel over z
#' (the primary z is always folded in, so the primary slice always exists), stacks the two output
#' tables with a `z_threshold` column, pads every z's labels to the union cluster set, and flags
#' per-cluster `label_stable` (auto_label constant across all swept z).
#'
#' @inheritParams score_at_z
#' @param z_vector Numeric vector of h2-reliability thresholds to sweep.
#' @param threads Worker count for the z-axis (uses `future` multicore/sequential).
#' @return A list with stacked `scores`/`labels`, the `primary` (z == primary) result, `zs`, and
#'   per-z `meta`.
#' @export
run_sensitivity <- function(df, ont, cfg, sroot, z_vector, threads = 1L,
                            trait_rg_override = NULL, emit = message) {
  z_primary <- as.numeric(cfg[["h2_z_threshold"]])
  zs <- sort(unique(c(as.numeric(z_vector), z_primary)))
  emit_quiet <- function(...) message(sprintf(...))     # per-z [vif]/WARN -> stderr, not the main log

  res <- parallel_lapply(zs, function(z)
    score_at_z(df, ont, cfg, sroot, z, trait_rg_override, emit_quiet), threads)
  names(res) <- as.character(zs)

  # stacked scores (eligible left logical; driver applies the .SENS_SCORE_COLS contract)
  scores <- do.call(rbind, lapply(zs, function(z) {
    r <- res[[as.character(z)]][["ranked"]]
    if (is.null(r) || !nrow(r)) return(NULL)
    r[["z_threshold"]] <- z; r
  }))

  # union cluster universe across the sweep, then pad each z's labels to it
  all_clusters <- sort(unique(unlist(lapply(zs, function(z)
    res[[as.character(z)]][["labels"]][["cluster_label"]]))))
  labels_by_z <- lapply(zs, function(z) {
    lab <- res[[as.character(z)]][["labels"]]
    present <- if (is.null(lab)) character(0) else lab[["cluster_label"]]
    miss <- setdiff(all_clusters, present)
    if (length(miss)) {
      pad <- do.call(rbind, lapply(miss, ambiguous_label_row, cfg = cfg))
      lab <- if (is.null(lab) || !nrow(lab)) pad else rbind(lab, pad)
    }
    if (is.null(lab)) return(NULL)
    lab <- lab[order(lab[["cluster_label"]]), , drop = FALSE]
    lab[["z_threshold"]] <- z
    lab
  })
  labels <- do.call(rbind, labels_by_z); rownames(labels) <- NULL

  # per-cluster label stability across all swept z (broadcast onto each z-row)
  if (!is.null(labels) && nrow(labels)) {
    stab <- tapply(labels[["auto_label"]], labels[["cluster_label"]],
                   function(v) length(unique(v)) == 1L)
    labels[["label_stable"]] <- as.logical(stab[labels[["cluster_label"]]])
  }

  meta <- lapply(zs, function(z)
    res[[as.character(z)]][c("z", "n_gated", "n_clusters", "source", "coverage")])
  list(scores = scores, labels = labels, primary = res[[as.character(z_primary)]],
       zs = zs, z_primary = z_primary, all_clusters = all_clusters, meta = meta)
}
