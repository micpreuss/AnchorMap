# sensitivity.R — Phase-3 parallel z-threshold sensitivity sweep.
#
# Re-runs the whole single-z engine (gate -> redundancy -> score -> label) at each h2-reliability
# threshold in the z-vector, in parallel over z via future.apply, then stacks the two output tables
# and flags per-cluster auto-label stability.
#
# Determinism is engineered, not hoped-for: each z-task re-seeds with cfg$random_seed and preserves
# the exact cluster x level loop order of anchor_map.R. Hence (a) the z == h2_z_threshold task
# reproduces the Phase-1/2 single-z output bit-for-bit (incl. perm_p), and (b) the whole sweep is
# invariant to worker count AND backend. future.seed=TRUE only silences the parallel-RNG warning;
# the inner set.seed dominates. score.R / label.R are deliberately untouched — perm_p stays serial
# so its RNG stream matches the legacy run; the sweep parallelizes the OUTER z axis only.

suppressPackageStartupMessages({ library(future.apply) })

# Map FUN over X with `threads` workers via future.apply. workers==1 -> sequential plan (plan
# selection, NOT an availability fallback; future.apply is a hard dependency). setDTthreads(1)
# avoids data.table x workers oversubscription.
#
# future.globals = FALSE: we use only `sequential` and `multicore` (fork) plans, so the worker sees
# the parent's globals directly (same process / fork-inherited memory) — no export needed. This also
# sidesteps future's globals scanner, which trips over the engine functions sourced into globalenv.
# (It would be unsafe only under multisession/cluster plans, which parallel_lapply never selects.)
parallel_lapply <- function(X, FUN, threads = 1L) {
  workers <- max(1L, min(as.integer(threads), length(X)))
  data.table::setDTthreads(1L)
  if (workers == 1L) future::plan("sequential")           # sequential takes no `workers` arg
  else               future::plan("multicore", workers = workers)
  on.exit(future::plan("sequential"), add = TRUE)
  future.apply::future_lapply(X, FUN, future.seed = TRUE, future.globals = FALSE)
}

# One all-ambiguous label row, byte-matching rank_and_label's empty-cluster branch (label.R:61-67).
.ambiguous_label_row <- function(cl) {
  data.frame(cluster_label = cl, auto_label = "ambiguous", anchor_shape = "weak",
             anchor_margin = NA_real_, anchor_focus = NA_real_, n_sig_domains = 0L,
             top_auc = NA_real_, top_q = NA_real_, top_pooled_rg = NA_real_,
             top_coherence = NA_real_, profile = "", stringsAsFactors = FALSE)
}

# Score the whole engine at a single reliability threshold z. Deterministic: set.seed(random_seed)
# then the exact cluster x level loop of anchor_map.R, so z == h2_z_threshold reproduces the legacy
# single-z output bit-for-bit. Returns ranked/labels (NULL when nothing scored) + per-z metadata.
score_at_z <- function(df, ont, cfg, sroot, z, trait_rg_override = NULL, emit = message) {
  # Pin the RNG to R's defaults BEFORE seeding: future.seed=TRUE switches the kind to L'Ecuyer-CMRG,
  # and set.seed(seed) with no `kind` keeps the current kind — so without this the sweep would draw
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

# Parallel z-sweep. Re-runs score_at_z over sort(unique(c(z_vector, h2_z_threshold))) (the primary z
# is always folded in, so the parity slice always exists), stacks the two tables with a z_threshold
# column, pads every z's labels to the union cluster set (so a cluster gated out at some z still shows
# ambiguous/weak there), and flags per-cluster label_stable (auto_label constant across all swept z).
# Returns the stacked frames (eligible/label_stable left logical for the driver's write contract), the
# z == primary slice (for the primary TSVs, from the SAME computation), and per-z metadata for logging.
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
      pad <- do.call(rbind, lapply(miss, .ambiguous_label_row))
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
