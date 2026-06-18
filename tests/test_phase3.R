#!/usr/bin/env Rscript
# tests/test_phase3.R — Phase-3 sensitivity-sweep tests (plain stopifnot; no testthat).
# Uses the synthetic GenomicSEM .rds panel fixture (no external big-file dependency), so it is fast
# and fully deterministic. Run: Rscript tests/test_phase3.R
#   (build fixtures first if missing: Rscript tests/fixtures/make_synthetic_ldsc.R)
#
# Covers: (1) primary-slice parity — run_sensitivity's z==primary slice == an independent inline
# single-z recomputation (the RNG/loop-order regression guard); (2) thread-invariance — identical
# output for threads 1 vs 4; (3) label_stable correctness (constancy invariant + a forced flip);
# (4) gate monotonicity in z; (5) sanity (n_eff<=n, vif>=1).

sdir <- dirname(dirname(normalizePath(sub("^--file=", "",
          grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))))
for (m in c("io.R","gate.R","redundancy.R","score.R","label.R","ingest_rds.R","sensitivity.R"))
  source(file.path(sdir, "R", m))

ok <- function(name, cond) {
  if (!isTRUE(cond)) stop(sprintf("FAIL: %s", name))
  cat(sprintf("  ok  %s\n", name))
}
quiet <- function(...) invisible(NULL)
same  <- function(a, b) {                       # data.frame equality incl. NA / perm_p, order-sensitive
  rownames(a) <- NULL; rownames(b) <- NULL
  isTRUE(all.equal(a, b, check.attributes = FALSE))
}

cat("== Phase-3 sensitivity-sweep tests ==\n")

# ---- build the engine inputs from the synthetic .rds config (mirrors the driver's .rds route) -----
cfg_path <- file.path(sdir, "configs", "synthetic_rds.yaml")
cfg      <- load_config(cfg_path)
sroot    <- stage_root_of(cfg_path)
route    <- read_rds_route(resolve_path(sroot, cfg[["rds"]]), cfg, sroot, quiet)
df       <- route[["df"]]
override <- route[["trait_rg"]]                 # .rds-derived trait x trait block (100% coverage)
ont      <- read_ontology(resolve_path(sroot, cfg[["ontology"]]), cfg[["ontology_key"]])
z_prim   <- as.numeric(cfg[["h2_z_threshold"]])

# ---- (1) primary-slice parity: independent inline single-z run == sweep's primary slice ------------
# Inline copy of the legacy (pre-Phase-3) single-z scoring — an independent implementation of the same
# spec; catches any RNG-seed / loop-order drift in score_at_z.
legacy_single_z <- function() {
  set.seed(as.integer(cfg[["random_seed"]]))
  g   <- apply_universe_gate(df, cfg)           # gates at cfg$h2_z_threshold == primary z
  g   <- attach_ontology(g, ont, cfg[["ontology_key"]], cfg[["levels"]])
  sel <- select_corr_source(g, cfg, sroot, override, quiet)
  rows <- list()
  for (cl in unique(g[["cluster_label"]])) {
    gc <- g[g[["cluster_label"]] == cl, , drop = FALSE]; rownames(gc) <- NULL
    for (lvl in cfg[["levels"]]) rows <- c(rows, score_cluster_level(gc, lvl, sel[["corr"]], cfg))
  }
  rank_and_label(do.call(rbind, rows), cfg)[["ranked"]]
}
leg <- legacy_single_z()
sw1 <- run_sensitivity(df, ont, cfg, sroot, cfg[["z_vector"]], threads = 1L, override, quiet)
ok("primary-slice parity (sweep z==primary == inline single-z, incl perm_p)",
   same(leg, sw1[["primary"]][["ranked"]]))
ok("primary slice present in stacked scores",
   z_prim %in% sw1[["scores"]][["z_threshold"]])

# ---- (2) thread-invariance: threads 1 vs 4 produce identical stacked tables ------------------------
sw4 <- run_sensitivity(df, ont, cfg, sroot, cfg[["z_vector"]], threads = 4L, override, quiet)
ok("thread-invariance: scores identical (1 vs 4 workers)", same(sw1[["scores"]], sw4[["scores"]]))
ok("thread-invariance: labels identical (1 vs 4 workers)", same(sw1[["labels"]], sw4[["labels"]]))

# ---- (3) label_stable correctness ------------------------------------------------------------------
sl <- sw1[["labels"]]
ok("label_stable column present", "label_stable" %in% names(sl))
ok("every cluster appears at every swept z",
   all(table(sl[["cluster_label"]]) == length(sw1[["zs"]])))
# constancy invariant: label_stable iff auto_label constant across that cluster's z-rows
inv_ok <- all(vapply(unique(sl[["cluster_label"]]), function(cl) {
  r <- sl[sl[["cluster_label"]] == cl, ]
  unique(r[["label_stable"]]) %in% list(TRUE, FALSE) &&  # single value per cluster
    all(r[["label_stable"]] == (length(unique(r[["auto_label"]])) == 1L))
}, logical(1)))
ok("label_stable matches auto_label constancy for every cluster", inv_ok)
# forced flip: a very high z gates everything out -> ambiguous there -> any earlier real label flips
swf <- run_sensitivity(df, ont, cfg, sroot, c(z_prim, 1e3), threads = 1L, override, quiet)
slf <- swf[["labels"]]
flipped <- slf[slf[["z_threshold"]] == z_prim & slf[["auto_label"]] != "ambiguous", "cluster_label"]
ok("forced flip yields label_stable == FALSE",
   length(flipped) == 0 || all(!slf[["label_stable"]][slf[["cluster_label"]] %in% flipped]))
ok("clusters gated out at high z still appear (ambiguous/weak)",
   all(slf[["auto_label"]][slf[["z_threshold"]] == 1e3] == "ambiguous"))

# ---- (4) gate monotonicity: higher z -> not more rows ----------------------------------------------
g3 <- nrow(apply_universe_gate(df, cfg, 3)); g5 <- nrow(apply_universe_gate(df, cfg, 5))
ok("gate monotonic in z (n(z=5) <= n(z=3))", g5 <= g3)

# ---- (5) sanity: n_eff <= n and vif >= 1 across every swept z --------------------------------------
sc <- sw1[["scores"]]
ok("n_eff <= n at every z", all(sc[["n_eff"]] <= sc[["n"]]))
ok("vif >= 1 at every z",   all(sc[["vif"]] >= 1))

cat("\nALL PHASE-3 TESTS PASSED\n")
