# Sensitivity-sweep tests (uses the synthetic .rds panel fixture via the shipped config).

quiet <- function(...) invisible(NULL)
same  <- function(a, b) {            # data.frame equality incl. NA / perm_p, order-sensitive
  rownames(a) <- NULL; rownames(b) <- NULL
  isTRUE(all.equal(a, b, check.attributes = FALSE))
}

cfg_path <- system.file("configs/synthetic_rds.yaml", package = "anchormap")
cfg      <- load_config(cfg_path)
sroot    <- stage_root_of(cfg_path)
route    <- read_rds_route(resolve_path(sroot, cfg[["rds"]]), cfg, sroot, quiet)
df       <- route[["df"]]
override <- route[["trait_rg"]]
ont      <- read_ontology(resolve_path(sroot, cfg[["ontology"]]), cfg[["ontology_key"]])
z_prim   <- as.numeric(cfg[["h2_z_threshold"]])

# Inline single-z scoring — an independent implementation of the same spec; catches RNG-seed /
# loop-order drift in score_at_z.
legacy_single_z <- function() {
  set.seed(as.integer(cfg[["random_seed"]]))
  g   <- apply_universe_gate(df, cfg)
  g   <- attach_ontology(g, ont, cfg[["ontology_key"]], cfg[["levels"]])
  sel <- select_corr_source(g, cfg, sroot, override, quiet)
  rows <- list()
  for (cl in unique(g[["cluster_label"]])) {
    gc <- g[g[["cluster_label"]] == cl, , drop = FALSE]; rownames(gc) <- NULL
    for (lvl in cfg[["levels"]]) rows <- c(rows, score_cluster_level(gc, lvl, sel[["corr"]], cfg))
  }
  rank_and_label(do.call(rbind, rows), cfg)[["ranked"]]
}

test_that("primary-slice parity (sweep z==primary == inline single-z, incl perm_p)", {
  leg <- legacy_single_z()
  sw1 <- run_sensitivity(df, ont, cfg, sroot, cfg[["z_vector"]], threads = 1L, override, quiet)
  expect_true(same(leg, sw1[["primary"]][["ranked"]]))
  expect_true(z_prim %in% sw1[["scores"]][["z_threshold"]])
})

test_that("thread-invariance: threads 1 vs 4 produce identical stacked tables", {
  sw1 <- run_sensitivity(df, ont, cfg, sroot, cfg[["z_vector"]], threads = 1L, override, quiet)
  sw4 <- run_sensitivity(df, ont, cfg, sroot, cfg[["z_vector"]], threads = 4L, override, quiet)
  expect_true(same(sw1[["scores"]], sw4[["scores"]]))
  expect_true(same(sw1[["labels"]], sw4[["labels"]]))
})

test_that("label_stable correctness (constancy invariant + forced flip)", {
  sw1 <- run_sensitivity(df, ont, cfg, sroot, cfg[["z_vector"]], threads = 1L, override, quiet)
  sl <- sw1[["labels"]]
  expect_true("label_stable" %in% names(sl))
  expect_true(all(table(sl[["cluster_label"]]) == length(sw1[["zs"]])))
  inv_ok <- all(vapply(unique(sl[["cluster_label"]]), function(cl) {
    r <- sl[sl[["cluster_label"]] == cl, ]
    length(unique(r[["label_stable"]])) == 1L &&
      all(r[["label_stable"]] == (length(unique(r[["auto_label"]])) == 1L))
  }, logical(1)))
  expect_true(inv_ok)

  swf <- run_sensitivity(df, ont, cfg, sroot, c(z_prim, 1e3), threads = 1L, override, quiet)
  slf <- swf[["labels"]]
  flipped <- slf[slf[["z_threshold"]] == z_prim & slf[["auto_label"]] != "ambiguous", "cluster_label"]
  expect_true(length(flipped) == 0 || all(!slf[["label_stable"]][slf[["cluster_label"]] %in% flipped]))
  expect_true(all(slf[["auto_label"]][slf[["z_threshold"]] == 1e3] == "ambiguous"))
})

test_that("gate monotonicity: higher z -> not more rows", {
  g3 <- nrow(apply_universe_gate(df, cfg, 3)); g5 <- nrow(apply_universe_gate(df, cfg, 5))
  expect_lte(g5, g3)
})

test_that("sanity: n_eff <= n and vif >= 1 across every swept z", {
  sw1 <- run_sensitivity(df, ont, cfg, sroot, cfg[["z_vector"]], threads = 1L, override, quiet)
  sc <- sw1[["scores"]]
  expect_true(all(sc[["n_eff"]] <= sc[["n"]]))
  expect_true(all(sc[["vif"]] >= 1))
})
