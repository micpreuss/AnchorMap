# GenomicSEM .rds ingestion + redundancy auto-fallback (uses the shipped synthetic fixtures).

test_that("vech indexing (column-major lower triangle)", {
  expect_identical(vech_index(3)[["2,1"]], 2L)
  expect_identical(vech_index(3)[["3,3"]], 6L)
  expect_identical(length(vech_index(5)), 15L)
  expect_identical(vech_index(5)[["5,4"]], 14L)
})

test_that("reader + analytic standardization", {
  x <- read_ldsc_rds(fx("synthetic_ldsc.rds"))
  S <- x[["S"]]; k <- nrow(S)
  S_std <- standardize_S(S)
  expect_true(approx(S_std, stats::cov2cor(S), 1e-12))
  expect_true(approx(diag(S_std), rep(1, k), 1e-12))
  expect_true(approx(S_std, t(S_std), 1e-12))

  tmpnl <- tempfile(fileext = ".rds"); saveRDS(1:3, tmpnl)
  expect_error(read_ldsc_rds(tmpnl))
  bad <- list(S = matrix(1, 3, 3, dimnames = list(letters[1:3], letters[1:3])), V = matrix(1, 5, 5))
  tmpbad <- tempfile(fileext = ".rds"); saveRDS(bad, tmpbad)
  expect_error(read_ldsc_rds(tmpbad))
})

test_that("delta-method rg_se == numeric difference (the key gate)", {
  x <- read_ldsc_rds(fx("synthetic_ldsc.rds"))
  S <- x[["S"]]; V <- x[["V"]]; k <- nrow(S)
  rg_se <- rg_se_matrix(S, V)
  vm  <- vech_index(k)
  pos <- function(a, b) vm[[sprintf("%d,%d", max(a, b), min(a, b))]]
  r_fun <- function(Sii, Sjj, Sij) Sij / sqrt(Sii * Sjj)
  check_pair <- function(i, j) {
    Sii <- S[i, i]; Sjj <- S[j, j]; Sij <- S[i, j]; eps <- 1e-6
    g_num <- c(
      (r_fun(Sii, Sjj, Sij + eps) - r_fun(Sii, Sjj, Sij - eps)) / (2 * eps),
      (r_fun(Sii + eps, Sjj, Sij) - r_fun(Sii - eps, Sjj, Sij)) / (2 * eps),
      (r_fun(Sii, Sjj + eps, Sij) - r_fun(Sii, Sjj - eps, Sij)) / (2 * eps))
    idx3 <- c(pos(i, j), pos(i, i), pos(j, j))
    se_num <- sqrt(as.numeric(t(g_num) %*% V[idx3, idx3] %*% g_num))
    abs(rg_se[i, j] - se_num)
  }
  for (pr in list(c(2,1), c(3,1), c(4,3), c(5,2), c(5,4)))
    expect_lt(check_pair(pr[1], pr[2]), 1e-6)
  expect_true(approx(rg_se, t(rg_se), 1e-12))
  expect_true(approx(diag(rg_se), rep(0, k), 1e-12))
})

test_that("partition (regex default + explicit override + empty error)", {
  S <- read_ldsc_rds(fx("synthetic_ldsc.rds"))[["S"]]
  p <- partition_S(rownames(S), list(cluster_factor_pattern = "^C[0-9]", cluster_factors = NULL))
  expect_setequal(p$factors, c("C0", "C5_sub0"))
  expect_setequal(p$panel, c("BMI", "WT", "HT"))
  p2 <- partition_S(rownames(S), list(cluster_factors = c("BMI"), cluster_factor_pattern = "^C[0-9]"))
  expect_setequal(p2$factors, "BMI")
  expect_setequal(p2$panel, c("C0","C5_sub0","WT","HT"))
  expect_error(partition_S(rownames(S), list(cluster_factor_pattern = "^ZZZ", cluster_factors = NULL)))
  expect_error(partition_S(rownames(S), list(cluster_factor_pattern = ".", cluster_factors = NULL)))
})

test_that("rds_to_long schema + negative-h2 guard", {
  x <- read_ldsc_rds(fx("synthetic_ldsc.rds"))
  S <- x[["S"]]; V <- x[["V"]]
  S_std <- standardize_S(S); rg_se <- rg_se_matrix(S, V)
  p <- partition_S(rownames(S), list(cluster_factor_pattern = "^C[0-9]", cluster_factors = NULL))
  h2_se <- h2_se_vector(V, rownames(S))
  cfg5  <- modifyList(default_config(), list(ontology_key = "trait_id", trait_group = "disease"))
  df5   <- rds_to_long(S, S_std, rg_se, h2_se, p$factors, p$panel, cfg5, trait_meta = NULL)
  expect_true(all(.LONG_REQUIRED %in% names(df5)))
  expect_identical(names(df5), .LONG_REQUIRED)
  expect_true(is.numeric(df5$rg) && is.numeric(df5$rg_se))
  expect_identical(nrow(df5), 6L)
  expect_true(all(df5$status == "success") && all(df5$ldsc_converged))
  expect_error(rds_to_long(S, S_std, rg_se, h2_se, p$factors, p$panel,
                           modifyList(cfg5, list(ontology_key = "trait_category")), NULL))

  xn <- read_ldsc_rds(fx("synthetic_ldsc_negh2.rds"))
  Sn <- xn[["S"]]; Vn <- xn[["V"]]
  rg_se_n <- rg_se_matrix(Sn, Vn)
  expect_true(is.na(rg_se_n["HT", "BMI"]) && is.na(rg_se_n["BMI", "HT"]))
  dfn <- rds_to_long(Sn, standardize_S(Sn), rg_se_n, h2_se_vector(Vn, rownames(Sn)),
                     p$factors, p$panel, cfg5, NULL)
  htr <- dfn[dfn$trait_id == "HT", ]
  expect_true(all(htr$negative_h2) && all(htr$status == "failed") && all(!htr$ldsc_converged))
  gn <- apply_universe_gate(dfn, cfg5)
  expect_false("HT" %in% gn$trait_id)
  expect_true(all(c("BMI","WT") %in% gn$trait_id))
})

test_that("round-trip: .rds route == equivalent TSV route (deterministic scores)", {
  cfgP <- modifyList(default_config(), list(
    ontology_key = "trait_id", levels = c("anthro_class"), primary_level = "anthro_class",
    trait_group = "disease", vif_correlation = "trait_rg", permutation_K = 2000, random_seed = 1))
  route <- read_rds_route(fx("synthetic_ldsc_panel.rds"), cfgP, ".", emit = function(...) invisible(NULL))
  df_rds   <- route$df
  trait_rg <- route$trait_rg
  ont <- read_ontology(fx("synthetic_panel_ontology.tsv"), "trait_id")

  tmp_long <- tempfile(fileext = ".tsv"); data.table::fwrite(df_rds, tmp_long, sep = "\t", na = "", quote = FALSE)
  df_tsv   <- read_long(tmp_long)
  panel    <- rownames(trait_rg)
  edges <- do.call(rbind, lapply(seq_len(length(panel) - 1L), function(a)
             do.call(rbind, lapply((a + 1L):length(panel), function(b)
               data.frame(p1 = panel[a], p2 = panel[b], rg = trait_rg[a, b], CONVERGED = "TRUE",
                          stringsAsFactors = FALSE)))))
  tmp_sum <- tempfile(fileext = ".tsv"); data.table::fwrite(edges, tmp_sum, sep = "\t", quote = FALSE)

  score_route <- function(df_in, corr_in) {
    g <- apply_universe_gate(df_in, cfgP)
    g <- attach_ontology(g, ont, "trait_id", c("anthro_class"))
    set.seed(as.integer(cfgP$random_seed))
    rows <- list()
    for (cl in unique(g$cluster_label)) {
      gc <- g[g$cluster_label == cl, , drop = FALSE]; rownames(gc) <- NULL
      rows <- c(rows, score_cluster_level(gc, "anthro_class", corr_in, cfgP))
    }
    rl <- rank_and_label(do.call(rbind, rows), cfgP)
    rl$ranked[order(rl$ranked$cluster_label, rl$ranked$category), ]
  }
  sc_rds <- score_route(df_rds, trait_rg)
  sc_tsv <- score_route(df_tsv, build_trait_rg_matrix(tmp_sum, unique(df_tsv$trait_id), TRUE))
  expect_true(nrow(sc_rds) == nrow(sc_tsv) && nrow(sc_rds) > 0)
  expect_identical(paste(sc_rds$cluster_label, sc_rds$category),
                   paste(sc_tsv$cluster_label, sc_tsv$category))
  for (col in setdiff(names(sc_rds), c("cluster_label","level","category"))) {
    a <- sc_rds[[col]]; b <- sc_tsv[[col]]
    expect_true(if (is.numeric(a)) eqnum(a, b, 1e-6) else identical(a, b))
  }
})

test_that("fallback selector: 3 branches under vif_correlation=auto", {
  cfgP <- modifyList(default_config(), list(
    ontology_key = "trait_id", levels = c("anthro_class"), primary_level = "anthro_class",
    trait_group = "disease", vif_correlation = "trait_rg", permutation_K = 2000, random_seed = 1))
  route <- read_rds_route(fx("synthetic_ldsc_panel.rds"), cfgP, ".", emit = function(...) invisible(NULL))
  df_rds <- route$df; trait_rg <- route$trait_rg
  ont <- read_ontology(fx("synthetic_panel_ontology.tsv"), "trait_id")
  gP   <- attach_ontology(apply_universe_gate(df_rds, cfgP), ont, "trait_id", c("anthro_class"))
  tids <- unique(gP$trait_id)
  cfgA <- modifyList(cfgP, list(vif_correlation = "auto", vif_coverage_min = 0.5))
  silent <- function(...) invisible(NULL)

  sa <- select_corr_source(gP, cfgA, ".", trait_rg_override = trait_rg, emit = silent)
  expect_true(sa$source == "trait_rg" && sa$coverage >= 0.5)
  sb <- select_corr_source(gP, cfgA, ".", trait_rg_override = identity_corr(tids), emit = silent)
  expect_true(sb$source == "cluster_profile" && approx(sb$coverage, 0))
  g2 <- gP[gP$cluster_label %in% c("C0", "C1"), , drop = FALSE]
  sc <- select_corr_source(g2, cfgA, ".", trait_rg_override = identity_corr(unique(g2$trait_id)), emit = silent)
  expect_identical(sc$source, "identity")
  expect_identical(
    select_corr_source(gP, modifyList(cfgP, list(vif_correlation="cluster_profile")), ".", emit = silent)$source,
    "cluster_profile")

  g2c  <- g2[g2$cluster_label == "C0", , drop = FALSE]; rownames(g2c) <- NULL
  sc_id <- do.call(rbind, score_cluster_level(g2c, "anthro_class", sc$corr, cfgA))
  expect_true(all(approx(sc_id$vif, 1)))
  expect_true(all(approx(sc_id$rho_bar, 0)))
  expect_true(all(approx(sc_id$n_eff, sc_id$n)))
})

test_that("VIF-invariance: AUC / ranks / pooled_rg / coherence identical across corr sources", {
  cfgP <- modifyList(default_config(), list(
    ontology_key = "trait_id", levels = c("anthro_class"), primary_level = "anthro_class",
    trait_group = "disease", vif_correlation = "trait_rg", permutation_K = 2000, random_seed = 1))
  route <- read_rds_route(fx("synthetic_ldsc_panel.rds"), cfgP, ".", emit = function(...) invisible(NULL))
  df_rds <- route$df; trait_rg <- route$trait_rg
  ont <- read_ontology(fx("synthetic_panel_ontology.tsv"), "trait_id")
  gP   <- attach_ontology(apply_universe_gate(df_rds, cfgP), ont, "trait_id", c("anthro_class"))
  tids <- unique(gP$trait_id)
  proxy <- build_trait_profile_corr(gP)
  ident <- identity_corr(tids)
  scores_with <- function(corr_in) {
    set.seed(as.integer(cfgP$random_seed))
    rows <- list()
    for (cl in unique(gP$cluster_label)) {
      gc <- gP[gP$cluster_label == cl, , drop = FALSE]; rownames(gc) <- NULL
      rows <- c(rows, score_cluster_level(gc, "anthro_class", corr_in, cfgP))
    }
    rl <- rank_and_label(do.call(rbind, rows), cfgP)
    rl$ranked[order(rl$ranked$cluster_label, rl$ranked$category), ]
  }
  s_rg <- scores_with(trait_rg); s_pr <- scores_with(proxy); s_id <- scores_with(ident)
  inv <- c("cluster_label","level","category","n","auc_abs","auc_signed","perm_p","q","rank",
           "pooled_rg","coherence","mean_abs_rg","mean_signed_rg","odds_ratio","fisher_p")
  cmp_inv <- function(a, b) all(vapply(inv, function(c) {
    if (is.numeric(a[[c]])) eqnum(a[[c]], b[[c]], 1e-6) else identical(a[[c]], b[[c]])
  }, logical(1)))
  expect_true(cmp_inv(s_rg, s_pr))
  expect_true(cmp_inv(s_rg, s_id))
  expect_true(all(approx(s_id$vif, 1)))
})
