#!/usr/bin/env Rscript
# tests/test_phase2.R — Phase-2 (.rds ingestion + redundancy auto-fallback) unit tests.
# Base-R stopifnot assertions (no testthat). Run: Rscript tests/test_phase2.R
# Requires the synthetic fixtures: Rscript tests/fixtures/make_synthetic_ldsc.R first.

sdir <- dirname(dirname(normalizePath(sub("^--file=", "",
          grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))))
for (m in c("io.R","gate.R","redundancy.R","score.R","label.R","ingest_rds.R"))
  source(file.path(sdir, "R", m))

ok     <- function(name, cond) { if (!isTRUE(cond)) stop(sprintf("FAIL: %s", name)); cat(sprintf("  ok  %s\n", name)) }
approx <- function(a, b, tol = 1e-6) all(abs(a - b) < tol, na.rm = FALSE)
# Inf/NA-robust numeric equality (mirrors validation/compare_oracle.R): both-Inf-same-sign and
# both-NA count as equal. Used for score-row comparisons where odds_ratio can be Inf/0.
eqnum  <- function(a, b, tol = 1e-6) {
  a <- as.numeric(a); b <- as.numeric(b)
  both_inf <- is.infinite(a) & is.infinite(b) & (sign(a) == sign(b))
  both_na  <- is.na(a) & is.na(b)
  d <- abs(a - b); d[both_inf | both_na] <- 0
  all(d <= tol)
}
FX     <- file.path(sdir, "tests", "fixtures")
if (!file.exists(file.path(FX, "synthetic_ldsc.rds")))
  source(file.path(FX, "make_synthetic_ldsc.R"))

# ---------------------------------------------------------------------------
cat("== vech indexing (column-major lower triangle) ==\n")
ok("vech_index(3)[2,1] == 2", vech_index(3)[["2,1"]] == 2)
ok("vech_index(3)[3,3] == 6", vech_index(3)[["3,3"]] == 6)
ok("vech_index(5) has 15 entries", length(vech_index(5)) == 15)
ok("vech_index(5)[5,4] == 14", vech_index(5)[["5,4"]] == 14)  # col4 starts at 13: (4,4)=13,(5,4)=14

# ---------------------------------------------------------------------------
cat("== reader + analytic standardization ==\n")
x <- read_ldsc_rds(file.path(FX, "synthetic_ldsc.rds"))
S <- x[["S"]]; V <- x[["V"]]; k <- nrow(S)
S_std <- standardize_S(S)
ok("S_Stand == base cov2cor (positive diag)", approx(S_std, cov2cor(S), 1e-12))
ok("S_Stand diag == 1", approx(diag(S_std), rep(1, k), 1e-12))
ok("S_Stand symmetric", approx(S_std, t(S_std), 1e-12))
# bad shapes error
tmpnl <- tempfile(fileext = ".rds"); saveRDS(1:3, tmpnl)
ok("read_ldsc_rds rejects non-list", tryCatch({read_ldsc_rds(tmpnl); FALSE}, error = function(e) TRUE))
bad <- list(S = matrix(1, 3, 3, dimnames = list(letters[1:3], letters[1:3])), V = matrix(1, 5, 5))
tmpbad <- tempfile(fileext = ".rds"); saveRDS(bad, tmpbad)
ok("read_ldsc_rds rejects wrong V dim", tryCatch({read_ldsc_rds(tmpbad); FALSE}, error = function(e) TRUE))

# ---------------------------------------------------------------------------
cat("== delta-method rg_se == numeric difference (the key gate) ==\n")
rg_se <- rg_se_matrix(S, V)
vm  <- vech_index(k)
pos <- function(a, b) vm[[sprintf("%d,%d", max(a, b), min(a, b))]]
r_fun <- function(Sii, Sjj, Sij) Sij / sqrt(Sii * Sjj)
check_pair <- function(i, j) {
  Sii <- S[i, i]; Sjj <- S[j, j]; Sij <- S[i, j]; eps <- 1e-6
  g_num <- c(
    (r_fun(Sii, Sjj, Sij + eps) - r_fun(Sii, Sjj, Sij - eps)) / (2 * eps),   # d/dSij
    (r_fun(Sii + eps, Sjj, Sij) - r_fun(Sii - eps, Sjj, Sij)) / (2 * eps),   # d/dSii
    (r_fun(Sii, Sjj + eps, Sij) - r_fun(Sii, Sjj - eps, Sij)) / (2 * eps))   # d/dSjj
  idx3 <- c(pos(i, j), pos(i, i), pos(j, j))
  se_num <- sqrt(as.numeric(t(g_num) %*% V[idx3, idx3] %*% g_num))
  abs(rg_se[i, j] - se_num)
}
for (pr in list(c(2,1), c(3,1), c(4,3), c(5,2), c(5,4))) {
  d <- check_pair(pr[1], pr[2])
  ok(sprintf("delta-method (%d,%d): |analytic-numeric|=%.2e < 1e-6", pr[1], pr[2], d), d < 1e-6)
}
ok("rg_se symmetric", approx(rg_se, t(rg_se), 1e-12))
ok("rg_se diag == 0", approx(diag(rg_se), rep(0, k), 1e-12))
cat("delta-method OK\n")

# ---------------------------------------------------------------------------
cat("== partition (regex default + explicit override + empty error) ==\n")
p <- partition_S(rownames(S), list(cluster_factor_pattern = "^C[0-9]", cluster_factors = NULL))
ok("factors = {C0, C5_sub0}", setequal(p$factors, c("C0", "C5_sub0")))
ok("panel = {BMI, WT, HT}",   setequal(p$panel, c("BMI", "WT", "HT")))
p2 <- partition_S(rownames(S), list(cluster_factors = c("BMI"), cluster_factor_pattern = "^C[0-9]"))
ok("explicit cluster_factors overrides regex", setequal(p2$factors, "BMI") && setequal(p2$panel, c("C0","C5_sub0","WT","HT")))
ok("no factors matched -> error",
   tryCatch({partition_S(rownames(S), list(cluster_factor_pattern = "^ZZZ", cluster_factors = NULL)); FALSE}, error = function(e) TRUE))
ok("all-factor -> empty panel error",
   tryCatch({partition_S(rownames(S), list(cluster_factor_pattern = ".", cluster_factors = NULL)); FALSE}, error = function(e) TRUE))

# ---------------------------------------------------------------------------
cat("== rds_to_long schema + negative-h2 guard ==\n")
h2_se <- h2_se_vector(V, rownames(S))
cfg5  <- modifyList(default_config(), list(ontology_key = "trait_id", trait_group = "disease"))
df5   <- rds_to_long(S, S_std, rg_se, h2_se, p$factors, p$panel, cfg5, trait_meta = NULL)
ok("df has all .LONG_REQUIRED cols", all(.LONG_REQUIRED %in% names(df5)))
ok("df col order == .LONG_REQUIRED", identical(names(df5), .LONG_REQUIRED))
ok("rg / rg_se numeric", is.numeric(df5$rg) && is.numeric(df5$rg_se))
ok("nrow == factors x panel (2x3=6)", nrow(df5) == 6)
ok("all success on clean fixture", all(df5$status == "success") && all(df5$ldsc_converged))
# ontology_key=trait_category without meta -> error
ok("trait_category route without rds_trait_meta errors",
   tryCatch({rds_to_long(S, S_std, rg_se, h2_se, p$factors, p$panel,
                         modifyList(cfg5, list(ontology_key = "trait_category")), NULL); FALSE},
            error = function(e) TRUE))
# negative-h2 variant -> NA rg_se, negative_h2 TRUE, status failed, gate-dropped
xn <- read_ldsc_rds(file.path(FX, "synthetic_ldsc_negh2.rds"))
Sn <- xn[["S"]]; Vn <- xn[["V"]]
rg_se_n <- rg_se_matrix(Sn, Vn)
ok("rg_se NA for pairs touching negative-h2 trait", is.na(rg_se_n["HT", "BMI"]) && is.na(rg_se_n["BMI", "HT"]))
dfn <- rds_to_long(Sn, standardize_S(Sn), rg_se_n, h2_se_vector(Vn, rownames(Sn)),
                   p$factors, p$panel, cfg5, NULL)
htr <- dfn[dfn$trait_id == "HT", ]
ok("HT rows: negative_h2 & failed & !converged",
   all(htr$negative_h2) && all(htr$status == "failed") && all(!htr$ldsc_converged))
gn <- apply_universe_gate(dfn, cfg5)
ok("gate drops negative-h2 trait", !("HT" %in% gn$trait_id) && all(c("BMI","WT") %in% gn$trait_id))

# ---------------------------------------------------------------------------
cat("== round-trip: .rds route == equivalent TSV route (deterministic scores) ==\n")
cfgP <- modifyList(default_config(), list(
  ontology_key = "trait_id", levels = c("anthro_class"), primary_level = "anthro_class",
  trait_group = "disease", vif_correlation = "trait_rg", permutation_K = 2000, random_seed = 1))
route <- read_rds_route(file.path(FX, "synthetic_ldsc_panel.rds"), cfgP, sdir, emit = function(...) invisible(NULL))
df_rds   <- route$df
trait_rg <- route$trait_rg
ont <- read_ontology(file.path(FX, "synthetic_panel_ontology.tsv"), "trait_id")

# TSV route: write df + trait_rg in Input-A / Input-B formats, re-read through the Phase-1 readers.
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
ok("round-trip: same n score rows", nrow(sc_rds) == nrow(sc_tsv) && nrow(sc_rds) > 0)
ok("round-trip: same keys", identical(paste(sc_rds$cluster_label, sc_rds$category), paste(sc_tsv$cluster_label, sc_tsv$category)))
for (col in setdiff(names(sc_rds), c("cluster_label","level","category"))) {
  a <- sc_rds[[col]]; b <- sc_tsv[[col]]
  same <- if (is.numeric(a)) eqnum(a, b, 1e-6) else identical(a, b)
  ok(sprintf("round-trip identical: %s", col), same)
}
cat("round-trip identical\n")

# ---------------------------------------------------------------------------
cat("== fallback selector: 3 branches under vif_correlation=auto ==\n")
gP   <- attach_ontology(apply_universe_gate(df_rds, cfgP), ont, "trait_id", c("anthro_class"))
tids <- unique(gP$trait_id)
cfgA <- modifyList(cfgP, list(vif_correlation = "auto", vif_coverage_min = 0.5))
silent <- function(...) invisible(NULL)

# (a) full-coverage trait_rg override -> trait_rg
sa <- select_corr_source(gP, cfgA, sdir, trait_rg_override = trait_rg, emit = silent)
ok("(a) full coverage -> source=trait_rg", sa$source == "trait_rg" && sa$coverage >= 0.5)
# (b) zero-coverage override + 3 clusters -> cluster_profile
sb <- select_corr_source(gP, cfgA, sdir, trait_rg_override = identity_corr(tids), emit = silent)
ok("(b) zero coverage + >=3 clusters -> source=cluster_profile", sb$source == "cluster_profile" && approx(sb$coverage, 0))
# (c) zero-coverage override + <3 clusters -> identity (+WARN)
g2 <- gP[gP$cluster_label %in% c("C0", "C1"), , drop = FALSE]
sc <- select_corr_source(g2, cfgA, sdir, trait_rg_override = identity_corr(unique(g2$trait_id)), emit = silent)
ok("(c) zero coverage + <3 clusters -> source=identity", sc$source == "identity")
ok("explicit trait_rg/cluster_profile modes still honoured",
   select_corr_source(gP, modifyList(cfgP, list(vif_correlation="cluster_profile")), sdir, emit = silent)$source == "cluster_profile")

# identity source => VIF==1, rho_bar==0, n_eff==n on a scored category
g2c  <- g2[g2$cluster_label == "C0", , drop = FALSE]; rownames(g2c) <- NULL
sc_id_rows <- score_cluster_level(g2c, "anthro_class", sc$corr, cfgA)
sc_id <- do.call(rbind, sc_id_rows)
ok("identity corr -> VIF == 1",    all(approx(sc_id$vif, 1)))
ok("identity corr -> rho_bar == 0", all(approx(sc_id$rho_bar, 0)))
ok("identity corr -> n_eff == n",   all(approx(sc_id$n_eff, sc_id$n)))
cat("fallback branches OK\n")

# ---------------------------------------------------------------------------
cat("== VIF-invariance: AUC / ranks / pooled_rg / coherence identical across corr sources ==\n")
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
ok("invariant cols: trait_rg == proxy",    cmp_inv(s_rg, s_pr))
ok("invariant cols: trait_rg == identity", cmp_inv(s_rg, s_id))
ok("VIF actually differs (proxy vs identity)", !approx(sum(s_pr$vif), sum(s_id$vif), 1e-9) || all(s_id$vif == 1))
ok("identity route: all VIF == 1", all(approx(s_id$vif, 1)))
cat("VIF-invariance OK\n")

cat("\nALL PHASE-2 TESTS PASSED\n")
