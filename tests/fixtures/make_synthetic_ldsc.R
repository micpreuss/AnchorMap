#!/usr/bin/env Rscript
# tests/fixtures/make_synthetic_ldsc.R — synthetic GenomicSEM ldsc() artifacts with known ground
# truth, for the Phase-2 .rds-ingestion tests (no real ldsc_output.rds + GenomicSEM not installed).
#
# Emits three fixtures + one ontology TSV (all deterministic via set.seed):
#   synthetic_ldsc.rds        5x5: 2 factors (C0, C5_sub0) + 3 traits (BMI, WT, HT). PD; positive
#                             diagonal => S_Stand == cov2cor(S). The plan's primary fixture; used by
#                             the delta-method, partition and analytic checks (q = nrow(V) = 15).
#   synthetic_ldsc_negh2.rds  the same 5x5 with HT's genetic variance set slightly negative, to
#                             exercise the negative-h2 / NA-rg_se guard and the gate-drop path.
#   synthetic_ldsc_panel.rds  11x11: 3 factors (C0, C1, C5_sub0) + 8 panel traits in 2 classes
#                             (anthro x4, cardio x4). Large enough to actually score (N>=6/cluster,
#                             >=3/category) and to reach all three `auto` fallback branches.
#   synthetic_panel_ontology.tsv  trait_id -> anthro_class (+ anchor_eligible) for the panel fixture.
#
# Each $V is a PD sampling-covariance of vech(S) (column-major lower triangle); built with realistic
# small off-diagonals so the delta-method numeric-difference test genuinely exercises the 3x3 path.

suppressPackageStartupMessages(library(data.table))

out_dir <- dirname(normalizePath(sub("^--file=", "",
             grep("^--file=", commandArgs(FALSE), value = TRUE)[1]), mustWork = FALSE))
if (!length(out_dir) || is.na(out_dir)) out_dir <- "tests/fixtures"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Column-major lower-triangle vech length and a PD sampling-covariance V for a k x k S.
vech_len <- function(k) k * (k + 1L) / 2L
make_V <- function(k, base_se = 0.025, seed = 1) {
  set.seed(seed)
  q <- vech_len(k)
  B <- matrix(rnorm(q * 4L, sd = 0.004), q, 4L)   # small correlated noise across vech elements
  V <- B %*% t(B)
  diag(V) <- diag(V) + base_se^2                   # ensure PD + a realistic per-element variance
  V
}

# Build S = diag(sqrt(h2)) %*% R0 %*% diag(sqrt(h2)) so S_Stand == R0 exactly (positive diag).
build_S <- function(names, h2, R0) {
  stopifnot(length(names) == length(h2), all(dim(R0) == length(names)))
  D <- diag(sqrt(h2))
  S <- D %*% R0 %*% D
  dimnames(S) <- list(names, names)
  S
}

# ---- 5x5 primary fixture ---------------------------------------------------
nm5 <- c("C0", "C5_sub0", "BMI", "WT", "HT")
R5 <- matrix(c(
  1.00, 0.10, 0.50, 0.45, 0.20,
  0.10, 1.00, 0.60, 0.55, 0.15,
  0.50, 0.60, 1.00, 0.70, 0.30,
  0.45, 0.55, 0.70, 1.00, 0.35,
  0.20, 0.15, 0.30, 0.35, 1.00), 5, 5, byrow = TRUE)
stopifnot(min(eigen(R5, symmetric = TRUE, only.values = TRUE)$values) > 0)  # PD
h2_5 <- c(0.30, 0.28, 0.32, 0.30, 0.34)
S5 <- build_S(nm5, h2_5, R5)
V5 <- make_V(5L, base_se = 0.025, seed = 11)
dimnames(V5) <- NULL
ldsc5 <- list(S = S5, V = V5, I = diag(5))
saveRDS(ldsc5, file.path(out_dir, "synthetic_ldsc.rds"))

# negative-h2 variant: make HT's genetic variance slightly negative (S unstandardized can be non-PD)
S5n <- S5; S5n["HT", "HT"] <- -0.02
ldsc5n <- list(S = S5n, V = V5, I = diag(5))
saveRDS(ldsc5n, file.path(out_dir, "synthetic_ldsc_negh2.rds"))

# ---- 11x11 panel fixture (scorable) ----------------------------------------
nm11 <- c("C0", "C1", "C5_sub0",
          "BMI", "WT", "HT", "WC",            # anthro
          "T2D", "CAD", "HTN", "AF")          # cardio
set.seed(202606)
A    <- matrix(rnorm(11L * 3L, sd = 0.65), 11L, 3L)        # random 3-factor loadings
Sig  <- A %*% t(A)
diag(Sig) <- diag(Sig) + runif(11L, 0.25, 0.5)            # guarantee PD
R11  <- cov2cor(Sig)
stopifnot(min(eigen(R11, symmetric = TRUE, only.values = TRUE)$values) > 0)
h2_11 <- runif(11L, 0.28, 0.45)                            # all positive => gate-passable
S11 <- build_S(nm11, h2_11, R11)
V11 <- make_V(11L, base_se = 0.020, seed = 23)             # se~0.02 => h2_z ~ 14-22 >> 4 gate
dimnames(V11) <- NULL
ldsc11 <- list(S = S11, V = V11, I = diag(11))
saveRDS(ldsc11, file.path(out_dir, "synthetic_ldsc_panel.rds"))

# ontology for the panel fixture (keyed on trait_id, single anthro_class level; all eligible)
ont <- data.frame(
  trait_id = c("BMI", "WT", "HT", "WC", "T2D", "CAD", "HTN", "AF"),
  anthro_class = c("anthro", "anthro", "anthro", "anthro", "cardio", "cardio", "cardio", "cardio"),
  anchor_eligible = "TRUE", stringsAsFactors = FALSE)
fwrite(ont, file.path(out_dir, "synthetic_panel_ontology.tsv"), sep = "\t", quote = FALSE)

cat(sprintf("wrote fixtures to %s:\n  synthetic_ldsc.rds (5x5, q=%d)\n  synthetic_ldsc_negh2.rds (5x5)\n  synthetic_ldsc_panel.rds (11x11, q=%d)\n  synthetic_panel_ontology.tsv\n",
            out_dir, vech_len(5L), vech_len(11L)))
