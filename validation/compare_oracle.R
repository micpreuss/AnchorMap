#!/usr/bin/env Rscript
# validation/compare_oracle.R — cross-language parity check (AnchorMap R vs Python reference).
#   Rscript validation/compare_oracle.R --r-out <scores.tsv> --oracle <scores.tsv> [--K 2000]
#
# Deterministic columns must match (rounded to their stored precision); perm_p/q are checked within
# Monte-Carlo tolerance AND for identical significance; the script exits non-zero on any
# deterministic mismatch.

suppressPackageStartupMessages(library(data.table))

a <- commandArgs(TRUE)
getv <- function(flag, default = NULL) { i <- match(flag, a); if (is.na(i)) default else a[i + 1L] }
r_path  <- getv("--r-out");  o_path <- getv("--oracle");  K <- as.numeric(getv("--K", "2000"))
if (is.null(r_path) || is.null(o_path)) stop("usage: --r-out <file> --oracle <file> [--K N]")

R <- fread(r_path, sep = "\t", na.strings = ""); O <- fread(o_path, sep = "\t", na.strings = "")
key <- c("cluster_label", "level", "category")
setkeyv(R, key); setkeyv(O, key)
M <- merge(O, R, by = key, suffixes = c(".o", ".r"))
if (nrow(M) != nrow(O) || nrow(O) != nrow(R))
  cat(sprintf("NOTE row counts: oracle=%d r=%d joined=%d\n", nrow(O), nrow(R), nrow(M)))

# tolerance by stored rounding precision. odds_ratio = 1.05e-3 because it is rounded to 3 decimals
# and compared across languages: at a round-half-to-even X.XXX5 boundary R 4.x and numpy can land on
# adjacent values (Δ=1e-3) from bit-identical inputs — a benign rounding artifact on the ORA
# diagnostic (not the ranker), not a logic difference (scipy's sample-OR formula is confirmed identical).
tol <- c(n = 0, n_hit = 0, n_eff = 5e-3, vif = 5e-3,
         rho_bar = 5e-4, vif_z = 5e-4, coherence = 5e-4, odds_ratio = 1.05e-3,
         auc_abs = 5e-5, auc_signed = 5e-5, pooled_rg = 5e-5,
         pooled_rg_ci_lo = 5e-5, pooled_rg_ci_hi = 5e-5, mean_abs_rg = 5e-5, mean_signed_rg = 5e-5,
         vif_p = 1e-6, fisher_p = 1e-6)

fails <- 0L
num_cmp <- function(col) {
  o <- suppressWarnings(as.numeric(M[[paste0(col, ".o")]]))
  r <- suppressWarnings(as.numeric(M[[paste0(col, ".r")]]))
  both_inf <- is.infinite(o) & is.infinite(r) & (sign(o) == sign(r))
  na_ok <- is.na(o) & is.na(r)
  d <- abs(o - r); d[both_inf | na_ok] <- 0
  mism <- which(d > tol[[col]] | is.na(d))
  status <- if (length(mism) == 0) "ok " else "FAIL"
  if (length(mism)) {
    fails <<- fails + length(mism)
    cat(sprintf("  %s %-16s maxΔ=%.3g  (%d mismatches; e.g. row %d: o=%s r=%s)\n",
                status, col, suppressWarnings(max(d, na.rm = TRUE)), length(mism), mism[1],
                format(o[mism[1]]), format(r[mism[1]])))
  } else {
    cat(sprintf("  %s %-16s maxΔ=%.3g\n", status, col, max(d)))
  }
}

cat("== deterministic columns ==\n")
for (col in names(tol)) num_cmp(col)

cat("== eligible (logical) ==\n")
elig_o <- as.logical(M[["eligible.o"]]); elig_r <- as.logical(M[["eligible.r"]])
if (all(elig_o == elig_r)) cat("  ok  eligible identical\n") else { fails <- fails + 1L; cat("  FAIL eligible differs\n") }

cat("== rank (where both non-NA) ==\n")
ro <- as.numeric(M[["rank.o"]]); rr <- as.numeric(M[["rank.r"]])
rank_bad <- sum((!is.na(ro) & !is.na(rr) & ro != rr) | (is.na(ro) != is.na(rr)))
if (rank_bad == 0) cat("  ok  rank identical\n") else cat(sprintf("  WARN rank differs in %d rows (perm_p MC noise near ties)\n", rank_bad))

cat("== stochastic columns (perm_p, q): MC tolerance + significance ==\n")
mc_check <- function(col, alpha = 0.05) {
  o <- as.numeric(M[[paste0(col, ".o")]]); r <- as.numeric(M[[paste0(col, ".r")]])
  mc <- 2 * sqrt(pmax(o * (1 - o), 1e-6) / K) + 1 / (K + 1)
  out_tol <- sum(abs(o - r) > pmax(mc, 0.02), na.rm = TRUE)
  sig_flip <- sum((o < alpha) != (r < alpha), na.rm = TRUE)
  cat(sprintf("  %-6s maxΔ=%.3g  | %d beyond MC band | %d significance flips\n",
              col, max(abs(o - r), na.rm = TRUE), out_tol, sig_flip))
}
mc_check("perm_p"); mc_check("q")

cat(sprintf("\n%s — %d deterministic mismatch(es) over %d rows × %d cols\n",
            if (fails == 0) "PASS" else "FAIL", fails, nrow(M), length(tol)))
quit(status = if (fails == 0) 0 else 1)
