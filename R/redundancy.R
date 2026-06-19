# redundancy.R - within-category correlation matrix + Li & Ji n_eff + mean pairwise rho.
# Ports load_trait_rg_matrix (L138-156), build_trait_profile_corr (L130-135),
# meff_li_ji (L159-173) and rho_bar (L176-184).
#
# n_eff: the engine uses a Python-matching implementation (eigen + clip negative eigenvalues),
# because the parity gate requires it on possibly non-PD matrices (NaN->0 off-diagonals).
# `poolr::meff(R,"liji")` computes the same formula and is asserted to agree on clean matrices in
# the test suite - but it does not clip, so it is a cross-check, not the parity primary.

# Actual trait x trait genetic-correlation matrix from an LDSC --rg summary
# (cols p1,p2,rg,CONVERGED). Symmetric trait_id x trait_id; missing pairs stay NA; diag = 1.
build_trait_rg_matrix <- function(path, traits, require_converged = TRUE) {
  traits <- unique(as.character(traits))
  d <- data.table::fread(path, sep = "\t", select = c("p1","p2","rg","CONVERGED"),
                         colClasses = list(character = c("p1","p2","CONVERGED")),
                         showProgress = FALSE)
  data.table::setDF(d)
  d <- d[d[["p1"]] %in% traits & d[["p2"]] %in% traits, , drop = FALSE]
  if (require_converged) d <- d[toupper(d[["CONVERGED"]]) == "TRUE", , drop = FALSE]
  d[["rg"]] <- pmin(pmax(suppressWarnings(as.numeric(d[["rg"]])), -1), 1)   # LDSC rg can exceed |1|

  idx <- sort(traits)
  M <- matrix(NA_real_, length(idx), length(idx), dimnames = list(idx, idx))
  if (nrow(d)) {
    agg <- aggregate(rg ~ p1 + p2, data = d, FUN = mean)                    # aggfunc="mean" on dups
    M[cbind(match(agg[["p1"]], idx), match(agg[["p2"]], idx))] <- agg[["rg"]]
    tr <- t(M)                                                              # combine_first(m.T): fill NA from transpose
    fill <- is.na(M) & !is.na(tr); M[fill] <- tr[fill]
  }
  diag(M) <- 1
  M
}

# Proxy: correlation of each trait's rg-profile across clusters (trait x trait). Fallback for
# tracks whose trait_id namespace is absent from the trait_rg matrix. min_periods=3 -> mask pairs
# sharing < 3 clusters.
build_trait_profile_corr <- function(g) {
  w <- tapply(g[["rg"]], list(g[["trait_id"]], g[["cluster_label"]]),
              function(x) mean(x, na.rm = TRUE))
  w <- as.matrix(w)
  C <- suppressWarnings(stats::cor(t(w), use = "pairwise.complete.obs"))
  present <- !is.na(w)
  overlap <- present %*% t(present)
  C[overlap < 3] <- NA
  diag(C) <- 1
  C
}

# Reindex a correlation matrix to `ids` x `ids` (missing -> NA), mirroring pandas .reindex.
reindex_corr <- function(corr, ids) {
  ids <- as.character(ids)
  M <- matrix(NA_real_, length(ids), length(ids), dimnames = list(ids, ids))
  common <- intersect(ids, rownames(corr))
  if (length(common)) M[common, common] <- corr[common, common]
  M
}

# Clean a correlation submatrix before the eigen step: non-finite -> 0, diag 1, symmetrize.
.clean_corr <- function(R) {
  R[!is.finite(R)] <- 0
  diag(R) <- 1
  (R + t(R)) / 2
}

# Li & Ji (2005) effective number of independent tests. Clips eigenvalues at 0 (non-PD inputs).
meff_liji <- function(R) {
  m <- nrow(R)
  if (is.null(m) || m < 2) return(as.numeric(if (is.null(m)) 0 else m))
  Rc <- .clean_corr(R)
  ev <- tryCatch(eigen(Rc, symmetric = TRUE, only.values = TRUE)[["values"]],
                 error = function(e) NULL)
  if (is.null(ev)) return(as.numeric(m))
  ev <- pmax(ev, 0)
  sum((ev >= 1) + (ev - floor(ev)))
}

# Mean of the finite off-diagonal (upper-triangle) correlations.
rho_bar <- function(R) {
  m <- nrow(R)
  if (is.null(m) || m < 2) return(0.0)
  vals <- R[upper.tri(R)]
  vals <- vals[is.finite(vals)]
  if (length(vals)) mean(vals) else 0.0
}

# Optional cross-check used by the test suite: poolr::meff on a clean (PD) matrix.
meff_poolr <- function(R) {
  if (!requireNamespace("poolr", quietly = TRUE)) return(NA_real_)
  Rc <- .clean_corr(R)
  tryCatch(as.numeric(poolr::meff(Rc, method = "liji")), error = function(e) NA_real_)
}

# ---- Phase 2: redundancy-source auto-selection -----------------------------

# Gated-trait coverage of a trait x trait correlation matrix: reindex to `tids`, NA the diagonal,
# then the fraction of traits with >= 1 finite off-diagonal correlation. Lifted verbatim from the
# Phase-1 anchor_map.R [vif] block so the explicit trait_rg path is byte-identical.
trait_rg_coverage <- function(corr, tids) {
  sub <- reindex_corr(corr, tids); diag(sub) <- NA
  mean(apply(sub, 1, function(r) any(is.finite(r))))
}

# No-deflation correlation source: unit diagonal, NA off-diagonal. -> rho_bar = 0 -> VIF = 1;
# meff_liji cleans NA->0 (identity) -> n_eff = n. The honest "VIF uncorrected" fallback.
identity_corr <- function(tids) {
  tids <- as.character(unique(tids))
  M <- matrix(NA_real_, length(tids), length(tids), dimnames = list(tids, tids))
  diag(M) <- 1
  M
}

# Select the within-category redundancy source, honouring cfg$vif_correlation:
#   "trait_rg"        -> actual trait x trait LDSC --rg matrix (or `trait_rg_override`); Phase-1 behaviour.
#   "cluster_profile" -> the rg-profile proxy; Phase-1 behaviour.
#   "auto"            -> trait_rg if coverage >= vif_coverage_min, else proxy if >=3 clusters, else
#                        identity (VIF=1) + loud WARN.
# Returns list(corr, source, coverage, reason). VIF affects only vif_p / CI width downstream - never
# the AUC, ranks, pooled_rg point estimate or coherence (asserted in the Phase-2 tests).
select_corr_source <- function(g, cfg, sroot, trait_rg_override = NULL, emit = message) {
  tids       <- unique(g[["trait_id"]])
  n_clusters <- length(unique(g[["cluster_label"]]))
  mode       <- cfg[["vif_correlation"]]
  cov_min    <- as.numeric(cfg[["vif_coverage_min"]])

  build_trait_rg <- function() {
    if (!is.null(trait_rg_override)) return(trait_rg_override)
    mpath <- resolve_path(sroot, cfg[["trait_rg_matrix"]])
    build_trait_rg_matrix(mpath, tids, isTRUE(cfg[["trait_rg_require_converged"]]))
  }

  if (identical(mode, "trait_rg")) {
    corr <- build_trait_rg(); cov <- trait_rg_coverage(corr, tids)
    emit("[vif] source=trait_rg coverage=%.0f%%", 100 * cov)
    if (cov < cov_min) emit("WARN trait_rg coverage %.0f%% < %.0f%% - VIF near-uncorrected", 100 * cov, 100 * cov_min)
    return(list(corr = corr, source = "trait_rg", coverage = cov, reason = "explicit trait_rg"))
  }

  if (identical(mode, "cluster_profile")) {
    corr <- build_trait_profile_corr(g)
    emit("[vif] source=cluster_profile proxy (trait x trait across clusters)")
    return(list(corr = corr, source = "cluster_profile", coverage = NA_real_, reason = "explicit cluster_profile"))
  }

  if (identical(mode, "auto")) {
    trait_rg <- tryCatch(build_trait_rg(), error = function(e) { emit("WARN trait_rg build failed: %s", conditionMessage(e)); NULL })
    cov <- if (!is.null(trait_rg)) trait_rg_coverage(trait_rg, tids) else 0
    if (!is.null(trait_rg) && cov >= cov_min) {
      emit("[vif] source=trait_rg (auto) coverage=%.0f%%", 100 * cov)
      return(list(corr = trait_rg, source = "trait_rg", coverage = cov,
                  reason = sprintf("auto: trait_rg coverage %.0f%% >= %.0f%%", 100 * cov, 100 * cov_min)))
    }
    if (n_clusters >= 3) {
      corr <- build_trait_profile_corr(g)
      emit("[vif] source=cluster_profile (auto fallback) trait_rg coverage=%.0f%% < %.0f%%", 100 * cov, 100 * cov_min)
      return(list(corr = corr, source = "cluster_profile", coverage = cov,
                  reason = sprintf("auto: trait_rg coverage %.0f%% < %.0f%%, %d clusters >=3 -> proxy", 100 * cov, 100 * cov_min, n_clusters)))
    }
    corr <- identity_corr(tids)
    emit("WARN [vif] source=identity (VIF=1 UNCORRECTED): trait_rg coverage=%.0f%% < %.0f%% AND only %d cluster(s) <3",
         100 * cov, 100 * cov_min, n_clusters)
    return(list(corr = corr, source = "identity", coverage = cov,
                reason = sprintf("auto: trait_rg coverage %.0f%% < %.0f%% and %d cluster(s) <3 -> VIF=1", 100 * cov, 100 * cov_min, n_clusters)))
  }

  stop(sprintf("select_corr_source: unknown vif_correlation mode '%s'", mode))
}
