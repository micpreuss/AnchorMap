# redundancy.R — within-category correlation matrix + Li & Ji n_eff + mean pairwise rho.
# Ports load_trait_rg_matrix (L138-156), build_trait_profile_corr (L130-135),
# meff_li_ji (L159-173) and rho_bar (L176-184).
#
# n_eff: the engine uses a Python-matching implementation (eigen + clip negative eigenvalues),
# because the parity gate requires it on possibly non-PD matrices (NaN->0 off-diagonals).
# `poolr::meff(R,"liji")` computes the same formula and is asserted to agree on clean matrices in
# the test suite — but it does not clip, so it is a cross-check, not the parity primary.

# Actual trait x trait genetic-correlation matrix from a FinnGen LDSC --rg summary.
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

# Clean a correlation submatrix the way numpy meff_li_ji does: non-finite -> 0, diag 1, symmetrize.
.clean_corr <- function(R) {
  R[!is.finite(R)] <- 0
  diag(R) <- 1
  (R + t(R)) / 2
}

# Li & Ji (2005) effective number of independent tests. Clips eigenvalues >= 0 to match numpy.
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
