# score.R - competitive enrichment per (cluster, level, category).
# Ports auc_from_ranks (L190-192), perm_null_sums (L195-202) and score_cluster_level (L205-282).

# Mann-Whitney U expressed via summed ranks -> AUC = P(in beats out).
auc_from_ranks <- function(ranks_in, n_in, n_out) {
  U <- sum(ranks_in) - n_in * (n_in + 1) / 2
  U / (n_in * n_out)
}

# K label-permutation null sums of n_in ranks drawn without replacement.
# perm_p is reproducible within a fixed RNG seed but Monte-Carlo in nature; deterministic anchoring
# relies on vif_p, not perm_p.
perm_null_sums <- function(rank_vec, n_in, K) {
  N <- length(rank_vec)
  vapply(seq_len(K), function(k) sum(rank_vec[sample.int(N, n_in)]), numeric(1))
}

# Score every category at `level` for one cluster's gated rows `gc`. Returns a list of 1-row frames.
score_cluster_level <- function(gc, level, corr, cfg) {
  N <- nrow(gc)
  min_n <- as.numeric(cfg[["min_category_n"]])
  if (N < min_n * 2) return(list())
  K <- as.integer(cfg[["permutation_K"]])
  emit_unc <- isTRUE(cfg[["emit_uncertainty"]])              # Phase 6: append AUC CI columns
  ci_method <- if (is.null(cfg[["auc_ci_method"]])) "delong" else cfg[["auc_ci_method"]]
  ci_level  <- if (is.null(cfg[["auc_ci_level"]]))  0.95     else as.numeric(cfg[["auc_ci_level"]])

  rv <- as.numeric(gc[[cfg[["rank_variable"]]]])
  ranks_abs    <- rank(rv,             ties.method = "average")
  ranks_signed <- rank(as.numeric(gc[["z"]]), ties.method = "average")

  alpha <- if (isTRUE(cfg[["hit_bonferroni"]])) 0.05 / N else 0.05
  hit <- (gc[["abs_rg"]] >= cfg[["hit_abs_rg"]]) & (gc[["p"]] < alpha)
  hit[is.na(hit)] <- FALSE
  n_hit_total <- sum(hit)

  cats <- as.character(gc[[level]]); cats[is.na(cats)] <- "NA"
  perm_cache <- list()
  rows <- list()
  for (cat in unique(cats)) {
    if (cat == "NA") next
    inmask <- cats == cat
    n_in <- sum(inmask); n_out <- N - n_in
    if (n_in < min_n || n_out < 1) next

    auc_abs    <- auc_from_ranks(ranks_abs[inmask],    n_in, n_out)
    auc_signed <- auc_from_ranks(ranks_signed[inmask], n_in, n_out)

    key <- as.character(n_in)                                  # permutation null cached by in-set size
    if (is.null(perm_cache[[key]])) perm_cache[[key]] <- perm_null_sums(ranks_abs, n_in, K)
    s_obs  <- sum(ranks_abs[inmask])
    perm_p <- (1 + sum(perm_cache[[key]] >= s_obs)) / (K + 1)

    corr_sub <- reindex_corr(corr, gc[["trait_id"]][inmask]) # within-category correlation -> VIF
    rb    <- max(rho_bar(corr_sub), as.numeric(cfg[["vif_min_rho"]]))
    m_eff <- meff_liji(corr_sub)
    vif   <- 1 + (m_eff - 1) * rb
    var0  <- (N + 1) / (12 * n_in * n_out)
    z_un  <- (auc_abs - 0.5) / sqrt(var0)
    vif_z <- if (vif > 0) z_un / sqrt(vif) else z_un
    vif_p <- pnorm(vif_z, lower.tail = FALSE)

    w  <- 1 / gc[["v"]][inmask]; yy <- gc[["y"]][inmask]       # inverse-variance pooled rg (Fisher-z)
    ybar <- sum(w * yy) / sum(w)
    var_ybar <- vif / sum(w)
    pooled_rg <- tanh(ybar)
    ci_lo <- tanh(ybar - 1.96 * sqrt(var_ybar))
    ci_hi <- tanh(ybar + 1.96 * sqrt(var_ybar))
    mean_abs    <- mean(gc[["abs_rg"]][inmask])
    mean_signed <- mean(gc[["rg"]][inmask])
    coherence <- if (mean_abs > 0) abs(mean_signed) / mean_abs else NA_real_

    tp <- sum(hit[inmask]); fn <- n_in - tp                     # Fisher over-representation (2x2)
    fp <- n_hit_total - tp; tn <- n_out - fp
    fisher_p <- stats::fisher.test(matrix(c(tp, fn, fp, tn), nrow = 2, byrow = TRUE),
                                   alternative = "greater")[["p.value"]]
    odds <- (tp * tn) / (fn * fp)                               # scipy sample OR, NOT conditional MLE
    odds_out <- if (is.finite(odds)) round(odds, 3) else Inf

    row <- data.frame(
      cluster_label = gc[["cluster_label"]][1], level = level, category = cat,
      eligible = all(gc[["anchor_eligible"]][inmask]),
      n = n_in, n_eff = round(m_eff, 2), n_hit = tp,
      rho_bar = round(rb, 3), vif = round(vif, 2),
      auc_abs = round(auc_abs, 4), auc_signed = round(auc_signed, 4),
      perm_p = perm_p, vif_z = round(vif_z, 3), vif_p = vif_p,
      pooled_rg = round(pooled_rg, 4),
      pooled_rg_ci_lo = round(ci_lo, 4), pooled_rg_ci_hi = round(ci_hi, 4),
      coherence = if (is.na(coherence)) NA_real_ else round(coherence, 3),
      mean_abs_rg = round(mean_abs, 4), mean_signed_rg = round(mean_signed, 4),
      odds_ratio = odds_out, fisher_p = fisher_p,
      stringsAsFactors = FALSE)

    if (emit_unc) {                                             # Phase 6 Part A: VIF-inflated DeLong CI
      ci <- auc_ci(ranks_abs, rv, inmask, n_in, n_out, auc_abs, vif, ci_method, ci_level)
      row[["auc_abs_se"]]    <- ci[["se"]]                      # full precision (matches perm_p/q style)
      row[["auc_abs_ci_lo"]] <- round(ci[["lo"]], 4)
      row[["auc_abs_ci_hi"]] <- round(ci[["hi"]], 4)
    }
    rows[[length(rows) + 1L]] <- row
  }
  rows
}
