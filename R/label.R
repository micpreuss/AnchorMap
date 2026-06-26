# label.R - BH-FDR, ranking, auto-label, and anchor-shape.
# Ports bh_fdr (L285-295), anchor_shape (L301-329) and rank_and_label (L332-374).

# Benjamini-Hochberg q-values (hand-rolled; equals stats::p.adjust(method = "BH")).
bh_fdr <- function(pvals) {
  p <- as.numeric(pvals); n <- length(p)
  if (n == 0) return(numeric(0))
  o <- order(p)                              # ascending
  q <- numeric(n); prev <- 1.0
  for (i in seq_len(n)) {                     # walk descending p (rankpos = n .. 1)
    rankpos <- n - i + 1L
    j <- o[rankpos]
    prev <- min(prev, p[j] * n / rankpos)
    q[j] <- prev
  }
  q
}

# Characterise the anchor profile shape from a cluster's eligible, rank-sorted domains. The verdict
# itself lives in decide_shape (uncertainty.R) so the point call and the Phase-6 MC draws share one
# ruleset; this just computes the summary quantities and rounds for output.
anchor_shape <- function(sub, cfg) {
  auc <- sub[["auc_abs"]]; q <- sub[["q"]]
  sig <- .sig_mask(auc, q, cfg)
  ss  <- shape_summary(auc, sig, cfg)
  list(n_sig_domains = ss[["n_sig"]],
       anchor_margin = if (is.na(ss[["margin"]])) NA_real_ else round(ss[["margin"]], 3),
       anchor_focus  = if (is.na(ss[["focus"]]))  NA_real_ else round(ss[["focus"]], 2),
       anchor_shape  = decide_shape(ss[["n_sig"]], ss[["margin"]], ss[["focus"]], cfg))
}

# One all-ambiguous label row (used for empty clusters here and padded clusters in the z-sweep). When
# emit_uncertainty is set it carries the Phase-6 columns as NA so every labels frame has one schema.
ambiguous_label_row <- function(cl, cfg = NULL) {
  base <- data.frame(
    cluster_label = cl, auto_label = "ambiguous", anchor_shape = "weak",
    anchor_margin = NA_real_, anchor_focus = NA_real_, n_sig_domains = 0L,
    top_auc = NA_real_, top_q = NA_real_, top_pooled_rg = NA_real_,
    top_coherence = NA_real_, profile = "", stringsAsFactors = FALSE)
  if (isTRUE(cfg[["emit_uncertainty"]])) {
    base[["shape_confidence"]]       <- NA_real_
    base[["anchor_focus_ci_lo"]]     <- NA_real_
    base[["anchor_focus_ci_hi"]]     <- NA_real_
    base[["shape_posterior"]]        <- NA_character_
    base[["shape_jackknife_stable"]] <- NA
  }
  base
}

# BH-FDR within (cluster, level); rank eligible categories by (q asc, AUC desc); then auto-label
# + shape per cluster at the primary level. Returns list(ranked=<df>, labels=<df>).
rank_and_label <- function(scores, cfg) {
  grp <- interaction(scores[["cluster_label"]], scores[["level"]], drop = TRUE, sep = "\r")
  parts <- split(scores, grp)
  ranked_parts <- lapply(parts, function(sub) {
    sub[["q"]] <- bh_fdr(sub[["perm_p"]])
    sub[["rank"]] <- NA_real_
    elig <- which(sub[["eligible"]])
    if (length(elig)) {
      ord <- elig[order(sub[["q"]][elig], -sub[["auc_abs"]][elig])]
      sub[["rank"]][ord] <- seq_along(ord)
    }
    sub
  })
  ranked <- do.call(rbind, ranked_parts); rownames(ranked) <- NULL

  prim <- ranked[ranked[["level"]] == cfg[["primary_level"]] & ranked[["eligible"]], , drop = FALSE]
  emit_unc <- isTRUE(cfg[["emit_uncertainty"]])
  sorted_clusters <- sort(unique(as.character(ranked[["cluster_label"]])))   # MC reseed key (order-invariant)
  label_rows <- list()
  for (cl in unique(ranked[["cluster_label"]])) {
    sub <- prim[prim[["cluster_label"]] == cl, , drop = FALSE]
    if (!nrow(sub)) {
      label_rows[[length(label_rows) + 1L]] <- ambiguous_label_row(cl, cfg)
      next
    }
    sub <- sub[order(sub[["rank"]]), , drop = FALSE]
    shp <- anchor_shape(sub, cfg)
    top <- sub[1, ]
    ok <- top[["q"]] < cfg[["label_q_max"]] && top[["auc_abs"]] >= cfg[["label_auc_min"]] &&
          top[["vif_z"]] > 0 && top[["vif_p"]] < 0.05 && top[["n"]] >= as.numeric(cfg[["min_category_n"]])
    auto <- if (isTRUE(ok)) top[["category"]] else "ambiguous"
    head8 <- utils::head(sub, 8)
    profile <- paste(sprintf("%s (AUC=%.2f q=%.1e rg=%.2f coh=%.2f n=%d)",
                             head8[["category"]], head8[["auc_abs"]], head8[["q"]],
                             head8[["pooled_rg"]], head8[["coherence"]], head8[["n"]]),
                     collapse = "; ")
    row <- data.frame(
      cluster_label = cl, auto_label = auto, anchor_shape = shp[["anchor_shape"]],
      anchor_margin = shp[["anchor_margin"]], anchor_focus = shp[["anchor_focus"]],
      n_sig_domains = shp[["n_sig_domains"]],
      top_auc = round(as.numeric(top[["auc_abs"]]), 3), top_q = as.numeric(top[["q"]]),
      top_pooled_rg = round(as.numeric(top[["pooled_rg"]]), 3),
      top_coherence = round(as.numeric(top[["coherence"]]), 3),
      profile = profile, stringsAsFactors = FALSE)
    if (emit_unc) {                            # Phase 6 Part B: MC shape support + deterministic jackknife
      mc <- shape_confidence_mc(sub, cl, sorted_clusters, cfg)
      row[["shape_confidence"]]       <- mc[["shape_confidence"]]
      row[["anchor_focus_ci_lo"]]     <- mc[["anchor_focus_ci_lo"]]
      row[["anchor_focus_ci_hi"]]     <- mc[["anchor_focus_ci_hi"]]
      row[["shape_posterior"]]        <- mc[["shape_posterior"]]
      row[["shape_jackknife_stable"]] <- shape_jackknife(sub, cfg)
    }
    label_rows[[length(label_rows) + 1L]] <- row
  }
  labels <- do.call(rbind, label_rows); rownames(labels) <- NULL
  list(ranked = ranked, labels = labels)
}
