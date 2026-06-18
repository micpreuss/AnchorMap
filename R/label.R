# label.R — BH-FDR, ranking, auto-label, and anchor-shape.
# Ports bh_fdr (L285-295), anchor_shape (L301-329) and rank_and_label (L332-374).

# Benjamini-Hochberg q-values (hand-rolled to bit-match the numpy reference; == p.adjust("BH")).
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

# Characterise the anchor profile shape from a cluster's eligible, rank-sorted domains.
anchor_shape <- function(sub, cfg) {
  auc <- sub[["auc_abs"]]; q <- sub[["q"]]
  sig <- (q < cfg[["label_q_max"]]) & (auc >= cfg[["label_auc_min"]])
  sig[is.na(sig)] <- FALSE
  n_sig <- sum(sig)
  margin <- if (length(auc) > 1) auc[1] - auc[2] else NA_real_
  w <- pmax(auc - 0.5, 0) * sig                            # positive enrichment, significant only
  focus <- if (sum(w) > 0) { pp <- w / sum(w); 1 / sum(pp^2) } else NA_real_   # inverse-Simpson
  shape <-
    if (n_sig == 0) "weak"
    else if (n_sig == 1 || (!is.na(margin) && margin >= cfg[["shape_margin_sharp"]])) "sharp"
    else if (!is.na(focus) && focus >= cfg[["shape_focus_diffuse"]] &&
             !is.na(margin) && margin < cfg[["shape_margin_diffuse"]]) "diffuse"
    else "focal"
  list(n_sig_domains = n_sig,
       anchor_margin = if (is.na(margin)) NA_real_ else round(margin, 3),
       anchor_focus  = if (is.na(focus))  NA_real_ else round(focus, 2),
       anchor_shape  = shape)
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
  label_rows <- list()
  for (cl in unique(ranked[["cluster_label"]])) {
    sub <- prim[prim[["cluster_label"]] == cl, , drop = FALSE]
    if (!nrow(sub)) {
      label_rows[[length(label_rows) + 1L]] <- data.frame(
        cluster_label = cl, auto_label = "ambiguous", anchor_shape = "weak",
        anchor_margin = NA_real_, anchor_focus = NA_real_, n_sig_domains = 0L,
        top_auc = NA_real_, top_q = NA_real_, top_pooled_rg = NA_real_,
        top_coherence = NA_real_, profile = "", stringsAsFactors = FALSE)
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
    label_rows[[length(label_rows) + 1L]] <- data.frame(
      cluster_label = cl, auto_label = auto, anchor_shape = shp[["anchor_shape"]],
      anchor_margin = shp[["anchor_margin"]], anchor_focus = shp[["anchor_focus"]],
      n_sig_domains = shp[["n_sig_domains"]],
      top_auc = round(as.numeric(top[["auc_abs"]]), 3), top_q = as.numeric(top[["q"]]),
      top_pooled_rg = round(as.numeric(top[["pooled_rg"]]), 3),
      top_coherence = round(as.numeric(top[["coherence"]]), 3),
      profile = profile, stringsAsFactors = FALSE)
  }
  labels <- do.call(rbind, label_rows); rownames(labels) <- NULL
  list(ranked = ranked, labels = labels)
}
