# plot.R — AnchorMap Phase 4 figures (ggplot2 port of the reference plot_*.py).
#
# Pure re-encoding of the scored TSVs (category_anchor_scores.tsv + cluster_anchor_labels.tsv):
# no statistic is recomputed. Ports three reference scripts into one module —
#   plot_anchors.py             -> fig_lollipops · fig_dotheatmap · fig_scatter
#   plot_specificity.py         -> specificity · distinctive_table · fig_specificity
#   plot_specificity_diagonal.py-> diagonal_column_order · fig_diagonal
# Four channels are load-bearing and kept distinct (AUC and pooled_rg are NOT redundant — they
# diverge at sign-split classes): AUC = x-position/size, signed pooled_rg = diverging colour,
# coherence = alpha, q<q_sig = ring/mask. Headless (cairo/ragg); deterministic (no RNG).

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

# Diverging colour endpoints (match matplotlib RdBu_r / PuOr_r tails, reversed -> low/high).
.RG_LOW   <- "#2166AC"; .RG_HIGH   <- "#B2182B"   # blue (-) -> red (+)   signed rg
.SPEC_LOW <- "#542788"; .SPEC_HIGH <- "#B35806"   # purple (-) -> orange (+)  specificity z
.NA_GREY  <- "grey90"

# ---- ordering ports --------------------------------------------------------
# Natural cluster order: C0, C1, C2_sub0, C2_sub1, … then noise_re0…N last, others last
# (mirror plot_anchors.natural_order; bare C5 sorts before C5_sub0 via sub index -1).
natural_order <- function(labels) {
  labels <- unique(as.character(labels))
  key <- function(l) {
    m <- regmatches(l, regexec("^C(\\d+)(?:_sub(\\d+))?$", l))[[1]]
    if (length(m)) return(c(0, as.integer(m[2]), if (m[3] == "") -1L else as.integer(m[3])))
    m <- regmatches(l, regexec("^noise_re(\\d+)$", l))[[1]]
    if (length(m)) return(c(1, as.integer(m[2]), -1L))
    c(2, 0L, -1L)
  }
  ks <- vapply(labels, key, numeric(3))                       # 3 x n
  labels[order(ks[1, ], ks[2, ], ks[3, ], labels)]
}

# Hierarchical-clustering leaf order of a (rows x cols) matrix; row names returned.
# NaN/NA -> 0, average linkage on euclidean distance (mirror plot_anchors.leaf_order).
leaf_order <- function(mat) {
  rn <- rownames(mat)
  if (is.null(rn) || nrow(mat) < 3) return(rn)
  m <- mat; m[!is.finite(m)] <- 0
  hc <- stats::hclust(stats::dist(m, method = "euclidean"), method = "average")
  rn[hc$order]
}

# Population SD (ddof=0), matching numpy/pandas .std(ddof=0) used in plot_specificity.specificity.
.sd_pop <- function(x) { x <- x[is.finite(x)]; if (length(x) < 1) return(NA_real_); sqrt(mean((x - mean(x))^2)) }

# ---- track loader ----------------------------------------------------------
# Read a scored TSV, coerce the pandas-style "True"/"False" `eligible` to logical, keep only the
# track's level + eligible rows, clamp coherence. Returns the per-track plotting frame + labels.
load_track <- function(t, stage_root) {
  s <- data.table::fread(resolve_path(stage_root, t[["scores"]]), sep = "\t",
                         na.strings = c("", "NA", "NaN"), showProgress = FALSE)
  s[["eligible"]] <- toupper(as.character(s[["eligible"]])) == "TRUE"
  s <- s[s[["level"]] == t[["level"]] & s[["eligible"]], ]
  s[["coherence"]] <- pmin(pmax(fifelse(is.na(s[["coherence"]]), 1, s[["coherence"]]), 0), 1)
  labels <- data.table::fread(resolve_path(stage_root, t[["labels"]]), sep = "\t",
                              na.strings = c("", "NA"), showProgress = FALSE)
  list(name = t[["name"]], level = t[["level"]], s = s, labels = labels)
}

# category (column) order from the rg-profile pivot, hierarchical leaf order
cat_order <- function(s) {
  if (!nrow(s)) return(character(0))
  M <- data.table::dcast(s, category ~ cluster_label, value.var = "pooled_rg")
  mat <- as.matrix(M[, -1, drop = FALSE]); rownames(mat) <- M[["category"]]
  leaf_order(mat)
}

# AUC -> point "size" channel: monotone (clip(auc-0.5,0,0.5)). Exact px is not the contract.
.auc_size <- function(auc) pmin(pmax(auc - 0.5, 0), 0.5)

# ---- diverging scales ------------------------------------------------------
make_rg_fill   <- function(cap) scale_fill_gradient2(low = .RG_LOW, mid = "white", high = .RG_HIGH,
  midpoint = 0, limits = c(-cap, cap), oob = scales::squish, name = "signed pooled rg")
make_rg_colour <- function(cap) scale_colour_gradient2(low = .RG_LOW, mid = "white", high = .RG_HIGH,
  midpoint = 0, limits = c(-cap, cap), oob = scales::squish, name = "signed pooled rg")
make_spec_fill <- function(cap) scale_fill_gradient2(low = .SPEC_LOW, mid = "white", high = .SPEC_HIGH,
  midpoint = 0, limits = c(-cap, cap), oob = scales::squish, na.value = .NA_GREY, name = "specificity z")

# ---- 1. lollipop small-multiples (per track) -------------------------------
# One panel per cluster: x = AUC (0.5->auc), colour = signed pooled_rg, alpha = coherence,
# black ring = q<q_sig, open star = auto-label. Assembled with patchwork at lollipop_ncols.
fig_lollipops <- function(track, row_order, cfg) {
  s <- track[["s"]]; labels <- track[["labels"]]
  q_sig <- cfg[["q_sig"]]; k <- cfg[["top_k"]]; cap <- cfg[["rg_cap"]]
  lab_idx <- labels; data.table::setkey(lab_idx, cluster_label)
  clusters <- row_order[row_order %in% unique(s[["cluster_label"]])]

  panels <- lapply(clusters, function(cl) {
    sub <- s[s[["cluster_label"]] == cl, ]
    sub <- sub[order(sub[["rank"]]), ][seq_len(min(k, nrow(sub))), ]
    auto <- lab_idx[cl][["auto_label"]]
    sub[["y"]]   <- rev(seq_len(nrow(sub)))                       # rank 1 at top
    sub[["a"]]   <- 0.30 + 0.70 * sub[["coherence"]]
    sub[["sig"]] <- sub[["q"]] < q_sig
    sub[["isauto"]] <- sub[["category"]] == auto
    lb <- lab_idx[cl]
    ttl <- sprintf("%s  —  %s [%s]", cl, lb[["auto_label"]], lb[["anchor_shape"]])
    ggplot(sub) +
      geom_vline(xintercept = 0.5, colour = "grey70", linewidth = 0.3) +
      geom_segment(aes(x = 0.5, xend = auc_abs, y = y, yend = y, colour = pooled_rg, alpha = a),
                   linewidth = 1.1) +
      geom_point(aes(x = auc_abs, y = y, colour = pooled_rg, alpha = a), size = 3) +
      geom_point(data = sub[sub[["sig"]], ], aes(x = auc_abs, y = y), shape = 21,
                 fill = NA, colour = "black", size = 3.4, stroke = 0.7) +
      geom_point(data = sub[sub[["isauto"]], ], aes(x = auc_abs, y = y), shape = 8,
                 colour = "black", size = 3, stroke = 0.6) +
      geom_text(aes(x = 0.508, y = y + 0.34, label = .short(category)),
                hjust = 0, vjust = 0.5, size = 2.0) +
      geom_text(aes(x = pmin(auc_abs + 0.008, 0.995), y = y, label = paste0("n", n)),
                hjust = 0, vjust = 0.5, size = 1.7, colour = "grey45") +
      scale_x_continuous(limits = c(0.5, 1.0), name = "AUC (in- vs out-category rank enrichment)") +
      scale_y_continuous(limits = c(0.4, nrow(sub) + 0.8)) +
      make_rg_colour(cap) + scale_alpha_identity() +
      ggtitle(ttl) +
      theme_minimal(base_size = 7) +
      theme(axis.text.y = element_blank(), axis.title.y = element_blank(),
            panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
            plot.title = element_text(size = 7.5), legend.position = "none")
  })
  pw <- patchwork::wrap_plots(panels, ncol = cfg[["lollipop_ncols"]]) +
    patchwork::plot_layout(guides = "collect") +
    patchwork::plot_annotation(title = sprintf(
      "Cluster anchoring — %s track (AUC lead; colour = signed rg, alpha = coherence, ring = q<%.2g)",
      track[["name"]], q_sig))
  pw & theme(legend.position = "right")
}

# ---- 2. dot-heatmap (all tracks side by side) ------------------------------
fig_dotheatmap <- function(tracks, row_order, cfg) {
  q_sig <- cfg[["q_sig"]]; cap <- cfg[["rg_cap"]]
  panels <- lapply(tracks, function(track) {
    s <- data.table::copy(track[["s"]]); labels <- track[["labels"]]
    cats <- cat_order(s)
    s <- s[s[["cluster_label"]] %in% row_order & s[["category"]] %in% cats, ]
    s[["cluster_label"]] <- factor(s[["cluster_label"]], levels = rev(row_order))   # C0 at top
    s[["category"]]      <- factor(s[["category"]], levels = cats)
    s[["sz"]]  <- .auc_size(s[["auc_abs"]])
    s[["sig"]] <- s[["q"]] < q_sig
    stars <- labels[labels[["auto_label"]] %in% cats & labels[["cluster_label"]] %in% row_order, ]
    if (nrow(stars)) {
      stars[["cluster_label"]] <- factor(stars[["cluster_label"]], levels = rev(row_order))
      stars[["category"]]      <- factor(stars[["auto_label"]], levels = cats)
    }
    ggplot(s, aes(x = category, y = cluster_label)) +
      geom_point(aes(size = sz, colour = pooled_rg, alpha = coherence)) +
      geom_point(data = s[s[["sig"]], ], aes(size = sz), shape = 21, fill = NA,
                 colour = "black", stroke = 0.6) +
      { if (nrow(stars)) geom_point(data = stars, shape = 8, colour = "black", size = 2.4,
                                    stroke = 0.6, inherit.aes = TRUE) } +
      make_rg_colour(cap) + scale_alpha_identity() +
      scale_size_area(max_size = 6, limits = c(0, 0.5),
                      breaks = c(0.1, 0.25, 0.4), labels = c("0.60", "0.75", "0.90"),
                      name = "AUC") +
      scale_x_discrete(labels = function(x) .short(x, 22)) +
      labs(title = track[["name"]], x = NULL, y = NULL) +
      theme_minimal(base_size = 8) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6.5),
            panel.grid = element_line(colour = "grey92", linewidth = 0.3))
  })
  pw <- patchwork::wrap_plots(panels, nrow = 1) +
    patchwork::plot_layout(guides = "collect") +
    patchwork::plot_annotation(title = sprintf(
      "Cluster anchoring — dot-heatmap (size = AUC, colour = signed rg, edge = q<%.2g, alpha = coherence, * = auto-label)",
      q_sig))
  pw
}

# ---- 3. AUC vs coherence scatter (diagnostic) ------------------------------
fig_scatter <- function(tracks, cfg) {
  q_sig <- cfg[["q_sig"]]; cap <- cfg[["rg_cap"]]
  panels <- lapply(tracks, function(track) {
    s <- data.table::copy(track[["s"]])
    s[["sz"]]  <- .auc_size(s[["auc_abs"]])
    s[["sig"]] <- s[["q"]] < q_sig
    flag <- s[s[["auc_abs"]] >= 0.65 & s[["coherence"]] <= 0.6, ]
    if (nrow(flag)) flag[["lab"]] <- paste0(flag[["cluster_label"]], "/", .short(flag[["category"]], 14))
    p <- ggplot(s, aes(x = auc_abs, y = coherence)) +
      annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0, ymax = 0.5, fill = "grey93") +
      geom_hline(yintercept = 0.6, linetype = "dashed", colour = "grey60", linewidth = 0.3) +
      geom_point(aes(size = sz, colour = pooled_rg), alpha = 0.85) +
      geom_point(data = s[s[["sig"]], ], aes(size = sz), shape = 21, fill = NA,
                 colour = "black", stroke = 0.4) +
      make_rg_colour(cap) +
      scale_size_area(max_size = 5, limits = c(0, 0.5), guide = "none") +
      scale_x_continuous(limits = c(0.45, 1.0), name = "AUC (enrichment)") +
      scale_y_continuous(limits = c(-0.03, 1.05),
                         name = "coherence (|mean signed rg| / mean |rg|)") +
      labs(title = track[["name"]]) +
      theme_minimal(base_size = 8)
    if (nrow(flag))
      p <- p + ggrepel::geom_text_repel(data = flag, aes(label = lab), size = 2,
                                        colour = "grey25", max.overlaps = 20, seed = 1)
    p
  })
  patchwork::wrap_plots(panels, nrow = 1) +
    patchwork::plot_layout(guides = "collect") +
    patchwork::plot_annotation(title =
      "AUC vs coherence — top-left = high enrichment but sign-split (opposition) classes")
}

# ---- specificity (cross-cluster distinctiveness) ---------------------------
# Within-category z of pooled_rg across clusters + a significance mask (plot_specificity.specificity).
specificity <- function(s, q_sig, rg_floor, min_clusters) {
  M  <- data.table::dcast(s, cluster_label ~ category, value.var = "pooled_rg")
  Q  <- data.table::dcast(s, cluster_label ~ category, value.var = "q")
  cl <- M[["cluster_label"]]
  Mm <- as.matrix(M[, -1, drop = FALSE]); rownames(Mm) <- cl
  Qm <- as.matrix(Q[, -1, drop = FALSE]); rownames(Qm) <- cl
  Qm <- Qm[, colnames(Mm), drop = FALSE]
  mu  <- apply(Mm, 2, function(x) mean(x[is.finite(x)]))
  sdc <- apply(Mm, 2, .sd_pop)
  Z   <- sweep(sweep(Mm, 2, mu, "-"), 2, sdc, "/")
  n_present <- apply(Mm, 2, function(x) sum(is.finite(x)))
  mask <- (Qm < q_sig) & (abs(Mm) >= rg_floor) &
          matrix(n_present >= min_clusters, nrow(Mm), ncol(Mm), byrow = TRUE)
  mask[is.na(mask)] <- FALSE
  list(M = Mm, Z = Z, mask = mask)
}

# Per-cluster most-distinctive (max |z|) significant cell -> distinctive table.
distinctive_table <- function(spec, track) {
  M <- spec[["M"]]; Z <- spec[["Z"]]; mask <- spec[["mask"]]
  rbindlist(lapply(rownames(M), function(cl) {
    z <- Z[cl, ]; m <- mask[cl, ]
    z[!m | !is.finite(z)] <- NA
    ord <- order(abs(z), decreasing = TRUE, na.last = NA)
    if (length(ord)) {
      c1 <- colnames(Z)[ord[1]]
      runner <- if (length(ord) > 1) sprintf("%s(z=%+.1f)", colnames(Z)[ord[2]], z[ord[2]]) else ""
      data.table(track = track, cluster_label = cl, distinctive_category = c1,
                 spec_z = round(z[ord[1]], 2), pooled_rg = round(M[cl, c1], 3), runner_up = runner)
    } else {
      data.table(track = track, cluster_label = cl, distinctive_category = "none",
                 spec_z = NA_real_, pooled_rg = NA_real_, runner_up = "")
    }
  }))
}

# Long frame for a (masked) z-heatmap: only significant cells carry a value, rest NA (grey).
.spec_long <- function(spec, row_order, cats) {
  Z <- spec[["Z"]]; mask <- spec[["mask"]]
  dt <- rbindlist(lapply(row_order, function(cl) {
    if (!cl %in% rownames(Z)) return(NULL)
    data.table(cluster_label = cl, category = colnames(Z), z = Z[cl, ], sig = mask[cl, ])
  }))
  dt <- dt[category %in% cats, ]
  dt[["zmask"]] <- fifelse(dt[["sig"]] & is.finite(dt[["z"]]), dt[["z"]], NA_real_)
  dt[["cluster_label"]] <- factor(dt[["cluster_label"]], levels = rev(row_order))
  dt[["category"]]      <- factor(dt[["category"]], levels = cats)
  dt
}

# Per cluster, the single max-|z| significant cell (box target). Returns a data.table.
.boxed <- function(spec, row_order) {
  Z <- spec[["Z"]]; mask <- spec[["mask"]]
  rbindlist(lapply(row_order, function(cl) {
    if (!cl %in% rownames(Z)) return(NULL)
    z <- Z[cl, ]; m <- mask[cl, ]; z[!m | !is.finite(z)] <- NA
    if (all(is.na(z))) return(NULL)
    dom <- colnames(Z)[which.max(abs(z))]
    data.table(cluster_label = cl, category = dom, z = z[dom])
  }))
}

# ---- 4. specificity heatmap (per track) ------------------------------------
fig_specificity <- function(spec, row_order, name, cap = 2.5) {
  cats <- leaf_order(t(spec[["Z"]]))   # columns clustered on Z (cluster x category) transposed
  cats <- cats[!is.na(cats)]
  dt   <- .spec_long(spec, row_order, cats)
  box  <- .boxed(spec, row_order)
  p <- ggplot(dt, aes(x = category, y = cluster_label)) +
    geom_tile(aes(fill = zmask), colour = "grey85", linewidth = 0.2) +
    make_spec_fill(cap)
  if (nrow(box)) {
    box[["cluster_label"]] <- factor(box[["cluster_label"]], levels = rev(row_order))
    box[["category"]]      <- factor(box[["category"]], levels = cats)
    p <- p + geom_tile(data = box, fill = NA, colour = "black", linewidth = 0.7)
  }
  p + scale_x_discrete(labels = function(x) .short(x, 22)) +
    labs(x = NULL, y = NULL, title = sprintf(
      "Cluster specificity — %s track\nz of pooled rg within each category across clusters (grey = not significant; box = most distinctive)",
      name)) +
    theme_minimal(base_size = 8) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6.5),
          panel.grid = element_blank())
}

# greedy diagonal column order: first appearance of each cluster's boxed domain top->bottom
diagonal_column_order <- function(box, row_order) {
  out <- character(0)
  for (cl in row_order) {
    d <- box[box[["cluster_label"]] == cl, ][["category"]]
    if (length(d) && !d[1] %in% out) out <- c(out, as.character(d[1]))
  }
  out
}

# ---- 5. specificity diagonal (per track) -----------------------------------
fig_diagonal <- function(spec, row_order, name, cap = 2.5) {
  box <- .boxed(spec, row_order)
  if (!nrow(box)) return(NULL)
  cols <- diagonal_column_order(box, row_order)
  if (!length(cols)) return(NULL)
  box[["cluster_label"]] <- factor(box[["cluster_label"]], levels = rev(row_order))
  box[["category"]]      <- factor(as.character(box[["category"]]), levels = cols)
  ggplot(box, aes(x = category, y = cluster_label)) +
    geom_tile(aes(fill = z), colour = "black", linewidth = 0.7) +
    make_spec_fill(cap) +
    scale_x_discrete(labels = function(x) .short(x, 22), drop = FALSE) +
    scale_y_discrete(drop = FALSE) +
    labs(x = NULL, y = NULL, title = sprintf(
      "Cluster specificity (diagonal) — %s track\nmost-distinctive significant domain per cluster; columns ordered for diagonal discrimination",
      name)) +
    theme_minimal(base_size = 8) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
          panel.grid = element_blank())
}

# ---- helpers ---------------------------------------------------------------
.short <- function(s, n = 24) {
  s <- as.character(s)
  ifelse(nchar(s) <= n, s, paste0(substr(s, 1, n - 1), "…"))
}

# Write a ggplot/patchwork object to PNG (+ PDF unless pdf=FALSE). Headless: ragg if present,
# else cairo. Returns the PNG path. The deliberate axis crops (lollipop xlim 0.5-1.0, scatter
# 0.45-1.0) drop depleted off-axis points exactly as the reference matplotlib set_xlim does — that
# "removed N rows outside the scale range" message is intended, so it is muffled here.
save_fig <- function(plot, png_path, width, height, pdf = TRUE) {
  png_dev <- if (requireNamespace("ragg", quietly = TRUE)) ragg::agg_png else
    function(filename, width, height, units, res, ...) grDevices::png(
      filename = filename, width = width, height = height, units = units, res = res, type = "cairo")
  suppressWarnings({
    ggplot2::ggsave(png_path, plot, device = png_dev, width = width, height = height,
                    units = "in", dpi = 200, limitsize = FALSE)
    if (pdf)
      ggplot2::ggsave(sub("\\.png$", ".pdf", png_path), plot, device = grDevices::cairo_pdf,
                      width = width, height = height, units = "in", limitsize = FALSE)
  })
  png_path
}
