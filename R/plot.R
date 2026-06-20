# plot.R - AnchorMap figures (ggplot2).
#
# Pure re-encoding of the scored TSVs (category_anchor_scores.tsv + cluster_anchor_labels.tsv):
# no statistic is recomputed. Three figure families in one module -
#   fig_lollipops / fig_dotheatmap / fig_scatter      (cluster-wise anchoring)
#   specificity / distinctive_table / fig_specificity (cross-cluster distinctiveness)
#   diagonal_column_order / fig_diagonal              (the diagonal reduction)
# Four channels are load-bearing and kept distinct (AUC and pooled_rg are NOT redundant - they
# diverge at sign-split classes): AUC = x-position/size, signed pooled_rg = diverging colour,
# coherence = alpha, q<q_sig = ring/mask. Headless (cairo/ragg); deterministic (no RNG).

# Diverging colour endpoints (RdBu_r / PuOr_r tails, reversed -> low/high).
.RG_LOW   <- "#2166AC"; .RG_HIGH   <- "#B2182B"   # blue (-) -> red (+)   signed rg
.SPEC_LOW <- "#542788"; .SPEC_HIGH <- "#B35806"   # purple (-) -> orange (+)  specificity z
.NA_GREY  <- "grey90"

# ---- ordering ports --------------------------------------------------------
# Natural cluster order: C0, C1, C2_sub0, C2_sub1, ... then noise_re0...N last, others last
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

# Population standard deviation (divide by N, not N-1) used in the cross-cluster specificity z.
.sd_pop <- function(x) { x <- x[is.finite(x)]; if (length(x) < 1) return(NA_real_); sqrt(mean((x - mean(x))^2)) }

# ---- track loader ----------------------------------------------------------
# Read a scored TSV, coerce the pandas-style "True"/"False" `eligible` to logical, keep only the
# track's level + eligible rows, clamp coherence. Returns the per-track plotting frame + labels.
load_track <- function(t, stage_root, in_dir = NULL) {
  # When in_dir is set, read the engine outputs from there (keeping each track's basename) so the
  # figures can target any engine --out-dir; otherwise resolve the config's paths against stage_root.
  loc <- function(p) if (!is.null(in_dir)) file.path(in_dir, basename(p)) else resolve_path(stage_root, p)
  track_name <- if (is.null(t[["name"]])) "<unnamed>" else as.character(t[["name"]])
  scores_path <- loc(t[["scores"]]); labels_path <- loc(t[["labels"]])
  for (path in c(scores_path, labels_path))
    if (!file.exists(path))
      stop(sprintf("track '%s': file not found: %s", track_name, path), call. = FALSE)

  s <- data.table::fread(file = scores_path, sep = "\t",
                         na.strings = c("", "NA", "NaN"), showProgress = FALSE)
  required_scores <- c("level", "eligible", "category", "cluster_label", "pooled_rg", "q",
                       "coherence", "auc_abs", "n", "rank")
  missing_scores <- setdiff(required_scores, names(s))
  if (length(missing_scores))
    stop(sprintf("track '%s': scores file missing required column(s): %s", track_name,
                 paste(missing_scores, collapse = ", ")), call. = FALSE)

  levels_present <- sort(unique(as.character(s[["level"]][!is.na(s[["level"]])])))
  s[["eligible"]] <- toupper(as.character(s[["eligible"]])) == "TRUE"
  s <- s[s[["level"]] == t[["level"]] & s[["eligible"]], ]
  if (!nrow(s))
    stop(sprintf("track '%s': no eligible rows at level '%s' in %s (levels present: %s)",
                 track_name, t[["level"]], scores_path,
                 if (length(levels_present)) paste(levels_present, collapse = ", ") else "<none>"),
         call. = FALSE)
  s[["coherence"]] <- pmin(pmax(fifelse(is.na(s[["coherence"]]), 1, s[["coherence"]]), 0), 1)
  labels <- data.table::fread(file = labels_path, sep = "\t",
                              na.strings = c("", "NA"), showProgress = FALSE)
  missing_labels <- setdiff(c("cluster_label", "auto_label", "anchor_shape"), names(labels))
  if (length(missing_labels))
    stop(sprintf("track '%s': labels file missing required column(s): %s", track_name,
                 paste(missing_labels, collapse = ", ")), call. = FALSE)
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
    ttl <- sprintf("%s  -  %s [%s]", cl, lb[["auto_label"]], lb[["anchor_shape"]])
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
      "Cluster anchoring - %s track (AUC lead; colour = signed rg, alpha = coherence, ring = q<%.2g)",
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
      "Cluster anchoring - dot-heatmap (size = AUC, colour = signed rg, edge = q<%.2g, alpha = coherence, * = auto-label)",
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
      "AUC vs coherence - top-left = high enrichment but sign-split (opposition) classes")
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
      "Cluster specificity - %s track\nz of pooled rg within each category across clusters (grey = not significant; box = most distinctive)",
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
      "Cluster specificity (diagonal) - %s track\nmost-distinctive significant domain per cluster; columns ordered for diagonal discrimination",
      name)) +
    theme_minimal(base_size = 8) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
          panel.grid = element_blank())
}

# ---- helpers ---------------------------------------------------------------
.short <- function(s, n = 24) {
  s <- as.character(s)
  ifelse(nchar(s) <= n, s, paste0(substr(s, 1, n - 1), "..."))
}

# Write a ggplot/patchwork object to PNG (+ PDF unless pdf=FALSE). Headless: ragg if present,
# else cairo. Returns every written path (PNG, and PDF unless pdf=FALSE) so callers can record the
# full artifact set. The deliberate axis crops (lollipop xlim 0.5-1.0, scatter 0.45-1.0) drop
# depleted off-axis points exactly as the reference matplotlib set_xlim does - that "removed N rows
# outside the scale range" message is intended, so it is muffled here.
save_fig <- function(plot, png_path, width, height, pdf = TRUE) {
  png_dev <- if (requireNamespace("ragg", quietly = TRUE)) ragg::agg_png else
    function(filename, width, height, units, res, ...) grDevices::png(
      filename = filename, width = width, height = height, units = units, res = res, type = "cairo")
  pdf_path <- sub("\\.png$", ".pdf", png_path)
  suppressWarnings({
    ggplot2::ggsave(png_path, plot, device = png_dev, width = width, height = height,
                    units = "in", dpi = 200, limitsize = FALSE)
    if (pdf)
      ggplot2::ggsave(pdf_path, plot, device = grDevices::cairo_pdf,
                      width = width, height = height, units = "in", limitsize = FALSE)
  })
  if (pdf) c(png_path, pdf_path) else png_path
}

# ---- driver ----------------------------------------------------------------
#' Render the AnchorMap figures from scored TSVs
#'
#' Reads a plot YAML config (`out_dir` + `tracks: [{name, level, scores, labels}]`) and writes the
#' cluster-anchoring lollipops, dot-heatmap, AUC-vs-coherence scatter, and cross-cluster specificity
#' heatmap + diagonal (PNG + PDF) plus `cluster_distinctive_categories.tsv` into `out_dir`. Does not
#' re-run the engine - it only re-encodes the scored TSVs.
#'
#' @param config_path Path to a plot YAML config (or a bare shipped-config name).
#' @param q_sig,rg_floor,min_clusters Optional overrides of the significance gate.
#' @param out_dir Optional output-directory override (else `cfg$out_dir`).
#' @param in_dir Optional single-track input-directory override: read the track's `scores`/`labels`
#'   from this directory (by basename) instead of the config paths. Multi-track configs must set
#'   each track's `scores`/`labels` paths in the config.
#' @return Invisibly, the character vector of written file paths.
#' @export
run_plots <- function(config_path, q_sig = NULL, rg_floor = NULL, min_clusters = NULL,
                      out_dir = NULL, in_dir = NULL) {
  options(bitmapType = "cairo")
  config_path <- resolve_config_path(config_path)
  cfg   <- yaml::read_yaml(config_path)
  sroot <- stage_root_of(config_path)
  defaults <- list(top_k = 8, lollipop_ncols = 3, scatter = TRUE, rg_cap = 0.55, q_sig = 0.05,
                   spec_rg_floor = 0.10, spec_min_clusters = 5)
  for (k in names(defaults)) if (is.null(cfg[[k]])) cfg[[k]] <- defaults[[k]]
  if (!is.null(q_sig))        cfg[["q_sig"]] <- q_sig
  if (!is.null(rg_floor))     cfg[["spec_rg_floor"]] <- rg_floor
  if (!is.null(min_clusters)) cfg[["spec_min_clusters"]] <- min_clusters

  in_dir <- if (!is.null(in_dir)) .abs_cwd(in_dir) else NULL
  if (!is.null(in_dir) && length(cfg[["tracks"]]) > 1L) {
    track_names <- vapply(cfg[["tracks"]], function(t)
      if (is.null(t[["name"]])) "<unnamed>" else as.character(t[["name"]]), character(1))
    stop(sprintf("--in-dir is a single-track convenience; config has %d tracks (%s). Set each ",
                 length(track_names), paste(track_names, collapse = ", ")),
         "track's scores/labels paths in the config instead.", call. = FALSE)
  }

  out_dir <- if (!is.null(out_dir)) .abs_cwd(out_dir) else resolve_path(sroot, cfg[["out_dir"]])
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  written <- character(0)
  emit <- function(p) { for (q in p) message(sprintf("[write] %s", q)); written <<- c(written, p) }

  tracks <- lapply(cfg[["tracks"]], load_track, stage_root = sroot, in_dir = in_dir)
  row_order <- natural_order(unlist(lapply(tracks, function(tr) tr[["s"]][["cluster_label"]])))

  # 1. lollipops (per track)
  for (tr in tracks) {
    nc   <- min(cfg[["lollipop_ncols"]], length(unique(tr[["s"]][["cluster_label"]])))
    n_cl <- length(intersect(row_order, unique(tr[["s"]][["cluster_label"]])))
    nrw  <- ceiling(n_cl / max(nc, 1))
    p <- fig_lollipops(tr, row_order, cfg)
    png <- file.path(out_dir, sprintf("anchor_lollipops_%s.png", tr[["name"]]))
    emit(save_fig(p, png, width = nc * 3.6 + 1.5, height = nrw * 2.3 + 0.6))
  }

  # 2. dot-heatmap (combined)
  totcats <- sum(vapply(tracks, function(tr) length(unique(tr[["s"]][["category"]])), integer(1)))
  png <- file.path(out_dir, "anchor_dotheatmap.png")
  emit(save_fig(fig_dotheatmap(tracks, row_order, cfg), png,
                width = totcats * 0.34 + 3.5, height = length(row_order) * 0.34 + 2.2))

  # 3. AUC-vs-coherence scatter (combined)
  if (isTRUE(cfg[["scatter"]])) {
    png <- file.path(out_dir, "anchor_auc_coherence.png")
    emit(save_fig(fig_scatter(tracks, cfg), png, width = length(tracks) * 4.6, height = 4.4))
  }

  # 4-5. specificity heatmap + diagonal (per track) + distinctive table
  dist_all <- list()
  for (tr in tracks) {
    spec <- specificity(tr[["s"]], cfg[["q_sig"]], cfg[["spec_rg_floor"]], cfg[["spec_min_clusters"]])
    ro   <- natural_order(rownames(spec[["M"]]))
    ncat <- ncol(spec[["Z"]])
    png  <- file.path(out_dir, sprintf("anchor_specificity_%s.png", tr[["name"]]))
    emit(save_fig(fig_specificity(spec, ro, tr[["name"]]), png,
                  width = ncat * 0.42 + 3.5, height = length(ro) * 0.34 + 2.0))
    pd <- fig_diagonal(spec, ro, tr[["name"]])
    if (is.null(pd)) {
      message(sprintf("[skip] anchor_specificity_diagonal_%s - no significant distinctive cells",
                      tr[["name"]]))
    } else {
      ncols <- length(diagonal_column_order(.boxed(spec, ro), ro))
      png <- file.path(out_dir, sprintf("anchor_specificity_diagonal_%s.png", tr[["name"]]))
      emit(save_fig(pd, png, width = ncols * 0.6 + 3.5, height = length(ro) * 0.34 + 2.0))
    }
    dist_all[[length(dist_all) + 1L]] <- distinctive_table(spec, tr[["name"]])
  }
  dist <- data.table::rbindlist(dist_all)
  dist_path <- file.path(out_dir, "cluster_distinctive_categories.tsv")
  data.table::fwrite(dist, dist_path, sep = "\t", na = "", quote = FALSE)
  emit(dist_path)

  message(sprintf("FINISHED ok | %d figure files + 1 table | out_dir: %s",
                  length(written) - 1L, out_dir))
  invisible(written)
}
