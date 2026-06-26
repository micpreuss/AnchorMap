# uncertainty.R - Phase-6 uncertainty quantification.
#
# Part A (deterministic): a redundancy-aware confidence interval for each competitive AUC.
#   DeLong, DeLong & Clarke-Pearson (1988) nonparametric variance via per-observation placement
#   values (Sun & Xu 2014 fast form from the midranks score.R already holds), VIF-inflated to match
#   the existing vif_z (so the CI shares the redundancy correction), Hanley-McNeil (1982) as the
#   strictly-positive fallback at perfect separation, and a logit transform (Newcombe 2006) that
#   keeps the bounds in (0,1). No RNG -> parity-safe, always on when emit_uncertainty is set.
#
# Part B (Monte-Carlo): a support score for the discrete anchor_shape call. The verdict is a
#   threshold function of (n_sig, margin, focus) crossing several thresholds at once, so a closed
#   form cannot represent the flips; we re-evaluate the shared ruleset (decide_shape) over B draws of
#   the AUCs from their logit-normal sampling distributions and report the fraction recovering the
#   point shape, the full posterior, and a focus CI. Reproducibility is pinned exactly as the z-sweep
#   (per-cluster Mersenne-Twister reseed), and it runs AFTER all perm_p draws so it never perturbs the
#   byte-for-byte perm_p parity. A deterministic jackknife flags single-domain dependence.

# ---- Part A: AUC confidence interval ---------------------------------------

# DeLong (alternative-hypothesis) variance of AUC from the midranks already in hand. NOTE this is the
# variance of the OBSERVED AUC, NOT the H0:AUC=0.5 null variance var0 the vif_p test uses.
# Sun & Xu (2014) placement values: for in-group i with pooled midrank R_i and within-in-group
# midrank Q_i, V10_i = (R_i - Q_i)/n_out; for out-group j, V01_j = 1 - (R_j - Q_j)/n_in.
auc_delong_var <- function(ranks_abs, rv, inmask, n_in, n_out) {
  R_in  <- ranks_abs[inmask]
  R_out <- ranks_abs[!inmask]
  tx <- rank(rv[inmask],  ties.method = "average")        # within-in-group midranks
  ty <- rank(rv[!inmask], ties.method = "average")        # within-out-group midranks
  V10 <- (R_in  - tx) / n_out
  V01 <- 1 - (R_out - ty) / n_in
  S10 <- if (n_in  > 1) stats::var(V10) else 0
  S01 <- if (n_out > 1) stats::var(V01) else 0
  S10 / n_in + S01 / n_out
}

# Hanley-McNeil (1982) closed-form variance; strictly positive for a clamped a in (0,1). Used as the
# fallback when the DeLong variance is 0 / non-finite (perfect separation, single-element group, all
# ties) and as the `auc_ci_method: hanley` option.
auc_hanley_var <- function(a, n_in, n_out) {
  a  <- min(max(a, 1e-6), 1 - 1e-6)
  Q1 <- a / (2 - a)
  Q2 <- 2 * a^2 / (1 + a)
  (a * (1 - a) + (n_in - 1) * (Q1 - a^2) + (n_out - 1) * (Q2 - a^2)) / (n_in * n_out)
}

# Logit-transformed CI from a (VIF-inflated) AUC-scale variance. Returns the AUC-scale SE plus the
# back-transformed bounds, clamped so the point estimate is always contained (a perfect-separation
# AUC=1 then yields a one-sided [lo, 1] interval rather than an empty/invalid one).
auc_ci_logit <- function(auc, var_adj, n_in, n_out, level = 0.95) {
  se  <- sqrt(var_adj)
  eps <- 1 / (2 * n_in * n_out)                            # boundary continuity clamp
  a   <- min(max(auc, eps), 1 - eps)
  se_logit <- se / (a * (1 - a))
  zc  <- stats::qnorm(1 - (1 - level) / 2)
  lo  <- stats::plogis(stats::qlogis(a) - zc * se_logit)
  hi  <- stats::plogis(stats::qlogis(a) + zc * se_logit)
  list(se = se, lo = min(lo, auc), hi = max(hi, auc))
}

# Full Part-A driver: pick the variance source, VIF-inflate, fall back to Hanley-McNeil when the
# chosen variance is degenerate, and return c(se, lo, hi). `vif >= 1` so the interval never shrinks.
auc_ci <- function(ranks_abs, rv, inmask, n_in, n_out, auc, vif, method = "delong", level = 0.95) {
  var_auc <- if (identical(method, "hanley")) auc_hanley_var(auc, n_in, n_out)
             else auc_delong_var(ranks_abs, rv, inmask, n_in, n_out)
  var_adj <- vif * var_auc
  if (!is.finite(var_adj) || var_adj <= 0)                 # perfect separation / degenerate group
    var_adj <- vif * auc_hanley_var(auc, n_in, n_out)
  auc_ci_logit(auc, var_adj, n_in, n_out, level)
}

# ---- Part B: shape confidence ----------------------------------------------

# The shape ruleset as a pure function of the three summary quantities, shared verbatim by the point
# label (anchor_shape) and every MC draw so the two can never diverge. Identical logic to the former
# inline anchor_shape body.
decide_shape <- function(n_sig, margin, focus, cfg) {
  if (n_sig == 0) "weak"
  else if (n_sig == 1 || (!is.na(margin) && margin >= cfg[["shape_margin_sharp"]])) "sharp"
  else if (!is.na(focus) && focus >= cfg[["shape_focus_diffuse"]] &&
           !is.na(margin) && margin < cfg[["shape_margin_diffuse"]]) "diffuse"
  else "focal"
}

# Summary quantities (n_sig, margin, focus) from a RANK-ORDERED AUC vector + its significance mask.
# `auc` must be in rank order (q asc, auc desc) so margin = auc[1]-auc[2] reproduces the point label.
# focus is the inverse-Simpson (Hill 2D) effective number of significant, positively-enriched domains.
shape_summary <- function(auc, sig, cfg) {
  n_sig  <- sum(sig)
  margin <- if (length(auc) > 1) auc[1] - auc[2] else NA_real_
  w <- pmax(auc - 0.5, 0) * sig
  focus <- if (sum(w) > 0) { pp <- w / sum(w); 1 / sum(pp^2) } else NA_real_
  list(n_sig = n_sig, margin = margin, focus = focus)
}

# Significance mask for a (cluster, primary-level) set: q below label_q_max AND AUC at/above
# label_auc_min. NA -> FALSE (mirrors anchor_shape).
.sig_mask <- function(auc, q, cfg) {
  sig <- (q < cfg[["label_q_max"]]) & (auc >= cfg[["label_auc_min"]])
  sig[is.na(sig)] <- FALSE
  sig
}

# Monte-Carlo shape support for one cluster. `sub` = the cluster's primary-level eligible rows,
# rank-sorted, carrying auc_abs, auc_abs_se, q, n. Draws the AUCs B times from independent Gaussian
# sampling distributions on the AUC scale (mean auc_abs, SD = the Part-A auc_abs_se), clamped to
# [0,1], re-evaluates decide_shape per draw, and returns the support fraction for the point shape, the
# posterior over {sharp,focal,diffuse,weak}, and a focus CI. Documented choices/approximations:
#  - AUC-scale (not logit-scale) draws: the logit delta-method SD se/(a(1-a)) degenerates at perfect
#    separation (a->0/1), exploding even when the AUC-scale SE is tiny; the AUC-scale Gaussian with the
#    same SE is robust there and is exactly the sampling distribution the Part-A CI summarises.
#  - q is held FIXED (only AUC perturbed) -> "support conditional on the significance calls".
#  - AUC draws are independent across categories (the true joint shares the complement set).
shape_confidence_mc <- function(sub, cl, sorted_clusters, cfg) {
  a  <- sub[["auc_abs"]]; se <- sub[["auc_abs_se"]]; q <- sub[["q"]]
  sig0 <- .sig_mask(a, q, cfg)
  ss0  <- shape_summary(a, sig0, cfg)
  point_shape <- decide_shape(ss0[["n_sig"]], ss0[["margin"]], ss0[["focus"]], cfg)

  B <- as.integer(cfg[["shape_confidence_B"]])
  na_out <- list(shape_confidence = NA_real_, shape_posterior = NA_character_,
                 anchor_focus_ci_lo = NA_real_, anchor_focus_ci_hi = NA_real_)
  if (is.na(B) || B < 1L || length(a) == 0L) return(na_out)

  # Per-cluster, order- and thread-invariant seed; pin the kind so the parallel sweep's
  # future.seed=TRUE (L'Ecuyer) cannot change the draws (exactly as score_at_z).
  set.seed(as.integer(cfg[["random_seed"]]) + match(cl, sorted_clusters),
           kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")

  k  <- length(a)
  sd_auc <- se
  sd_auc[!is.finite(sd_auc)] <- 0                         # NA/Inf SE -> hold that category fixed

  shapes  <- character(B)
  focuses <- numeric(B)
  for (b in seq_len(B)) {
    astar <- pmin(pmax(a + stats::rnorm(k, 0, sd_auc), 0), 1)
    sigb  <- .sig_mask(astar, q, cfg)
    o     <- order(q, -astar)                             # rank order (q fixed asc, astar desc)
    ss    <- shape_summary(astar[o], sigb[o], cfg)
    shapes[b]  <- decide_shape(ss[["n_sig"]], ss[["margin"]], ss[["focus"]], cfg)
    focuses[b] <- if (is.na(ss[["focus"]])) NA_real_ else ss[["focus"]]
  }

  lv   <- c("sharp", "focal", "diffuse", "weak")
  post <- table(factor(shapes, levels = lv)) / B
  shape_posterior <- paste(sprintf("%s=%.2f", lv, as.numeric(post[lv])), collapse = ";")

  fin    <- focuses[is.finite(focuses)]
  min_n  <- as.numeric(cfg[["shape_confidence_min_sig_frac"]]) * B
  focus_ci <- if (length(fin) >= max(1, min_n))
    stats::quantile(fin, c(0.025, 0.975), names = FALSE, type = 7) else c(NA_real_, NA_real_)

  list(shape_confidence   = round(as.numeric(post[point_shape]), 3),
       shape_posterior    = shape_posterior,
       anchor_focus_ci_lo = if (is.na(focus_ci[1])) NA_real_ else round(focus_ci[1], 2),
       anchor_focus_ci_hi = if (is.na(focus_ci[2])) NA_real_ else round(focus_ci[2], 2))
}

# Deterministic leave-one-domain-out stability: drop each SIGNIFICANT domain in turn, re-rank the
# remainder by (q asc, auc desc), recompute decide_shape; TRUE iff the verdict never changes. With
# fewer than two significant domains there is no alternative-domain dependence to probe, so it is
# trivially stable (a single-domain sharp is flagged by n_sig=1/shape_confidence, not here).
shape_jackknife <- function(sub, cfg) {
  a <- sub[["auc_abs"]]; q <- sub[["q"]]
  sig <- .sig_mask(a, q, cfg)
  ss0 <- shape_summary(a, sig, cfg)
  point <- decide_shape(ss0[["n_sig"]], ss0[["margin"]], ss0[["focus"]], cfg)
  idx <- which(sig)
  if (length(idx) < 2L) return(TRUE)
  for (j in idx) {
    keep <- seq_along(a) != j
    aj <- a[keep]; sigj <- sig[keep]
    o  <- order(q[keep], -aj)
    ssj <- shape_summary(aj[o], sigj[o], cfg)
    if (decide_shape(ssj[["n_sig"]], ssj[["margin"]], ssj[["focus"]], cfg) != point) return(FALSE)
  }
  TRUE
}
