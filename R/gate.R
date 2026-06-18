# gate.R — LDSC reliability gate, per-trait statistics, ontology attach.
# Ports anchor_categories.py:apply_universe_gate (L77-102) and attach_ontology (L105-124).

# Per-row reliability gate + per-trait anchoring statistics. `z` here is the h2-reliability
# threshold (h2_z = h2_trait / h2_trait_se), NOT a trait-relevance cut.
apply_universe_gate <- function(df, cfg) {
  m <- df[["trait_group"]] == cfg[["trait_group"]] &
       df[["status"]] == "success" &
       !is.na(df[["rg"]]) & !is.na(df[["rg_se"]]) & df[["rg_se"]] > 0 &
       !is.na(df[["h2_trait_se"]]) & df[["h2_trait_se"]] > 0
  m[is.na(m)] <- FALSE
  if (isTRUE(cfg[["require_ldsc_converged"]])) m <- m & df[["ldsc_converged"]]
  if (isTRUE(cfg[["drop_negative_h2"]]))       m <- m & !df[["negative_h2"]]
  m[is.na(m)] <- FALSE

  g <- df[m, , drop = FALSE]
  g[["h2_z"]] <- g[["h2_trait"]] / g[["h2_trait_se"]]
  g <- g[g[["h2_z"]] > as.numeric(cfg[["h2_z_threshold"]]), , drop = FALSE]

  rg_c <- pmin(pmax(g[["rg"]], -0.999), 0.999)          # clip for the Fisher-z transform
  g[["abs_rg"]] <- abs(g[["rg"]])
  g[["z"]]      <- g[["rg"]] / g[["rg_se"]]
  g[["abs_z"]]  <- abs(g[["z"]])
  g[["y"]]      <- atanh(rg_c)                           # Fisher-z
  g[["v"]]      <- (g[["rg_se"]]^2) / (1 - rg_c^2)^2     # delta-method variance of y
  rownames(g) <- NULL
  g
}

# Join the ontology onto gated rows on `key`; never clobber existing columns; coerce
# anchor_eligible (missing -> FALSE, mirroring the unmapped-row behaviour); alias the `native` level.
attach_ontology <- function(g, ont, key, levels) {
  bring <- names(ont)[names(ont) == key | !(names(ont) %in% names(g))]
  g <- merge(g, ont[, bring, drop = FALSE], by = key, all.x = TRUE, sort = FALSE)
  if ("anchor_eligible" %in% names(g)) {
    ae <- toupper(as.character(g[["anchor_eligible"]]))
    g[["anchor_eligible"]] <- !is.na(ae) & ae == "TRUE"
  } else {
    g[["anchor_eligible"]] <- TRUE
  }
  if ("native" %in% levels && !("native" %in% names(g))) g[["native"]] <- g[[key]]

  chk <- levels[levels %in% names(g)][1]
  if (!is.na(chk)) {
    unmapped <- sort(unique(as.character(g[[key]][is.na(g[[chk]])])))
    if (length(unmapped))
      message(sprintf("WARN %d %s without ontology mapping on '%s' (e.g. %s)",
                      length(unmapped), key, chk, paste(utils::head(unmapped, 5), collapse = ", ")))
  }
  g
}
