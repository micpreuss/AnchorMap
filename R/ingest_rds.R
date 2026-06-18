# ingest_rds.R — GenomicSEM ldsc() .rds ingestion for AnchorMap (Phase 2, Input C).
#
# Reads a GenomicSEM ldsc() artifact (named list: $S genetic covariance, $V sampling covariance
# of vech(S) in column-major lower-triangle order, $I intercepts), standardizes S -> genetic
# correlation, derives a per-pair rg_se by the *exact* delta method on the 3x3 V-submatrix,
# partitions the variables into cluster factors (-> long-table rows) vs panel traits (-> trait x
# trait redundancy matrix), and assembles the Phase-1 long-table + redundancy matrix so the rest
# of the engine (gate -> redundancy -> score -> label) is route-agnostic.
#
# The parent runs ldsc(..., stand = FALSE) (run_cluster_gpca.R L401-411, to avoid an internal
# crash on negative h2), so the .rds carries the *unstandardized* objects and AnchorMap must
# standardize itself (mirroring run_cluster_gpca.R L414-422).
#
# Depends on resolve_path() (io.R) for rds_trait_meta resolution; sourced after io.R.

# ---- vech indexing ---------------------------------------------------------
# Position map for vech(S) in *column-major lower-triangle* order, matching GenomicSEM's $V:
#   (1,1),(2,1),...,(k,1),(2,2),(3,2),...,(k,2),...,(k,k)
# Returns a named list keyed "i,j" (i >= j) -> 1-based position. Load-bearing: index $V wrong
# and every SE is wrong, so this is unit-tested directly against known positions.
vech_index <- function(k) {
  q   <- k * (k + 1L) / 2L
  idx <- vector("list", q)
  nm  <- character(q)
  p   <- 0L
  for (j in seq_len(k)) for (i in j:k) {
    p <- p + 1L
    nm[p]    <- sprintf("%d,%d", i, j)
    idx[[p]] <- p
  }
  names(idx) <- nm
  idx
}

# ---- reader ----------------------------------------------------------------
# readRDS + shape asserts. $S square + named; $V square with nrow == k(k+1)/2; $I warn-and-carry.
read_ldsc_rds <- function(path) {
  if (!file.exists(path)) stop(sprintf("read_ldsc_rds: file not found: %s", path))
  obj <- readRDS(path)
  if (!is.list(obj) || is.null(obj[["S"]]) || is.null(obj[["V"]]))
    stop("read_ldsc_rds: expected a named list with $S and $V (a GenomicSEM ldsc() artifact)")
  S <- as.matrix(obj[["S"]]); V <- as.matrix(obj[["V"]])
  k <- nrow(S)
  if (ncol(S) != k) stop(sprintf("read_ldsc_rds: $S must be square (got %dx%d)", k, ncol(S)))
  if (is.null(rownames(S)) && is.null(colnames(S)))
    stop("read_ldsc_rds: $S must carry dimnames (the variable names)")
  if (is.null(rownames(S))) rownames(S) <- colnames(S)
  if (is.null(colnames(S))) colnames(S) <- rownames(S)
  q <- k * (k + 1L) / 2L
  if (nrow(V) != q || ncol(V) != q)
    stop(sprintf("read_ldsc_rds: $V must be %dx%d (= k(k+1)/2 for k=%d), got %dx%d",
                 q, q, k, nrow(V), ncol(V)))
  if (is.null(obj[["I"]]))
    message("WARN read_ldsc_rds: $I (LDSC intercepts) absent; carrying on (not consumed by scoring)")
  obj[["S"]] <- S; obj[["V"]] <- V
  obj
}

# ---- standardization (mirror run_cluster_gpca.R L414-422) -------------------
# Clamp negative h2 diagonals to 0 *for the denominator only* so sqrt() is real; never overwrite S
# itself (h2 extraction needs the raw diagonal). Off-diagonals are NOT clipped here (the Phase-1 gate
# clips +/-0.999 for the Fisher-z transform; the matrix builder clips +/-1).
standardize_S <- function(S) {
  s_diag         <- diag(S)
  s_diag_clamped <- pmax(s_diag, 0)
  denom          <- outer(sqrt(s_diag_clamped), sqrt(s_diag_clamped))
  denom[denom == 0] <- 1                       # avoid 0/0 for truly-zero h2 traits
  S_Stand        <- S / denom
  diag(S_Stand)  <- 1
  dimnames(S_Stand) <- dimnames(S)
  S_Stand
}

# ---- delta-method rg_se ----------------------------------------------------
# For r_ij = S_ij / sqrt(S_ii * S_jj), r depends only on (S_ij, S_ii, S_jj), so its variance is
# *exact* from the 3x3 V-submatrix:
#   g    = ( dr/dS_ij , dr/dS_ii , dr/dS_jj )
#        = ( 1/sqrt(S_ii*S_jj) , -r/(2*S_ii) , -r/(2*S_jj) )
#   idx3 = vech positions of (i,j),(i,i),(j,j)   [same order as g]
#   Var  = g' V[idx3,idx3] g ;  rg_se = sqrt(Var)
# Guard: S_ii<=0 or S_jj<=0 (or non-finite) -> rg_se = NA (the gate then drops the row). Diagonal = 0.
rg_se_matrix <- function(S, V) {
  k   <- nrow(S)
  vm  <- vech_index(k)
  pos <- function(a, b) vm[[sprintf("%d,%d", max(a, b), min(a, b))]]
  s   <- diag(S)
  M   <- matrix(NA_real_, k, k, dimnames = dimnames(S))
  diag(M) <- 0
  if (k >= 2) for (j in seq_len(k - 1L)) for (i in (j + 1L):k) {
    Sii <- s[i]; Sjj <- s[j]; Sij <- S[i, j]
    if (!is.finite(Sii) || !is.finite(Sjj) || Sii <= 0 || Sjj <= 0) next
    r    <- Sij / sqrt(Sii * Sjj)
    g    <- c(1 / sqrt(Sii * Sjj), -r / (2 * Sii), -r / (2 * Sjj))
    idx3 <- c(pos(i, j), pos(i, i), pos(j, j))
    v    <- as.numeric(t(g) %*% V[idx3, idx3] %*% g)
    if (is.finite(v) && v >= 0) { M[i, j] <- M[j, i] <- sqrt(v) }
  }
  M
}

# Per-variable h2 SE = sqrt of the V diagonal entry for (t,t). NA where that variance is bad.
h2_se_vector <- function(V, varnames) {
  k  <- length(varnames)
  vm <- vech_index(k)
  se <- vapply(seq_len(k), function(t) {
    v <- V[vm[[sprintf("%d,%d", t, t)]], vm[[sprintf("%d,%d", t, t)]]]
    if (is.finite(v) && v >= 0) sqrt(v) else NA_real_
  }, numeric(1))
  names(se) <- varnames
  se
}

# ---- partition: cluster factors vs panel traits ----------------------------
# Explicit cfg$cluster_factors (character vector) overrides the cfg$cluster_factor_pattern regex.
# Everything not a factor is a panel trait. Error if either set is empty.
partition_S <- function(varnames, cfg) {
  varnames <- as.character(varnames)
  if (!is.null(cfg[["cluster_factors"]])) {
    factors <- intersect(as.character(cfg[["cluster_factors"]]), varnames)
  } else {
    pat <- cfg[["cluster_factor_pattern"]]; if (is.null(pat)) pat <- "^C[0-9]"
    factors <- varnames[grepl(pat, varnames)]
  }
  panel <- setdiff(varnames, factors)
  if (!length(factors))
    stop(sprintf("partition_S: no cluster factors matched (pattern='%s', explicit=%s); names: %s",
                 cfg[["cluster_factor_pattern"]],
                 if (is.null(cfg[["cluster_factors"]])) "NULL" else paste(cfg[["cluster_factors"]], collapse = ","),
                 paste(utils::head(varnames, 10), collapse = ", ")))
  if (!length(panel))
    stop("partition_S: no panel traits left after removing cluster factors (every variable matched).")
  list(factors = factors, panel = panel)
}

# ---- assemble the Phase-1 long-table contract ------------------------------
# One row per (factor f, panel trait t), emitting exactly .LONG_REQUIRED (io.R). status="failed"
# rows survive into df so the gate (not the reader) drops them, matching the TSV route.
rds_to_long <- function(S, S_Stand, rg_se, h2_se, factors, panel, cfg, trait_meta = NULL) {
  okey    <- cfg[["ontology_key"]]
  cat_map <- NULL
  if (!is.null(trait_meta)) {
    if (!all(c("trait_id", "trait_category") %in% names(trait_meta)))
      stop("rds_to_long: rds_trait_meta must have columns trait_id, trait_category")
    cat_map <- setNames(as.character(trait_meta[["trait_category"]]),
                        as.character(trait_meta[["trait_id"]]))
  } else if (identical(okey, "trait_category")) {
    stop(paste("rds_to_long: ontology_key=='trait_category' requires an rds_trait_meta",
               "(trait_id -> trait_category) map to join the disease ontology."))
  }
  tg   <- cfg[["trait_group"]]
  rows <- vector("list", length(factors) * length(panel))
  n    <- 0L
  for (f in factors) for (t in panel) {
    n    <- n + 1L
    rg   <- S_Stand[f, t]; se <- rg_se[f, t]
    h2   <- S[t, t];       h2se <- h2_se[[t]]
    conv <- is.finite(rg) && is.finite(se) && is.finite(h2) && is.finite(h2se)
    tcat <- if (!is.null(cat_map)) unname(cat_map[t]) else NA_character_
    rows[[n]] <- data.frame(
      cluster_label  = f, trait_id = t,
      trait_category = if (length(tcat) && !is.na(tcat)) as.character(tcat) else NA_character_,
      trait_group    = tg,
      rg             = as.numeric(rg), rg_se = as.numeric(se),
      p              = if (is.finite(rg) && is.finite(se) && se > 0) 2 * pnorm(-abs(rg / se)) else NA_real_,
      h2_trait       = as.numeric(h2), h2_trait_se = as.numeric(h2se),
      ldsc_converged = conv, negative_h2 = is.finite(h2) && h2 < 0,
      status         = if (conv) "success" else "failed",
      stringsAsFactors = FALSE)
  }
  do.call(rbind, rows)
}

# Trait x trait redundancy matrix = the panel-trait block of S_Stand, clipped to [-1,1], diag 1.
# Drop-in for build_trait_rg_matrix's output (same shape; reindex_corr/rho_bar/meff_liji consume it).
rds_to_trait_rg <- function(S_Stand, panel) {
  M <- S_Stand[panel, panel, drop = FALSE]
  M <- pmin(pmax(M, -1), 1)
  diag(M) <- 1
  M
}

# ---- driver-facing orchestrator --------------------------------------------
# read -> partition -> standardize -> delta-method -> long-table + trait_rg block.
# Returns list(df, trait_rg, n_factors, n_panel) for anchor_map.R.
read_rds_route <- function(path, cfg, sroot, emit = message) {
  obj      <- read_ldsc_rds(path)
  S        <- obj[["S"]]; V <- obj[["V"]]
  varnames <- rownames(S)
  part     <- partition_S(varnames, cfg)
  S_Stand  <- standardize_S(S)
  rg_se    <- rg_se_matrix(S, V)
  h2_se    <- h2_se_vector(V, varnames)
  trait_meta <- NULL
  if (!is.null(cfg[["rds_trait_meta"]])) {
    mpath <- resolve_path(sroot, cfg[["rds_trait_meta"]])
    trait_meta <- data.table::fread(mpath, sep = "\t", colClasses = "character",
                                    na.strings = c("", "NA"), showProgress = FALSE)
    data.table::setDF(trait_meta)
    emit("[ingest] trait meta %s (%d rows)", mpath, nrow(trait_meta))
  }
  df       <- rds_to_long(S, S_Stand, rg_se, h2_se, part[["factors"]], part[["panel"]], cfg, trait_meta)
  trait_rg <- rds_to_trait_rg(S_Stand, part[["panel"]])
  list(df = df, trait_rg = trait_rg,
       n_factors = length(part[["factors"]]), n_panel = length(part[["panel"]]))
}
