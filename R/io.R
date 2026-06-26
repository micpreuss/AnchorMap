# io.R - config + standardized-TSV readers + schema asserts.
# The config + standardized-TSV / ontology readers. The YAML schema is stable, so configs are
# portable across runs.

# ---- config ----------------------------------------------------------------
# Engine defaults; a YAML config overlays these.
default_config <- function() {
  list(
    trait_group = "disease", require_ldsc_converged = TRUE, drop_negative_h2 = TRUE,
    h2_z_threshold = 4.0, levels = c("native","domain","icd_chapter"),
    primary_level = "domain", ontology_key = "trait_category", min_category_n = 3,
    rank_variable = "abs_z", permutation_K = 2000, random_seed = 1, vif_min_rho = 0.0,
    vif_correlation = "cluster_profile", trait_rg_require_converged = TRUE,
    hit_abs_rg = 0.2, hit_bonferroni = TRUE, label_q_max = 0.05, label_auc_min = 0.60,
    shape_margin_sharp = 0.10, shape_margin_diffuse = 0.05, shape_focus_diffuse = 3.0,
    # ---- Phase 2 (additive; do NOT change vif_correlation's default) ----
    # vif_coverage_min: trait_rg coverage below which `auto` falls back to the proxy.
    # cluster_factor_pattern / cluster_factors: how the .rds route splits factors vs panel traits.
    # rds / rds_trait_meta: optional GenomicSEM .rds input route + its trait_id->trait_category map.
    vif_coverage_min = 0.5, cluster_factor_pattern = "^C[0-9]", cluster_factors = NULL,
    rds = NULL, rds_trait_meta = NULL,
    # ---- Phase 3 (additive) ----
    # z_vector: the h2-reliability thresholds swept for the sensitivity TSVs. The primary z
    # (h2_z_threshold) is always folded in, so the primary slice always exists.
    z_vector = c(3, 4, 5),
    # ---- Phase 6 (additive uncertainty quantification) ----
    # emit_uncertainty: master toggle. TRUE appends the AUC CI (auc_abs_se/ci_lo/ci_hi) and the
    #   shape-support columns (shape_confidence, anchor_focus_ci_lo/hi, shape_posterior,
    #   shape_jackknife_stable); FALSE reproduces the Phase-5 contracts byte-for-byte.
    # auc_ci_method: delong (default, nonparametric) | hanley (closed form). auc_ci_level: CI level.
    # shape_confidence_B: MC draws for the shape support. shape_confidence_min_sig_frac: min fraction
    #   of draws with >=1 significant domain required to report a focus CI (else NA).
    emit_uncertainty = TRUE, auc_ci_method = "delong", auc_ci_level = 0.95,
    shape_confidence_B = 2000, shape_confidence_min_sig_frac = 0.05)
}

# Read a YAML config and overlay it on the defaults. `levels` is flattened to a character vector.
load_config <- function(path) {
  cfg <- yaml::read_yaml(path)
  dflt <- default_config()
  for (k in names(dflt)) if (is.null(cfg[[k]])) cfg[[k]] <- dflt[[k]]
  cfg[["levels"]]   <- as.character(unlist(cfg[["levels"]]))
  cfg[["z_vector"]] <- as.numeric(unlist(cfg[["z_vector"]]))
  validate_config(cfg)
  cfg
}

# Fail-early source-independent config validation: scientific software should reject a malformed
# run, not silently produce garbage or no scores. Input-source checks run after CLI overrides.
# Returns `cfg` invisibly on success; otherwise stop()s with an actionable message.
validate_config <- function(cfg) {
  bad <- character(0)
  add <- function(...) bad[[length(bad) + 1L]] <<- sprintf(...)

  in_unit  <- function(x) is.numeric(x) && length(x) == 1L && is.finite(x) && x >= 0 && x <= 1
  pos_int  <- function(x) is.numeric(x) && length(x) == 1L && is.finite(x) && x >= 1 && x == round(x)
  pos_num  <- function(x) is.numeric(x) && length(x) == 1L && is.finite(x) && x > 0

  if (!cfg[["vif_correlation"]] %in% c("auto", "trait_rg", "cluster_profile"))
    add("vif_correlation must be one of {auto, trait_rg, cluster_profile}, got '%s'", cfg[["vif_correlation"]])
  if (!cfg[["rank_variable"]] %in% c("abs_z", "abs_rg"))
    add("rank_variable must be one of {abs_z, abs_rg}, got '%s'", cfg[["rank_variable"]])

  if (!length(cfg[["levels"]]) || any(!nzchar(cfg[["levels"]])))
    add("levels must be a non-empty list of level names")
  else if (!cfg[["primary_level"]] %in% cfg[["levels"]])
    add("primary_level '%s' is not in levels {%s}", cfg[["primary_level"]], paste(cfg[["levels"]], collapse = ", "))

  if (!pos_int(cfg[["permutation_K"]]))  add("permutation_K must be an integer >= 1, got %s", cfg[["permutation_K"]])
  if (!pos_int(cfg[["min_category_n"]])) add("min_category_n must be an integer >= 1, got %s", cfg[["min_category_n"]])
  if (!pos_num(cfg[["h2_z_threshold"]])) add("h2_z_threshold must be > 0, got %s", cfg[["h2_z_threshold"]])
  if (!in_unit(cfg[["label_auc_min"]]))  add("label_auc_min must be in [0,1], got %s", cfg[["label_auc_min"]])
  if (!in_unit(cfg[["label_q_max"]]))    add("label_q_max must be in [0,1], got %s", cfg[["label_q_max"]])
  if (!in_unit(cfg[["hit_abs_rg"]]))     add("hit_abs_rg must be in [0,1], got %s", cfg[["hit_abs_rg"]])
  if (!in_unit(cfg[["vif_coverage_min"]])) add("vif_coverage_min must be in [0,1], got %s", cfg[["vif_coverage_min"]])
  if (!(is.numeric(cfg[["vif_min_rho"]]) && length(cfg[["vif_min_rho"]]) == 1L &&
        is.finite(cfg[["vif_min_rho"]]) && cfg[["vif_min_rho"]] >= 0))
    add("vif_min_rho must be a finite value >= 0, got %s", cfg[["vif_min_rho"]])

  if (!length(cfg[["z_vector"]]) || any(!is.finite(cfg[["z_vector"]])) || any(cfg[["z_vector"]] <= 0))
    add("z_vector must be non-empty and all entries finite and > 0")

  # ---- Phase 6 uncertainty block ----
  if (!is.logical(cfg[["emit_uncertainty"]]) || length(cfg[["emit_uncertainty"]]) != 1L)
    add("emit_uncertainty must be a single TRUE/FALSE, got %s", cfg[["emit_uncertainty"]])
  if (!cfg[["auc_ci_method"]] %in% c("delong", "hanley"))
    add("auc_ci_method must be one of {delong, hanley}, got '%s'", cfg[["auc_ci_method"]])
  if (!in_unit(cfg[["auc_ci_level"]]) || cfg[["auc_ci_level"]] <= 0 || cfg[["auc_ci_level"]] >= 1)
    add("auc_ci_level must be in (0,1), got %s", cfg[["auc_ci_level"]])
  if (!pos_int(cfg[["shape_confidence_B"]]))
    add("shape_confidence_B must be an integer >= 1, got %s", cfg[["shape_confidence_B"]])
  if (!in_unit(cfg[["shape_confidence_min_sig_frac"]]))
    add("shape_confidence_min_sig_frac must be in [0,1], got %s", cfg[["shape_confidence_min_sig_frac"]])

  if (length(bad))
    stop("Invalid config:\n  - ", paste(bad, collapse = "\n  - "), call. = FALSE)
  invisible(cfg)
}

# Source-dependent validation runs only after CLI overrides have been applied. `rds_active` includes
# both cfg$rds and the `--rds`/`rds=` override resolved by run_anchormap().
validate_config_sources <- function(cfg, rds_active = FALSE) {
  if (identical(cfg[["vif_correlation"]], "trait_rg") &&
      is.null(cfg[["trait_rg_matrix"]]) && !isTRUE(rds_active))
    stop("Invalid config:\n  - vif_correlation: trait_rg requires a trait_rg_matrix or an input ",
         "supplied via --trait-rg / --rds", call. = FALSE)
  invisible(cfg)
}

# Resolve a --config argument: a real file path is used as-is; otherwise a *bare name* is looked up
# among the configs shipped with the installed package (so `--config example_anthro` / `synthetic_rds`
# work out of the box). Returns the path unchanged when nothing matches, so the caller errors clearly.
resolve_config_path <- function(path) {
  if (is.null(path) || !nzchar(path)) return(path)
  if (file.exists(path)) return(path)
  for (cand in c(file.path("configs", path),
                 file.path("configs", paste0(path, ".yaml")),
                 file.path("configs", paste0(path, ".yml")))) {
    p <- system.file(cand, package = "anchormap")
    if (nzchar(p) && file.exists(p)) return(p)
  }
  path
}

# stage_root = config dir, or the directory above it when the config lives in a `configs/` subdir.
stage_root_of <- function(config_path) {
  base <- dirname(normalizePath(config_path, mustWork = FALSE))
  if (basename(base) == "configs") dirname(base) else base
}

resolve_path <- function(base_dir, path) {
  if (startsWith(path, "/")) path else file.path(base_dir, path)
}

# Resolve a CLI-supplied path relative to the current working directory (absolute paths pass
# through). Used for --out-dir / --rg-long / --ontology overrides so they behave like shell paths,
# independent of where the config lives.
.abs_cwd <- function(path) {
  if (is.null(path)) return(NULL)
  if (startsWith(path, "/")) path else file.path(getwd(), path)
}

# ---- rg long-table (Input A) -----------------------------------------------
# Read as character then coerce, mirroring load_long (py L67-74). Missing booleans -> FALSE
# (pandas `.get(...,"")` semantics). Precomputed z/abs_rg columns are ignored; the engine recomputes.
.LONG_REQUIRED <- c("cluster_label","trait_id","trait_category","trait_group",
                    "rg","rg_se","p","h2_trait","h2_trait_se",
                    "ldsc_converged","negative_h2","status")

read_long <- function(path) {
  df <- data.table::fread(path, sep = "\t", colClasses = "character",
                          na.strings = c("", "NA", "NaN", "NULL", "None"),
                          showProgress = FALSE)
  data.table::setDF(df)
  miss <- setdiff(.LONG_REQUIRED, names(df))
  if (length(miss)) stop(sprintf("rg long-table %s missing required column(s): %s",
                                 path, paste(miss, collapse = ", ")))
  for (col in c("rg","rg_se","p","h2_trait","h2_trait_se"))
    df[[col]] <- suppressWarnings(as.numeric(df[[col]]))
  cv <- toupper(as.character(df[["ldsc_converged"]])); df[["ldsc_converged"]] <- !is.na(cv) & cv == "TRUE"
  nv <- toupper(as.character(df[["negative_h2"]]));    df[["negative_h2"]]    <- !is.na(nv) & nv == "TRUE"
  df[["status"]] <- ifelse(is.na(df[["status"]]), "", as.character(df[["status"]]))

  # Fail early on duplicate (cluster_label, trait_id) rows: a duplicate silently changes a cluster's
  # N and corrupts the rank-based AUC, so reject rather than score garbage.
  dup_key <- paste(df[["cluster_label"]], df[["trait_id"]], sep = "\r")
  if (anyDuplicated(dup_key)) {
    ex <- unique(dup_key[duplicated(dup_key)])
    stop(sprintf("rg long-table %s has %d duplicate (cluster_label, trait_id) row(s), e.g. %s",
                 path, length(ex),
                 paste(sub("\r", "/", utils::head(ex, 5), fixed = TRUE), collapse = ", ")),
         call. = FALSE)
  }
  df
}

# ---- ontology (Input D) ----------------------------------------------------
read_ontology <- function(path, key) {
  ont <- data.table::fread(path, sep = "\t", colClasses = "character",
                           na.strings = c("", "NA"), showProgress = FALSE)
  data.table::setDF(ont)
  if (!(key %in% names(ont)))
    stop(sprintf("ontology %s lacks join column '%s'", path, key))
  # Fail early on a non-unique join key: a duplicated key fans out the merge in attach_ontology
  # (one gated trait -> many rows), silently expanding the universe.
  if (anyDuplicated(ont[[key]])) {
    ex <- unique(ont[[key]][duplicated(ont[[key]])])
    stop(sprintf("ontology %s has %d duplicate '%s' key(s), e.g. %s",
                 path, length(ex), key, paste(utils::head(ex, 5), collapse = ", ")),
         call. = FALSE)
  }
  ont
}
