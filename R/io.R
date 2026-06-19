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
    z_vector = c(3, 4, 5))
}

# Read a YAML config and overlay it on the defaults. `levels` is flattened to a character vector.
load_config <- function(path) {
  cfg <- yaml::read_yaml(path)
  dflt <- default_config()
  for (k in names(dflt)) if (is.null(cfg[[k]])) cfg[[k]] <- dflt[[k]]
  cfg[["levels"]]   <- as.character(unlist(cfg[["levels"]]))
  cfg[["z_vector"]] <- as.numeric(unlist(cfg[["z_vector"]]))
  cfg
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
  df
}

# ---- ontology (Input D) ----------------------------------------------------
read_ontology <- function(path, key) {
  ont <- data.table::fread(path, sep = "\t", colClasses = "character",
                           na.strings = c("", "NA"), showProgress = FALSE)
  data.table::setDF(ont)
  if (!(key %in% names(ont)))
    stop(sprintf("ontology %s lacks join column '%s'", path, key))
  ont
}
