# Fail-early input validation: malformed config / data must stop() before scoring; valid shipped
# inputs must still pass unchanged (guards the Phase-1 oracle parity).

# ---- config validation ------------------------------------------------------

# Mirror load_config's coercion so validate_config sees the same shapes it does in production.
default_config_levels_chr <- function() {
  cfg <- default_config()
  cfg[["levels"]]   <- as.character(unlist(cfg[["levels"]]))
  cfg[["z_vector"]] <- as.numeric(unlist(cfg[["z_vector"]]))
  cfg
}

test_that("default config is valid", {
  expect_silent(validate_config(default_config_levels_chr()))
})

test_that("validate_config rejects bad enumerations and ranges", {
  base <- default_config_levels_chr()
  expect_error(validate_config(modifyList(base, list(vif_correlation = "nope"))), "vif_correlation")
  expect_error(validate_config(modifyList(base, list(rank_variable = "rg"))),     "rank_variable")
  expect_error(validate_config(modifyList(base, list(label_q_max = 2))),          "label_q_max")
  expect_error(validate_config(modifyList(base, list(label_auc_min = -0.1))),     "label_auc_min")
  expect_error(validate_config(modifyList(base, list(permutation_K = 0))),        "permutation_K")
  expect_error(validate_config(modifyList(base, list(min_category_n = 0))),       "min_category_n")
  expect_error(validate_config(modifyList(base, list(h2_z_threshold = 0))),       "h2_z_threshold")
  expect_error(validate_config(modifyList(base, list(primary_level = "ghost"))),  "primary_level")
  expect_error(validate_config(modifyList(base, list(z_vector = c(3, -1)))),      "z_vector")
})

test_that("validate_config_sources checks trait_rg after overrides", {
  base <- default_config_levels_chr()
  base[["vif_correlation"]] <- "trait_rg"
  expect_silent(validate_config(base))
  expect_error(validate_config_sources(base), "--trait-rg / --rds")
  expect_silent(validate_config_sources(modifyList(base, list(trait_rg_matrix = "x.tsv"))))
  expect_silent(validate_config_sources(base, rds_active = TRUE))
})

test_that("run_anchormap accepts CLI trait_rg and rds source overrides", {
  # TSV route: explicit trait_rg mode is satisfied by trait_rg= even though the YAML has no matrix.
  pkg_file <- function(...) system.file(..., package = "anchormap", mustWork = TRUE)
  cfg_tsv <- yaml::read_yaml(pkg_file("configs", "example_disease.yaml"))
  cfg_tsv[["rg_long"]] <- pkg_file("fixtures", "example_rg_long.tsv")
  cfg_tsv[["ontology"]] <- pkg_file("fixtures", "example_ontology.tsv")
  cfg_tsv[["vif_correlation"]] <- "trait_rg"
  cfg_tsv[["trait_rg_matrix"]] <- NULL
  cfg_tsv_path <- tempfile(fileext = ".yaml"); yaml::write_yaml(cfg_tsv, cfg_tsv_path)

  tids <- unique(data.table::fread(cfg_tsv[["rg_long"]])[["trait_id"]])
  edges <- expand.grid(p1 = tids, p2 = tids, stringsAsFactors = FALSE)
  edges[["rg"]] <- ifelse(edges[["p1"]] == edges[["p2"]], 1, 0.1)
  edges[["CONVERGED"]] <- "TRUE"
  trait_rg_path <- tempfile(fileext = ".tsv")
  data.table::fwrite(edges, trait_rg_path, sep = "\t")
  expect_no_error(run_anchormap(cfg_tsv_path, trait_rg = trait_rg_path,
                                out_dir = tempfile(), z_vector = 4))

  # RDS route: rds= similarly satisfies explicit trait_rg mode after config loading.
  cfg_rds <- yaml::read_yaml(pkg_file("configs", "synthetic_rds.yaml"))
  cfg_rds[["rds"]] <- NULL
  cfg_rds[["ontology"]] <- pkg_file("fixtures", "synthetic_panel_ontology.tsv")
  cfg_rds[["vif_correlation"]] <- "trait_rg"
  cfg_rds_path <- tempfile(fileext = ".yaml"); yaml::write_yaml(cfg_rds, cfg_rds_path)
  expect_no_error(run_anchormap(
    cfg_rds_path, rds = pkg_file("fixtures", "synthetic_ldsc_panel.rds"),
    out_dir = tempfile(), z_vector = 4))
})

test_that("shipped configs load and validate", {
  for (nm in c("example_disease", "example_anthro", "synthetic_rds")) {
    p <- system.file(file.path("configs", paste0(nm, ".yaml")), package = "anchormap")
    skip_if(!nzchar(p))
    expect_silent(load_config(p))
  }
})

# ---- data uniqueness --------------------------------------------------------

.long_header <- function() paste(
  "cluster_label", "trait_id", "trait_category", "trait_group", "rg", "rg_se", "p",
  "h2_trait", "h2_trait_se", "ldsc_converged", "negative_h2", "status", sep = "\t")

.long_row <- function(cl, tid, cat = "A") paste(
  cl, tid, cat, "disease", "0.3", "0.05", "0.01", "0.2", "0.02", "TRUE", "FALSE", "success", sep = "\t")

test_that("read_long rejects duplicate (cluster_label, trait_id) rows", {
  tmp <- tempfile(fileext = ".tsv")
  writeLines(c(.long_header(), .long_row("C0", "T1"), .long_row("C0", "T1")), tmp)
  expect_error(read_long(tmp), "duplicate")
  # distinct rows are fine
  writeLines(c(.long_header(), .long_row("C0", "T1"), .long_row("C0", "T2")), tmp)
  expect_silent(read_long(tmp))
})

test_that("read_ontology rejects a non-unique join key", {
  tmp <- tempfile(fileext = ".tsv")
  writeLines(c("trait_id\tanthro_class\tanchor_eligible",
               "T1\tsize\tTRUE", "T1\tlength\tTRUE"), tmp)
  expect_error(read_ontology(tmp, "trait_id"), "duplicate")
  writeLines(c("trait_id\tanthro_class\tanchor_eligible",
               "T1\tsize\tTRUE", "T2\tlength\tTRUE"), tmp)
  expect_silent(read_ontology(tmp, "trait_id"))
})

# ---- ontology levels --------------------------------------------------------

test_that("attach_ontology errors on a configured level absent from the ontology", {
  g <- data.frame(trait_id = c("T1", "T2"), cluster_label = "C0", stringsAsFactors = FALSE)
  ont <- data.frame(trait_id = c("T1", "T2"), anthro_class = c("a", "b"),
                    anchor_eligible = "TRUE", stringsAsFactors = FALSE)
  expect_error(attach_ontology(g, ont, "trait_id", c("anthro_class", "ghost_level")), "ghost_level")
  expect_silent(attach_ontology(g, ont, "trait_id", c("anthro_class")))
})

test_that("attach_ontology errors when nothing matches the ontology", {
  g <- data.frame(trait_id = c("X1", "X2"), cluster_label = "C0", stringsAsFactors = FALSE)
  ont <- data.frame(trait_id = c("T1", "T2"), anthro_class = c("a", "b"),
                    anchor_eligible = "TRUE", stringsAsFactors = FALSE)
  expect_error(attach_ontology(g, ont, "trait_id", c("anthro_class")), "matched the ontology")
})

test_that("attach_ontology: the ontology overrides stale columns carried in the rg table", {
  # An adversarial rg table carries its own (stale) domain columns that contradict the curated
  # ontology; the ontology must win and the override must be announced, not silent.
  g <- data.frame(trait_id = c("T1", "T2"), cluster_label = "C0",
                  anthro_class = c("WRONG", "WRONG"), anchor_eligible = c("FALSE", "FALSE"),
                  stringsAsFactors = FALSE)
  ont <- data.frame(trait_id = c("T1", "T2"), anthro_class = c("size", "length"),
                    anchor_eligible = c("TRUE", "TRUE"), stringsAsFactors = FALSE)
  expect_message(out <- attach_ontology(g, ont, "trait_id", c("anthro_class")),
                 "ontology overrides")
  expect_setequal(out[["anthro_class"]], c("size", "length"))   # ontology values, not "WRONG"
  expect_true(all(out[["anchor_eligible"]]))                    # ontology TRUE, not the stale FALSE
  expect_false(any(out[["anthro_class"]] == "WRONG"))
})

# ---- plot input validation --------------------------------------------------

.write_plot_fixture <- function(dir) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  scores <- data.frame(
    level = "anthro_class", eligible = "True", category = "anthro", cluster_label = "C0",
    pooled_rg = 0.4, q = 0.01, coherence = 1, auc_abs = 0.9, n = 3, rank = 1)
  labels <- data.frame(cluster_label = "C0", auto_label = "anthro", anchor_shape = "sharp")
  data.table::fwrite(scores, file.path(dir, "category_anchor_scores.tsv"), sep = "\t")
  data.table::fwrite(labels, file.path(dir, "cluster_anchor_labels.tsv"), sep = "\t")
  invisible(dir)
}

.plot_track <- function(level = "anthro_class", name = "synthetic")
  list(name = name, level = level, scores = "category_anchor_scores.tsv",
       labels = "cluster_anchor_labels.tsv")

test_that("load_track accepts a single-track in_dir containing spaces", {
  d <- .write_plot_fixture(file.path(tempdir(), "AnchorMap plot input with spaces"))
  tr <- load_track(.plot_track(), stage_root = tempdir(), in_dir = d)
  expect_equal(nrow(tr[["s"]]), 1L)
  expect_identical(tr[["s"]][["category"]], "anthro")
})

test_that("load_track reports missing files and empty requested levels clearly", {
  missing_dir <- tempfile("missing plot inputs ")
  expect_error(load_track(.plot_track(), tempdir(), missing_dir), "file not found")

  d <- .write_plot_fixture(file.path(tempdir(), "AnchorMap plot level fixture"))
  expect_error(load_track(.plot_track(level = "ghost"), tempdir(), d),
               "no eligible rows at level 'ghost'.*levels present: anthro_class")
})

test_that("run_plots rejects global in_dir for multi-track configs", {
  d <- .write_plot_fixture(file.path(tempdir(), "AnchorMap multi track inputs"))
  cfg <- list(out_dir = tempfile(), tracks = list(.plot_track(name = "disease"),
                                                   .plot_track(name = "anthro")))
  cfg_path <- tempfile(fileext = ".yaml"); yaml::write_yaml(cfg, cfg_path)
  expect_error(run_plots(cfg_path, in_dir = d),
               "--in-dir is a single-track convenience.*disease, anthro")
})
