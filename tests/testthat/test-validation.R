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

test_that("validate_config rejects trait_rg without a matrix source", {
  base <- default_config_levels_chr()
  base[["vif_correlation"]] <- "trait_rg"
  expect_error(validate_config(base), "trait_rg_matrix")
  expect_silent(validate_config(modifyList(base, list(trait_rg_matrix = "x.tsv"))))
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
