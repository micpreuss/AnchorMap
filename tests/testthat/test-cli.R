# CLI argument parsing (Part B): help, input/output overrides, unknown-flag errors.

test_that("engine parser captures config, output, and input overrides", {
  a <- parse_engine_args(c("--config", "x.yaml", "--out-dir", "out", "--run-label", "r1",
                           "--threads", "4", "--rds", "a.rds", "--rg-long", "b.tsv",
                           "--trait-rg", "c.tsv", "--ontology", "o.tsv", "--z-vector", "2,3,4"))
  expect_identical(a$config, "x.yaml")
  expect_identical(a$out_dir, "out")
  expect_identical(a$run_label, "r1")
  expect_identical(a$threads, 4L)
  expect_identical(a$rds, "a.rds")
  expect_identical(a$rg_long, "b.tsv")
  expect_identical(a$trait_rg, "c.tsv")
  expect_identical(a$ontology, "o.tsv")
  expect_equal(a$z_vector, c(2, 3, 4))
  expect_false(a$help)
})

test_that("engine parser handles --help / -h and unknown flags", {
  expect_true(parse_engine_args("--help")$help)
  expect_true(parse_engine_args("-h")$help)
  expect_error(parse_engine_args(c("--bogus", "1")), "unknown option")
  expect_error(parse_engine_args("--config"), "needs a value")
})

test_that("plots parser captures overrides + help + unknown flags", {
  a <- parse_plots_args(c("--config", "p.yaml", "--out-dir", "fig", "--q-sig", "0.01",
                          "--rg-floor", "0.2", "--min-clusters", "3"))
  expect_identical(a$config, "p.yaml")
  expect_identical(a$out_dir, "fig")
  expect_equal(a$q_sig, 0.01)
  expect_equal(a$rg_floor, 0.2)
  expect_identical(a$min_clusters, 3L)
  expect_true(parse_plots_args("--help")$help)
  expect_error(parse_plots_args("--nope"), "unknown option")
})

test_that("usage text mentions the new flags", {
  u <- engine_usage()
  expect_true(grepl("--out-dir", u))
  expect_true(grepl("--rg-long", u))
  expect_true(grepl("--help", u))
})
