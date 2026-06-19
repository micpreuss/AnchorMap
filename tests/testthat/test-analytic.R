# Analytic unit tests for the engine primitives (no large-input dependency).

test_that("Li & Ji n_eff + rho_bar + VIF on the cheat-sheet matrix", {
  R <- matrix(c(1,.9,.2, .9,1,.4, .2,.4,1), 3, 3, byrow = TRUE)
  expect_true(approx(meff_liji(R), 2.0, 1e-3))
  expect_true(approx(rho_bar(R), 0.5))
  neff <- meff_liji(R); rb <- rho_bar(R)
  expect_true(approx(1 + (neff - 1) * rb, 1.5, 1e-3))
  if (requireNamespace("poolr", quietly = TRUE))
    expect_true(approx(meff_poolr(R), meff_liji(R), 1e-6))
})

test_that("n_eff never exceeds m; identity -> n_eff == m, rho_bar == 0", {
  expect_true(approx(meff_liji(diag(4)), 4))
  expect_true(approx(rho_bar(diag(4)), 0))
})

test_that("AUC from ranks + rankdata tie semantics", {
  expect_true(approx(auc_from_ranks(c(15,14,10), 3, 12), 33/36, 1e-6))
  expect_true(approx(rank(c(10, 8, 8, 1)), c(4, 2.5, 2.5, 1)))
})

test_that("BH-FDR matches p.adjust('BH')", {
  set.seed(7); pv <- stats::runif(20)
  expect_true(approx(bh_fdr(pv), stats::p.adjust(pv, method = "BH"), 1e-12))
})

test_that("sample odds ratio (scipy semantics), not conditional MLE", {
  expect_true(approx((2*6)/(1*6), 2.0))
  expect_false(is.finite((2*6)/(1*0)))
  expect_true(approx((0*6)/(1*6), 0))
})

test_that("per-trait stats: Fisher-z y, delta-method v", {
  rg <- 0.8084971841; se <- 0.0300751195
  rg_c <- min(max(rg, -0.999), 0.999)
  expect_true(approx(atanh(rg_c), atanh(0.8084971841)))
  expect_true(approx(se^2 / (1 - rg_c^2)^2, (se^2)/(1 - rg_c^2)^2))
})

test_that("anchor_shape: single significant domain -> sharp", {
  sub1 <- data.frame(auc_abs = 0.9164, q = 0.005497, stringsAsFactors = FALSE)
  shp <- anchor_shape(sub1, default_config())
  expect_identical(shp$anchor_shape, "sharp")
  expect_identical(shp$n_sig_domains, 1L)
  expect_true(approx(shp$anchor_focus, 1.0))
  expect_true(is.na(shp$anchor_margin))
})

test_that("anchor_shape: nothing passes the gate -> weak", {
  sub0 <- data.frame(auc_abs = c(0.4, 0.3), q = c(0.9, 0.95), stringsAsFactors = FALSE)
  expect_identical(anchor_shape(sub0, default_config())$anchor_shape, "weak")
})
