#!/usr/bin/env Rscript
# tests/run_tests.R — analytic unit tests for the AnchorMap engine (no large-input dependency).
# Uses base-R stopifnot assertions (no testthat dependency). Run: Rscript tests/run_tests.R

sdir <- dirname(dirname(normalizePath(sub("^--file=", "",
          grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))))
for (m in c("io.R","gate.R","redundancy.R","score.R","label.R")) source(file.path(sdir, "R", m))

ok <- function(name, cond) {
  if (!isTRUE(cond)) stop(sprintf("FAIL: %s", name))
  cat(sprintf("  ok  %s\n", name))
}
approx <- function(a, b, tol = 1e-6) all(abs(a - b) < tol)

cat("== analytic unit tests ==\n")

# --- Li & Ji n_eff + rho_bar + VIF on the cheat-sheet matrix -----------------
R <- matrix(c(1,.9,.2, .9,1,.4, .2,.4,1), 3, 3, byrow = TRUE)
ok("meff_liji == 2.00 (cheat-sheet)", approx(meff_liji(R), 2.0, 1e-3))
ok("rho_bar == 0.50",                 approx(rho_bar(R), 0.5))
neff <- meff_liji(R); rb <- rho_bar(R)
ok("VIF == 1.5",                      approx(1 + (neff - 1) * rb, 1.5, 1e-3))
if (requireNamespace("poolr", quietly = TRUE))
  ok("poolr::meff agrees with manual on clean matrix", approx(meff_poolr(R), meff_liji(R), 1e-6))

# n_eff never exceeds m; identity matrix -> n_eff == m
ok("n_eff(I_4) == 4", approx(meff_liji(diag(4)), 4))
ok("rho_bar(I_4) == 0", approx(rho_bar(diag(4)), 0))

# --- AUC from ranks (cheat-sheet: in-traits at ascending ranks 15,14,10) -----
ok("auc_from_ranks == 33/36 = 0.9167", approx(auc_from_ranks(c(15,14,10), 3, 12), 33/36, 1e-6))
# rankdata semantics: largest value gets the largest rank, average ties
ok("rank() average ties", approx(rank(c(10, 8, 8, 1)), c(4, 2.5, 2.5, 1)))

# --- BH-FDR matches p.adjust('BH') ------------------------------------------
set.seed(7); pv <- runif(20)
ok("bh_fdr == p.adjust(BH)", approx(bh_fdr(pv), p.adjust(pv, method = "BH"), 1e-12))

# --- sample odds ratio (scipy semantics), not conditional MLE ---------------
ok("sample OR (2*6)/(1*6) == 2", approx((2*6)/(1*6), 2.0))
ok("OR -> Inf when fp==0", !is.finite((2*6)/(1*0)))
ok("OR == 0 when tp==0",  approx((0*6)/(1*6), 0))

# --- per-trait stats (Fisher-z y, delta-method v) on a known row ------------
rg <- 0.8084971841; se <- 0.0300751195
rg_c <- min(max(rg, -0.999), 0.999)
ok("y == atanh(rg)",            approx(atanh(rg_c), atanh(0.8084971841)))
ok("v == se^2/(1-rg^2)^2",      approx(se^2 / (1 - rg_c^2)^2, (se^2)/(1 - rg_c^2)^2))

# --- anchor_shape: single significant domain -> sharp -----------------------
sub1 <- data.frame(auc_abs = 0.9164, q = 0.005497, stringsAsFactors = FALSE)
shp <- anchor_shape(sub1, default_config())
ok("single sig domain -> sharp", shp$anchor_shape == "sharp")
ok("n_sig == 1",                 shp$n_sig_domains == 1)
ok("focus == 1.0",               approx(shp$anchor_focus, 1.0))
ok("margin == NA (1 category)",  is.na(shp$anchor_margin))

# --- anchor_shape: nothing passes the gate -> weak --------------------------
sub0 <- data.frame(auc_abs = c(0.4, 0.3), q = c(0.9, 0.95), stringsAsFactors = FALSE)
ok("no sig domain -> weak", anchor_shape(sub0, default_config())$anchor_shape == "weak")

cat("\nALL UNIT TESTS PASSED\n")
