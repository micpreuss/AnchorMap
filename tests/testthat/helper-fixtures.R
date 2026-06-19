# Shared test helpers: locate the shipped synthetic fixtures + numeric comparison utilities.

# Path to a synthetic fixture shipped under inst/fixtures (resolved whether the package is installed
# or loaded via devtools/load_all).
fx <- function(name) {
  p <- system.file("fixtures", name, package = "anchormap")
  if (!nzchar(p)) skip(sprintf("fixture not found: %s", name))
  p
}

approx <- function(a, b, tol = 1e-6) all(abs(a - b) < tol)

# Inf/NA-robust numeric equality: both-Inf-same-sign and both-NA count as equal. Used for score-row
# comparisons where odds_ratio can be Inf/0.
eqnum <- function(a, b, tol = 1e-6) {
  a <- as.numeric(a); b <- as.numeric(b)
  both_inf <- is.infinite(a) & is.infinite(b) & (sign(a) == sign(b))
  both_na  <- is.na(a) & is.na(b)
  d <- abs(a - b); d[both_inf | both_na] <- 0
  all(d <= tol)
}
