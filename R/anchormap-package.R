#' @keywords internal
"_PACKAGE"

## Package-level imports. The modules call data.table and ggplot2 idioms unqualified
## (the figure builders are a ggplot2 DSL); everything else is reached via `::`.
#' @import data.table
#' @import ggplot2
#' @importFrom stats aggregate pnorm cor fisher.test cov2cor hclust dist median setNames
#' @importFrom utils head
#' @importFrom grDevices png cairo_pdf
NULL

## Quiet R CMD check "no visible binding" notes for the column-name symbols referenced
## inside ggplot2 aes() and data.table reshapes (they are columns, not globals).
utils::globalVariables(c(
  "auc_abs", "pooled_rg", "coherence", "category", "cluster_label", "rank", "q",
  "y", "yend", "a", "n", "sz", "sig", "z", "zmask", "lab", "auto_label", "isauto",
  "ax", "auc_abs_ci_lo", "auc_abs_ci_hi"
))
