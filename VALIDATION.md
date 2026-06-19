# Validation

AnchorMap is validated three ways: an analytic unit suite, cross-implementation parity on the
deterministic outputs, and built-in positive/negative controls.

## Analytic unit tests (`testthat`)

`tests/testthat/` exercises the engine primitives against closed-form expectations, including:

- **Li & Ji** effective number of tests and mean pairwise correlation, and the resulting VIF, on a
  known correlation matrix (cross-checked against `poolr::meff` on positive-definite inputs);
- **Mann–Whitney AUC** from summed ranks and average-tie rank semantics;
- **BH-FDR** equal to `stats::p.adjust(method = "BH")` to machine precision;
- the **sample odds ratio** `(a·d)/(b·c)` (not the conditional-MLE estimate);
- the **delta-method `rg_se`** from a GenomicSEM `$V` submatrix, checked against a numeric finite
  difference of the standardization map;
- `.rds`-route ↔ TSV-route **round-trip** identity on every deterministic score column;
- the `vif_correlation: auto` **fallback** branches (trait_rg → cluster-profile proxy → VIF = 1);
- **VIF-invariance**: AUC, ranks, `pooled_rg`, and coherence are unchanged across redundancy sources
  (VIF affects only `vif_z` / `vif_p` / CI width);
- **sensitivity-sweep** determinism: the `z == primary` slice reproduces the single-z primary output
  bit-for-bit, the sweep is invariant to thread count, and `label_stable` matches auto-label constancy.

Run them with:

```r
testthat::test_local()
```

## Cross-implementation parity

The deterministic columns (everything except the Monte-Carlo `perm_p` / `q`) match an independent
reference implementation to machine precision, with identical auto-labels, anchor shapes, and ranks.
`perm_p` / `q` agree within Monte-Carlo error (different RNG streams), so the parity gate anchors on
the deterministic `vif_p` and on label stability rather than on `perm_p`. The cross-cluster
specificity z used by the figures reproduces the reference `cluster_distinctive_categories.tsv`
exactly.

## Controls

- **Positive control:** a cluster that genuinely anchors to a single domain is recovered as that
  domain with a *sharp* shape and a significant, redundancy-deflated test. The shipped synthetic
  example reproduces this (`C5_sub0 → anthro [sharp]`), and it is one of the two build-time self-tests
  baked into the Docker image.
- **Negative control / forbidden false-positive:** a category flagged `anchor_eligible = FALSE` may
  be *scored* but can **never** become a cluster's auto-label, and a cluster with no real enrichment
  is labelled *ambiguous* with a *weak* shape.

## Build-time self-tests

`docker build` runs two checks and fails on regression: (1) the synthetic engine run recovers the
positive control and writes a `FINISHED` log; (2) the ggplot2 / ragg / cairo figure stack renders a
PNG. A bad dependency, engine change, or figure regression therefore breaks the build.
