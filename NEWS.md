# anchormap 0.1.2

- Made ontology input authoritative when same-named columns are also present in the correlation table,
  preventing stale domain or eligibility annotations from being retained.
- Hardened `--z-vector` validation to reject empty, non-finite, and non-positive values.
- Included both PNG and PDF outputs in the plotting return value and figure count.
- Corrected the R package MIT license layout and clarified the Docker image reproducibility scope.

# anchormap 0.1.1

- Added push/pull-request CI for the test suite and end-to-end engine/figure smoke tests, plus a
  separate `R CMD check` job that fails on warnings or errors.
- Added strict configuration, input-uniqueness, and ontology validation with actionable failures.
- Fixed source validation sequencing so `--trait-rg` and `--rds` overrides satisfy explicit
  `trait_rg` redundancy mode before source-dependent checks run.
- Added the single-track `--in-dir` plotting convenience, an explicit multi-track guard, and clear
  score/label file, schema, and requested-level validation.
- Declared minimum supported `future` and `future.apply` versions.
- Improved public-project metadata, contribution guidance, citation metadata, and MIT license
  detection on GitHub.
