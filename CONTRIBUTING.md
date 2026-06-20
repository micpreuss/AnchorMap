# Contributing to AnchorMap

Thanks for your interest in improving AnchorMap. This is a small R package validated by
cross-language parity against a Python reference; the bar for changes is that the **deterministic
outputs stay reproducible**.

## Development setup

Requires R ≥ 4.4. Install the package and its dependencies from the repository root:

```bash
R CMD INSTALL .
# or, in R:  remotes::install_deps(dependencies = TRUE)
```

## Before opening a pull request

Run the same checks CI runs (`.github/workflows/ci.yml`):

1. **Test suite** — the `testthat` suite must pass:

   ```bash
   Rscript -e 'testthat::test_local(stop_on_failure = TRUE)'
   ```

2. **Engine + figures smoke run** — the synthetic example must still recover the positive control
   (`C5_sub0 -> anthro [sharp]`) and render figures:

   ```bash
   Rscript inst/scripts/anchor_map.R   --config synthetic_rds        --out-dir results/ci --threads 2
   Rscript inst/scripts/plot_anchors.R --config synthetic_rds_plots  --in-dir results/ci --out-dir results/ci/figures
   ```

3. **Package check:**

   ```bash
   R CMD build . && R CMD check anchormap_*.tar.gz
   ```

CI runs the test suite, smoke run, and package check on every push and pull request, so these must be
green. The full MIT text intentionally stays in `LICENSE` so GitHub detects it; R therefore reports
one tolerated `License stub is invalid DCF` NOTE. A future CRAN submission would restore R's two-line
MIT stub.

## Conventions

- **Determinism is load-bearing.** Anything touching the z-sweep, RNG, or scoring must preserve the
  byte-identical primary-slice parity (see `CLAUDE.md` and the tests). `perm_p`/`q` are compared
  distributionally, never bit-for-bit, across languages.
- **Config-over-CLI.** New parameters belong in the YAML config (with a default in
  `default_config()` and a check in `validate_config()`), not as new CLI flags.
- **Commits** use [Conventional Commits](https://www.conventionalcommits.org/) tags and carry no
  `Co-Authored-By` footer (match the existing history).
- Read real input headers before parsing; prefer long-format over wide; keep comments "smart and slim".

## Reporting issues

Please open an issue with a minimal reproducer (config + a small synthetic input where possible) and
the `anchormap.log` from the run.
