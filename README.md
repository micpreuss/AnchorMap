# AnchorMap

**AnchorMap** scores — competitively, size-aware, and correlation-aware — **which ontology domain
each latent cluster factor anchors to**, how confidently, and whether the anchor is *sharp* or
*diffuse*. Given cluster factors and their genetic correlations (`rg`) to a trait panel, it ranks
each cluster's affinity to every ontology category, calibrates significance against a permutation
null and a redundancy-deflated test, pools a signed effect size, and emits an auto-label per cluster
plus publication-ready figures.

It ships as an installable **R package** (`anchormap`) with a library API and two command-line
entry points, and as a pinned, self-validating **Docker image**.

## What it does

For each `(cluster, ontology level, category)` it computes:

- a competitive **Mann–Whitney AUC** (in-category vs out-of-category rank enrichment) with a
  label-permutation **p-value**;
- a redundancy-aware significance via **Li & Ji** effective number of tests and a **CAMERA**-style
  variance-inflation factor (so correlated traits don't inflate confidence);
- an inverse-variance **Fisher-z pooled `rg`** (signed, with a CI) and a **coherence** score;
- a Fisher **over-representation** test;
- **BH-FDR** `q`, a per-cluster **auto-label**, and an **anchor shape** (weak / sharp / diffuse / focal).

A parallel **sensitivity sweep** re-runs the whole pipeline across reliability thresholds and flags
whether each cluster's label is stable.

## Install

**R package** (engine + figures):

```r
# install.packages("remotes")
remotes::install_github("micpreuss/AnchorMap")
```

Requires R ≥ 4.4. Dependencies (`data.table`, `yaml`, `future`, `future.apply`, `ggplot2`,
`patchwork`, `scales`, `ggrepel`) install automatically; `poolr` and `ragg` are optional.

**Docker image** (no R setup; pinned, reproducible):

```bash
docker build -t anchormap:0.1.0 -f docker/Dockerfile .
```

## Quick start

A self-contained synthetic example ships with the package — no external data needed:

```r
library(anchormap)
run_anchormap("synthetic_rds", out_dir = "results/demo", threads = 2)
#> ... C5_sub0 -> anthro [sharp] ...
run_plots("synthetic_rds_plots", out_dir = "results/demo/figures")
```

Same thing on the command line:

```bash
anchor_map  --config synthetic_rds --out-dir results/demo --threads 2
plot_anchors --config synthetic_rds_plots --out-dir results/demo/figures
```

(If you installed only the package, the shims are `Rscript -e 'anchormap:::cli_anchor_map()' ...`; the
Docker image puts `anchor_map` / `plot_anchors` on the `PATH`.)

## Running on your own data

Copy a shipped template and edit the paths, or pass everything on the CLI:

```bash
# config-driven (recommended): copy example_anthro.yaml / example_disease.yaml and edit
anchor_map --config myrun.yaml --out-dir results/run1 --threads 4

# or point at files directly, no YAML editing:
anchor_map --config example_anthro.yaml \
  --rg-long data/cluster_trait_rg.tsv --ontology data/ontology.tsv \
  --out-dir results/run1 --threads 4
```

### Via Docker (mount your data and an output dir)

```bash
docker run --rm \
  -v "$PWD/data:/data:ro" -v "$PWD/out:/out" \
  anchormap:0.1.0 \
  anchor_map --config /data/myrun.yaml --out-dir /out --threads 4
```

### Command-line help

Both entry points print full help with `--help`:

```bash
anchor_map --help        # all flags, defaults, and an example
plot_anchors --help
```

Key flags: `--config` (a YAML path **or** a bare shipped-config name); input overrides `--rds`,
`--rg-long`, `--trait-rg`, `--ontology`; output control `--out-dir`, `--run-label`; engine
`--threads`, `--z-vector`.

## Inputs

Two interchangeable routes converge on the same engine:

- **Route A — two TSVs:**
  - a **cluster × trait `rg` long-table** with columns `cluster_label, trait_id, trait_category,
    trait_group, rg, rg_se, p, h2_trait, h2_trait_se, ldsc_converged, negative_h2, status`
    (one row per cluster × trait);
  - optionally a **trait × trait `rg` matrix** (an LDSC `--rg` summary: `p1, p2, rg, …, CONVERGED`)
    as the within-category redundancy source.
- **Route C — a GenomicSEM `ldsc()` `.rds`** (`--rds`): AnchorMap standardizes `$S` → `rg`, derives
  `rg_se` by the delta method on `$V`, splits cluster factors from panel traits, and derives both the
  long-table and the trait × trait block — so the rest of the engine is route-agnostic.

Both routes also need an **ontology TSV** mapping traits/categories to domains, with an
`anchor_eligible` flag (ineligible categories can be scored but never become an auto-label).

## Outputs

Written to `--out-dir`:

| File | Contents |
| --- | --- |
| `category_anchor_scores.tsv` | one row per `(cluster, level, category)`: `n, n_eff, n_hit, rho_bar, vif, auc_abs, auc_signed, perm_p, vif_z, vif_p, pooled_rg [ci], coherence, odds_ratio, fisher_p, q, rank` |
| `cluster_anchor_labels.tsv` | one row per cluster: `auto_label, anchor_shape, anchor_margin, anchor_focus, n_sig_domains, top_*`, profile |
| `sensitivity_z_scores.tsv` / `sensitivity_z_labels.tsv` | the two tables stacked across the reliability-threshold sweep (`z_threshold` column; labels carry a per-cluster `label_stable` flag) |
| `anchormap.log` | timestamped steps ending in a `FINISHED` line |
| `figures/` (from `plot_anchors`) | lollipop small-multiples, a cluster × category dot-heatmap, an AUC-vs-coherence diagnostic, and a cross-cluster specificity heatmap + diagonal (PNG + PDF), plus `cluster_distinctive_categories.tsv` |

## Library API

```r
run_anchormap(config_path, threads, rds, z_vector, out_dir, run_label, rg_long, trait_rg, ontology)
run_sensitivity(df, ont, cfg, sroot, z_vector, threads, trait_rg_override)
score_at_z(df, ont, cfg, sroot, z, trait_rg_override)
read_ldsc_rds(path); read_rds_route(path, cfg, sroot)
run_plots(config_path, q_sig, rg_floor, min_clusters, out_dir)
```

See `?run_anchormap` for the full reference.

## Reproducibility

The Docker image pins R (`rocker/r-ver:4.6.0`) and every package to a single dated CRAN snapshot, and
**self-validates at build time** — the build fails if the engine or the figure stack regresses.
Reference the image by version tag, never `latest`.

## Validation

AnchorMap's deterministic outputs are checked against an independent reference implementation and an
analytic test suite (`testthat`). See [VALIDATION.md](VALIDATION.md).

## License

MIT — see [LICENSE.md](LICENSE.md).
