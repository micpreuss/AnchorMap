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

### A minimal example

Suppose one cluster `C0` has these genetic correlations (`rg`) to six diseases. Each disease
(`trait_id`) carries a category (`trait_category`) that rolls up into a `domain`:

```text
trait_id     rg      trait_category   →  domain
T2D          0.45    E4_DM2           →  Endocrine-metabolic
OBESITY      0.40    E4_OBESITY       →  Endocrine-metabolic
HYPERLIPID   0.38    E4_HYPERCHOL     →  Endocrine-metabolic
HTN          0.08    I9_HYPTENS       →  Circulatory
CHD         -0.05    I9_CHD           →  Circulatory
AFIB         0.06    I9_AF            →  Circulatory
```

AnchorMap groups `C0`'s diseases by `domain` and asks **which domain its correlations concentrate in**
(ranking in-domain vs out-of-domain traits — the AUC):

```text
domain                n   mean|rg|   AUC
Endocrine-metabolic   3     0.41     1.00   ← C0's rg ranks top here  →  anchor
Circulatory           3     0.06     0.00
```

`C0`'s correlations sit almost entirely on the endocrine diseases (its 3 endocrine `rg` values
out-rank every circulatory one, so AUC = 1.00), and the margin over the runner-up is large, so the
auto-label is:

```text
C0  →  Endocrine-metabolic  [sharp]
```

That is the anchor: the domain a cluster is most enriched for, with an **anchor shape**
(`sharp` here — a clean, dominant winner) and an FDR `q`-value for confidence. Note the roles —
`trait_id` identifies each disease (and drives the redundancy correction); `trait_category → domain`
is what the cluster is mapped onto. The shipped example pair
([Running on your own data](#running-on-your-own-data)) is this same logic over 3 clusters × 9 diseases —
run it to reproduce `C0 → Endocrine-metabolic [sharp]`.

## Install

**R package** (engine + figures):

```r
# install.packages("remotes")
remotes::install_github("micpreuss/AnchorMap")
```

Requires R ≥ 4.4. Dependencies (`data.table`, `yaml`, `optparse`, `future`, `future.apply`, `ggplot2`,
`patchwork`, `scales`, `ggrepel`) install automatically; `poolr` and `ragg` are optional.

**Docker image**
The only two things you install are [git](https://git-scm.com/downloads) and
[Docker Desktop](https://www.docker.com/products/docker-desktop/) — **R and every package live inside
the image**, so there is nothing else to manage. Download the source and build the image once:

```bash
git clone https://github.com/micpreuss/AnchorMap.git    # download the source
cd AnchorMap                                             # enter the project folder
docker build -t anchormap:0.1.0 -f docker/Dockerfile .  # build the image (a few minutes)
```

```bash
docker run --rm -v "$PWD:/work" -w /work anchormap:0.1.0 \
  anchor_map --config synthetic_rds --out-dir results/demo --threads 2
#> ... C5_sub0 -> anthro [sharp] ...      (results written to ./results/demo)
```

Then point it at your own inputs — see [Running on your own data](#running-on-your-own-data) below.
(The repository must be reachable from your machine.)

## Quick start

A self-contained synthetic example ships with the tool.

### R

```r
library(anchormap)

# score each cluster's anchoring (also runs the reliability z-sweep)
run_anchormap("synthetic_rds", out_dir = "results/demo", threads = 2)

# render the figures
run_plots("synthetic_rds_plots", out_dir = "results/demo/figures")
```

### Docker

```bash
# score each cluster's anchoring (also runs the reliability z-sweep)
docker run --rm -v "$PWD:/work" -w /work anchormap:0.1.0 \
  anchor_map --config synthetic_rds --out-dir results/demo --threads 2

# render the figures
docker run --rm -v "$PWD:/work" -w /work anchormap:0.1.0 \
  plot_anchors --config synthetic_rds_plots --out-dir results/demo/figures
```

## Running on your own data

AnchorMap is config-driven. Two fully-commented templates ship with the package —
[`example_anthro.yaml`](inst/configs/example_anthro.yaml) (single-category track) and
[`example_disease.yaml`](inst/configs/example_disease.yaml) (multi-level disease track). Copy one,
point its `rg_long` / `ontology` / `out_dir` at your files, and run; see
[Configuration](#configuration) for what each field means.

To see a complete run first, a coherent example pair ships with the package — no editing, no external
data (3 clusters × 9 disease traits; the `trait_category → domain` join in action). Run it from your
cloned `AnchorMap` folder (that's where the two `inst/fixtures/...` files live); `--config
example_disease` is a bare name that the tool resolves from its built-in configs:

```bash
anchor_map --config example_disease \
  --rg-long inst/fixtures/example_rg_long.tsv --ontology inst/fixtures/example_ontology.tsv \
  --out-dir results/example_disease --threads 4
#> C0 -> Endocrine-metabolic [sharp] | C1 -> Circulatory [sharp] | C2 -> ambiguous
```

Then for your own data:

```bash
# config-driven (recommended): copy a template, edit the paths
anchor_map --config myrun.yaml --out-dir results/run1 --threads 4

# or keep the template as-is and override the inputs on the CLI (no YAML editing):
anchor_map --config example_anthro.yaml \
  --rg-long data/cluster_trait_rg.tsv --ontology data/ontology.tsv \
  --out-dir results/run1 --threads 4
```

To grab a copy of a template (the shipped configs install alongside the package):

```r
file.copy(system.file("configs/example_disease.yaml", package = "anchormap"), "myrun.yaml")
```

### Via Docker

The shipped configs live **inside the image**, so you refer to them by bare name (`synthetic_rds`,
`example_disease`) — there's nothing to locate or mount. The example *pair's* data TSVs, though, are
files in your cloned `AnchorMap` folder, so run that example from there: `-v "$PWD:/work"` mounts the
current folder (the repo) into the container as `/work`, and `-w /work` makes it the working dir, so
`inst/fixtures/...` resolves:

```bash
docker run --rm -v "$PWD:/work" -w /work anchormap:0.1.0 \
  anchor_map --config example_disease \
    --rg-long inst/fixtures/example_rg_long.tsv --ontology inst/fixtures/example_ontology.tsv \
    --out-dir results/example_disease
```

**For your own data**, the container still can't see your files unless you *mount* their folder. The
`-v host:container` flag maps a folder on your computer to one inside the container — here `./data`
(your inputs, read-only via `:ro`) becomes `/data`, and `./out` (where results land) becomes `/out`:

```bash
docker run --rm \
  -v "$PWD/data:/data:ro" -v "$PWD/out:/out" \
  anchormap:0.1.0 \
  anchor_map --config /data/myrun.yaml --out-dir /out --threads 4
```

The paths in the command (`/data/...`, `/out`) refer to the folders *inside* the container — i.e. the
right-hand side of each `-v` mapping. A template to start from is
[`inst/configs/example_disease.yaml`](inst/configs/example_disease.yaml) in the cloned repo (or copy
one out with the `file.copy(...)` snippet above); put your edited copy in `./data` and mount it.

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

- **TSV route — two files:**
  - a **cluster × trait `rg` long-table** with columns `cluster_label, trait_id, trait_category,
    trait_group, rg, rg_se, p, h2_trait, h2_trait_se, ldsc_converged, negative_h2, status`
    (one row per cluster × trait) — see [`example_rg_long.tsv`](inst/fixtures/example_rg_long.tsv)
    for a small, ready-to-read example (it pairs with
    [`example_ontology.tsv`](inst/fixtures/example_ontology.tsv); see below to run it);
  - optionally a **trait × trait `rg` matrix** (an LDSC `--rg` summary: `p1, p2, rg, …, CONVERGED`)
    as the within-category redundancy source.
- **`.rds` route — a GenomicSEM `ldsc()` object** (`--rds`): AnchorMap standardizes `$S` → `rg`, derives
  `rg_se` by the delta method on `$V`, splits cluster factors from panel traits, and derives both the
  long-table and the trait × trait block — so the rest of the engine is route-agnostic.

Both routes also need an **ontology TSV** mapping each trait to its grouping(s), with an
`anchor_eligible` flag (ineligible groups can be scored but never become an auto-label).

**What gets mapped to what.** For the disease track the chain is, **per cluster**:

```text
long-table trait_category  ──join──►  ontology row  ──►  domain  (──►  icd_chapter)
        (e.g. E4_DM2)                                  (Endocrine-metabolic)
```

For each cluster the engine ranks its `rg` across all traits, groups those traits by `domain` (looked
up from each trait's `trait_category`), and asks which domain is enriched — **that** is the cluster's
anchor. The auto-label is the winning `domain` (or whichever level is set as `primary_level`). So the
column that *drives the anchor* is `trait_category` → `domain`; `trait_id` is the trait's unique
identifier (used to compute the redundancy/VIF correction and to key one row per `(cluster, trait)`),
**not** part of the mapping itself.

The join key is `ontology_key` — `trait_category` for the disease track (above), or `trait_id` for
the anthro/lab track (where traits map directly to a single class level, e.g. `anthro_class`). On the
disease track the levels form a **fine → coarse hierarchy** — `trait_category` (finest, = the join
key) rolls up into `domain`, which rolls up into `icd_chapter` — so `levels: [native, domain,
icd_chapter]` anchors at three zoom levels of the same ontology (`native` aliases the join key).

## Configuration

A run is fully described by one YAML file. The shipped templates
([`example_anthro.yaml`](inst/configs/example_anthro.yaml),
[`example_disease.yaml`](inst/configs/example_disease.yaml)) carry an inline comment on every field;
the table below is the reference. Any field can be overridden on the CLI (e.g. `--out-dir`,
`--threads`, `--z-vector`, `--rds`, `--rg-long`, `--trait-rg`, `--ontology`).

| Group | Field | Meaning |
| --- | --- | --- |
| **Inputs** | `run_label` | name for this run (used in log lines / output labelling) |
| | `rg_long` | path to the cluster × trait `rg` long-table (TSV route) |
| | `trait_rg_matrix` | *(optional)* path to the trait × trait LDSC `--rg` summary (redundancy source) |
| | `rds` | *(optional)* path to a GenomicSEM `ldsc()` `.rds` (`.rds` route; replaces `rg_long`) |
| | `ontology` | path to the ontology TSV |
| | `out_dir` | where TSVs, the log, and `figures/` are written |
| **Reliability gate** | `trait_group` | which `trait_group` rows to score (e.g. `disease`) |
| | `require_ldsc_converged` | drop rows where LDSC did not converge |
| | `drop_negative_h2` | drop rows with negative trait h² |
| | `h2_z_threshold` | keep a trait only if `h2_trait / h2_trait_se >` this (default `4.0`) |
| **Ontology** | `ontology_key` | join key — `trait_id` (anthro/lab) or `trait_category` (disease) |
| | `levels` | ontology levels to score (e.g. `[native, domain, icd_chapter]`) |
| | `primary_level` | the level used for the per-cluster auto-label |
| | `min_category_n` | minimum traits in a category for it to be scored |
| **Enrichment** | `rank_variable` | statistic ranked for the Mann–Whitney AUC (`abs_z`) |
| | `permutation_K` | label-permutation draws for `perm_p` (default `2000`) |
| | `random_seed` | RNG seed (determinism) |
| | `vif_correlation` | redundancy source: `auto` → `trait_rg` → `cluster_profile` → `VIF=1` |
| | `vif_coverage_min` | min trait×trait coverage before `auto` falls back from `trait_rg` (`0.5`) |
| | `trait_rg_require_converged` | use only converged trait×trait pairs |
| | `vif_min_rho` | floor applied to correlations entering the VIF |
| **Over-representation** | `hit_abs_rg` | `\|rg\|` threshold for a trait to count as a Fisher "hit" |
| | `hit_bonferroni` | additionally require Bonferroni-significant `p` for a hit |
| **Auto-label** | `label_q_max` | max FDR `q` for a category to be label-eligible (`0.05`) |
| | `label_auc_min` | min AUC for a category to be label-eligible (`0.60`) |
| **Anchor shape** | `shape_margin_sharp` | AUC margin (top vs next) above which the anchor is *sharp* |
| | `shape_margin_diffuse` | margin below which it is *diffuse* |
| | `shape_focus_diffuse` | spread of significant domains above which it is *diffuse* |

The **figures** config is separate ([`example_plots.yaml`](inst/configs/example_plots.yaml)): it lists
one or more `tracks` (each a `name`, `level`, and the `scores`/`labels` TSVs the engine wrote) plus
plot knobs (`top_k`, `rg_cap`) and the cross-cluster specificity gate (`q_sig`, `spec_rg_floor`,
`spec_min_clusters`, overridable via `--q-sig` / `--rg-floor` / `--min-clusters`).

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
