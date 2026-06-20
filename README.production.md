# AnchorMap — production / development notes

> **Internal doc.** This is the detailed engineering record (provenance, phase history, validation
> against the reference implementation). The public-facing front page is [README.md](README.md); this
> file is linked from [CLAUDE.md](CLAUDE.md) and is where ongoing dev notes live.
>
> **Restructured into an installable R package (`anchormap`).** The engine is now `R/*.R` exporting
> `run_anchormap()` / `run_plots()` / `run_sensitivity()`; CLI entry points are the `anchor_map` /
> `plot_anchors` shims (or `inst/scripts/*.R`), no longer a top-level `anchor_map.R`. Example configs +
> ontologies + fixtures live under `inst/`; the machine-specific Carey/FinnGen configs moved to the
> gitignored `local/configs/`. The run commands below have been updated to this package layout;
> [README.md](README.md) is the public quick-start.

Top-level index for the project. AnchorMap is a portable, reproducible **R + Docker** tool (with a
Nextflow container-validation harness) that generalizes the `cluster_anchoring` method: given latent **cluster factors** and their
genetic correlations (`rg`) to a trait panel, it scores — competitively, size-aware, and
correlation-aware — **which ontology domain each cluster anchors to**, how confidently, and whether
the anchor is *sharp* or *diffuse*.

Read [ANALYSIS_DESIGN.md](ANALYSIS_DESIGN.md) (the source-of-truth ADD) and [CLAUDE.md](CLAUDE.md)
(conventions + load-bearing schemas) first; then the engine modules in the order under
[Workflow order](#workflow-order).

## Purpose

Brownfield port of a validated Python/numpy reference
(`../UKBB_CLUSTER_GWAS/scripts/cluster_anchoring/anchor_categories.py`) into a stand-alone R engine
that consumes standardized GenomicSEM/LDSC inputs — either the two-TSV route **or** a GenomicSEM
`ldsc()` `.rds` directly (Phase 2) — selects its within-category redundancy source automatically
(trait×trait `rg` → cluster-profile proxy → no-deflation `VIF=1`), and emits detailed TSVs + a step
log. The headline requirement is **bit-for-bit cross-language parity** on the deterministic outputs,
not unit tests alone.

## Status

**Phases 1–5 complete — the R engine is a validated drop-in for the Python reference, ingests a
GenomicSEM `.rds` directly, runs a parallel reliability-threshold sensitivity sweep, renders the
publication-ready figures, and ships as a pinned, self-validating Docker image** (with a Nextflow
harness that proves the image runs under Nextflow):

1. **Phase 1 — R engine port + fixture** ✅ — the `R/` modules + the `anchor_map` CLI
   ([inst/scripts/anchor_map.R](inst/scripts/anchor_map.R)), validated by cross-language parity on the
   anthro + disease tracks.
2. **Phase 2 — input generalization** ✅ — GenomicSEM `ldsc()` `.rds` ingestion
   ([R/ingest_rds.R](R/ingest_rds.R): rg = S/√(diag·diag); `rg_se` via the exact delta-method on the
   3×3 `V`-submatrix, column-major vech indexing) via `--rds`/`cfg$rds`, plus the
   `vif_correlation: auto` redundancy auto-fallback (trait×trait → cluster-profile proxy → `VIF=1`).
   Validated by the [tests/testthat/](tests/testthat/) suite (delta-method numeric-diff, `.rds`↔TSV
   round-trip, fallback branches, VIF-invariance); Phase-1 oracle parity preserved byte-for-byte.
3. **Phase 3 — sensitivity + parallelism** ✅ — parallel h²-reliability **z-sweep**
   ([R/sensitivity.R](R/sensitivity.R)): `run_sensitivity` re-runs the whole engine at each z in
   `cfg$z_vector` (default `{3,4,5}`, the primary z always folded in) in parallel via `future.apply`,
   emitting `sensitivity_z_scores.tsv` + `sensitivity_z_labels.tsv` (with a per-cluster `label_stable`
   flag) alongside the unchanged primaries. `--z-vector` / `--threads` on the CLI. Validated by
   the [tests/testthat/](tests/testthat/) suite (primary-slice parity incl. `perm_p`, thread-invariance,
   `label_stable`, gate monotonicity). Determinism is engineered: each z-task re-seeds with
   `random_seed` **and pins the RNG kind to Mersenne-Twister** (`future.seed` flips it to L'Ecuyer),
   so the z = `h2_z_threshold` slice is byte-identical to the Phase-1/2 single-z primaries and the
   sweep is thread-count- and backend-invariant; `perm_p` stays serial to protect that parity.
4. **Phase 4 — visualization** ✅ — publication-ready figures via [R/plot.R](R/plot.R) (ggplot2) + the
   `plot_anchors` CLI ([inst/scripts/plot_anchors.R](inst/scripts/plot_anchors.R)) + a plots config:
   per-track lollipop small-multiples, a cluster×category dot-heatmap, an AUC-vs-coherence diagnostic,
   and the cross-cluster specificity heatmap + its diagonal reduction — PNG + PDF, config-driven and
   headless (ragg if present, else cairo), reading only the scored TSVs. The four channels stay distinct
   (AUC = x/size, signed `pooled_rg` = diverging colour, coherence = alpha, `q<q_sig` = ring/mask) so the
   AUC↔rg divergence at sign-split classes survives. The only recomputation — the cross-cluster
   specificity z — is **byte-identical to the Python reference's `cluster_distinctive_categories.tsv`** on
   the disease track (15/15 clusters). Ports `plot_anchors.py` / `plot_specificity{,_diagonal}.py`.
5. **Phase 5 — Docker + Nextflow harness** ✅ — the pinned, self-validating **Docker image** is the tool
   and primary run interface ([docker/Dockerfile](docker/Dockerfile)): `rocker/r-ver:4.6.0` (= the
   validated host R) + a single dated P3M snapshot reproducing the validated `future.apply 1.20.2` /
   `ggplot2 4.0.3`; `procps` / `USER root` / `ENTRYPOINT []` fixes; **no GenomicSEM** (the engine reads
   `ldsc()` `.rds` with base `readRDS`); two build-time self-tests (synthetic `.rds` → C5_sub0 anthro
   [sharp]; figure render). **Nextflow is a container-validation harness, not orchestration**
   ([nextflow/main.nf](nextflow/main.nf)): a single `ANCHORMAP_SMOKE` process — `test` (local) is the CI
   gate (entrypoint / procps / output capture), `gcp` (Google Batch, spot) is a one-time `USER root` /
   GCS-FUSE check. **Diverges from ADD §7.3 deliberately**: base `4.6.0` not `4.4.2` (match the validated
   env), Nextflow scoped to validation. See [docker/README.md](docker/README.md).

The project is **under git** (GitHub: `micpreuss/AnchorMap`, public); the vendored
`claude-science-scaffold/` subdir is gitignored (it is its own repo).

## Inputs and outputs

- **Inputs** — two interchangeable routes converging on the same engine: the two TSVs (A + B) **or** a
  single GenomicSEM `.rds` (C):
  - **A — cluster×trait `rg` long-table (TSV):** required `(cluster_label, trait_id, trait_category,
    trait_group, rg, rg_se, p, h2_trait, h2_trait_se, ldsc_converged, negative_h2, status)`; one row per
    (cluster, trait). The engine **recomputes** any precomputed `z`/`abs_rg`.
  - **B — trait×trait `rg` (LDSC `--rg` summary):** `(p1, p2, rg, …, CONVERGED)` long edge-list →
    symmetrized, clipped to [−1, 1], diag = 1. Default redundancy source.
  - **C — GenomicSEM `ldsc()` `.rds`** (Phase 2; `--rds`/`cfg$rds`): list `$S` (genetic covariance),
    `$V` (sampling cov of vech(S), column-major lower triangle), `$I` (intercepts). [R/ingest_rds.R](R/ingest_rds.R)
    standardizes `S`→rg, derives `rg_se` by delta-method, partitions cluster factors vs panel traits
    (`cluster_factor_pattern`, default `^C[0-9]`), and emits both the long-table (A) and the trait×trait
    block (B) so the rest of the engine is route-agnostic.
  - **D — ontology TSV:** disease joins on `trait_category`; anthro/lab join on `trait_id`.
    `anchor_eligible=FALSE` categories may be scored but never labeled (the forbidden-FP gate).
- **Terminal outputs** (`results/<run_label>/`):
  - `category_anchor_scores.tsv` — per `(cluster_label, level, category)`: `n, n_eff, n_hit, rho_bar,
    vif, auc_abs, auc_signed, perm_p, vif_z, vif_p, pooled_rg [ci], coherence, odds_ratio, fisher_p, q, rank`.
  - `cluster_anchor_labels.tsv` — per cluster: `auto_label, anchor_shape, anchor_margin, anchor_focus,
    n_sig_domains, top_*`, profile.
  - `sensitivity_z_scores.tsv` / `sensitivity_z_labels.tsv` *(Phase 3)* — the two tables above stacked
    across the z-sweep (each + a `z_threshold` column; labels + a per-cluster `label_stable` flag). The
    z = `h2_z_threshold` slice equals the primaries byte-for-byte.
  - `anchormap.log` — timestamped steps (incl. per-z gate counts + the `label-stable` summary) ending
    in a `FINISHED` line (status, elapsed, output manifest).
  - `figures/` *(Phase 4)* — per-track lollipop small-multiples, cluster×category dot-heatmap,
    AUC-vs-coherence diagnostic, and cross-cluster specificity heatmap + diagonal (PNG + PDF), plus
    `cluster_distinctive_categories.tsv` (`track, cluster_label, distinctive_category, spec_z, pooled_rg,
    runner_up`). Written by the `plot_anchors` CLI from the scored TSVs, into `results/<run>/figures/`.

  See [CLAUDE.md](CLAUDE.md) for the full column contracts, units, sign, and rounding rules.

## Canonical datasets and runs

| Dataset | Canonical run | Notes |
| --- | --- | --- |
| Carey RINT-15 **anthro** track | [local/configs/carey_rint15_anthro.yaml](local/configs/carey_rint15_anthro.yaml) | **The positive control.** C5_sub0 → "Anthropometric [sharp]". Single `anthro_class` level vs the full disease universe. |
| Carey RINT-15 **disease** track | [local/configs/carey_rint15.yaml](local/configs/carey_rint15.yaml) | Comprehensive oracle: all clusters × {native, domain, icd_chapter}, primary level = `domain`. |
| Input `rg` long-table (A) | `../UKBB_CLUSTER_GWAS/scripts/FinnGen_PheWAS_RG/results/carey_rint_tuned_15clusters_neff_max_empirical_covz/rg/cluster_trait_rg_long_with_p.tsv` | 42,795 rows. Clusters scored independently → safe to subset whole clusters for fast fixtures. |
| Input trait×trait `rg` (B) | `../UKBB_CLUSTER_GWAS/data/finngen_rg/finngen_R12_FIN.ldsc.summary.tsv` | FinnGen R12 FIN LDSC `--rg` summary; 100% disease-endpoint coverage. |
| Synthetic `.rds` route smoke | [inst/configs/synthetic_rds.yaml](inst/configs/synthetic_rds.yaml) | End-to-end Phase-2 `.rds` ingestion on a self-contained synthetic fixture ([inst/fixtures/](inst/fixtures/), built by `make_synthetic_ldsc.R`); `vif_correlation: auto` → trait_rg @ 100% coverage. |

Configs reuse the parent `cluster_anchoring` configs verbatim — only paths are adapted (big inputs →
absolute parent-repo paths; `ontology`/`out_dir` → AnchorMap-relative).

## Workflow order

The engine is the `anchor_map` CLI ([inst/scripts/anchor_map.R](inst/scripts/anchor_map.R)) over the
`R/` modules below; read them in execution order to catch up on the method (they port the reference 1:1,
single z).

1. [R/io.R](R/io.R) — config (identical YAML to the parent) + standardized-TSV readers + schema asserts.
2. [R/ingest_rds.R](R/ingest_rds.R) *(Phase 2)* — GenomicSEM `.rds` route: `vech_index`,
   `standardize_S`, delta-method `rg_se_matrix`, factor/panel `partition_S`, and `rds_to_long` /
   `rds_to_trait_rg` (emit the Input-A/B contracts). Skipped on the TSV route.
3. [R/gate.R](R/gate.R) — reliability gate `h2_z = h2_trait/h2_trait_se > z`; per-trait stats
   (`abs_z`, Fisher-z `y = atanh(rg)`, delta-var `v = rg_se²/(1−rg²)²`); ontology in/out join.
4. [R/redundancy.R](R/redundancy.R) — within-category redundancy: `n_eff` (Li & Ji 2005 via clipped
   eigenvalues; `poolr::meff` cross-check) and mean pairwise `ρ̄`. `select_corr_source` picks the
   matrix per `vif_correlation`: explicit `trait_rg`/`cluster_profile` (Phase-1, verbatim) or `auto`
   (Phase 2: trait×trait if coverage ≥ `vif_coverage_min`, else proxy if ≥3 clusters, else `VIF=1`+WARN).
5. [R/score.R](R/score.R) — competitive Mann–Whitney **AUC**, label-permutation **perm_p**, CAMERA
   **VIF** deflation, IVW Fisher-z **pooled_rg** + coherence, Fisher **ORA**.
6. [R/label.R](R/label.R) — **BH-FDR** `q`, the **auto_label** gate, and **anchor_shape**
   (weak / sharp / diffuse / focal).
7. [R/sensitivity.R](R/sensitivity.R) *(Phase 3)* — wraps steps 3–6 in `score_at_z` (one full re-run
   per reliability threshold) and maps it over `cfg$z_vector` via `parallel_lapply` (`future.apply`,
   multicore/sequential); stacks the two tables with `z_threshold` and flags per-cluster `label_stable`.
8. [R/plot.R](R/plot.R) *(Phase 4)* — figure builders consuming only the scored TSVs: `natural_order` /
   `leaf_order` (cluster row + category column ordering), `specificity` + `distinctive_table` (the
   cross-cluster z), and `fig_lollipops` / `fig_dotheatmap` / `fig_scatter` / `fig_specificity` /
   `fig_diagonal` with diverging RdBu (rg) / PuOr (specificity) scales and a `save_fig` (ragg→cairo,
   PNG+PDF). Driven by the separate `plot_anchors` CLI ([inst/scripts/plot_anchors.R](inst/scripts/plot_anchors.R))
   — it does **not** re-run the engine.

Drivers:

- **Engine** — the `anchor_map` CLI ([inst/scripts/anchor_map.R](inst/scripts/anchor_map.R)):
  `anchor_map --config <yaml|name> [--out-dir DIR] [--threads N] [--rds <file>] [--z-vector 3,4,5]` → the
  two primary TSVs + the two sensitivity TSVs + `anchormap.log`.
- **Figures** — the `plot_anchors` CLI ([inst/scripts/plot_anchors.R](inst/scripts/plot_anchors.R)):
  `plot_anchors --config <plots.yaml|name> [--out-dir DIR] [--q-sig N] [--rg-floor N] [--min-clusters N]`
  → `results/<run>/figures/` (PNG+PDF + the distinctive TSV).

## Data flow at a glance

```text
rg long-table (TSV) ─┐
LDSC --rg summary  ──┼─► io.R ─► gate.R ─────► redundancy.R ──────► score.R ───────────────► label.R ─► category_anchor_scores.tsv
ontology TSV  ───────┘  (read)   (h2_z>z;       (n_eff, ρ̄; corr     (AUC, perm_p, VIF,        (BH-FDR q,  cluster_anchor_labels.tsv
GenomicSEM .rds ─► ingest_rds.R   per-trait       source: trait_rg→    IVW pooled_rg,           auto_label, anchormap.log
  (S,V,I → A + B)                 stats)          proxy→VIF=1 auto)    Fisher ORA)              shape)
```

- **Must NOT cross:** `anchor_eligible=FALSE` categories (Quantitative/Lab/administrative) are scored
  but can **never** become a cluster's `auto_label` — the forbidden-FP gate.
- **VIF affects only** `vif_z`/`vif_p`/CI width — never the AUC, ranks, `pooled_rg` point estimate, or
  coherence.
- **Sign is load-bearing:** `rg` ∈ [−1, 1] is signed and the sign must survive pooling (coherence
  depends on it).

## Requirements

R ≥ 4.4 with `data.table`, `yaml`, `future`, **`future.apply` ≥ 1.20** (required; the last two drive the
Phase-3 z-sweep — note 1.11.x has a `future.globals` regression that breaks the sweep) and `poolr`
(optional, `n_eff` cross-check). The Phase-4 figures additionally need `ggplot2`, `patchwork`, `scales`,
`ggrepel` (and optionally `ragg` — `R/plot.R` falls back to cairo when it is absent). The CLI uses
`optparse`; `testthat` is needed only to run the test suite. **Or skip host R entirely and use the pinned
image** (`docker run anchormap:0.1.1 …`), which carries the exact validated versions (R 4.6.0,
`future.apply 1.20.2`, `ggplot2 4.0.3`).

## Run

```bash
R CMD INSTALL .                                                 # install the anchormap package (host)

Rscript inst/scripts/anchor_map.R  --config local/configs/carey_rint15_anthro.yaml --out-dir results/carey_rint15_anthro --threads 4  # TSV route
Rscript inst/scripts/anchor_map.R  --config local/configs/carey_rint15.yaml --rds <ldsc_output.rds> --out-dir results/carey_rint15      # .rds route (override input)
Rscript inst/scripts/anchor_map.R  --config local/configs/carey_rint15.yaml --z-vector 2,3,4,5,6 --out-dir results/carey_rint15         # override the z-sweep
Rscript inst/scripts/plot_anchors.R --config local/configs/carey_rint15_plots.yaml --out-dir results/carey_rint15/figures               # Phase-4 figures

Rscript inst/scripts/anchor_map.R  --config synthetic_rds --out-dir results/synthetic_rds   # .rds-route end-to-end smoke (shipped config + fixture)
Rscript -e 'testthat::test_local()'                            # full testthat suite (analytic + parity + sensitivity + CLI)
bash    validation/run_oracle.sh                               # full cross-language parity vs Python
```

(A source install does not add PATH shims; use the bundled `Rscript inst/scripts/*.R` wrappers on the
host. The pinned image puts `anchor_map` / `plot_anchors` on `PATH`. `--config` accepts a YAML path
**or** a bare shipped-config name resolved from `inst/configs/`.)

### Phase 5 — pinned image + Nextflow harness

```bash
# Build the image (THE tool). The build runs two self-tests; it fails if either regresses.
docker build -t anchormap:0.1.1 -f docker/Dockerfile .          # release: add --platform linux/amd64

# Run AnchorMap reproducibly via the image (primary interface; mount cwd as /work):
docker run --rm -v "$PWD:/work" -w /work anchormap:0.1.1 \
  anchor_map --config local/configs/carey_rint15_anthro.yaml --out-dir results/carey_rint15_anthro --threads 4

# Validate the image runs flawlessly UNDER Nextflow (not how you run AnchorMap):
nextflow run nextflow/main.nf -profile test -params-file nextflow/params/test.yaml   # local CI gate
# nextflow run nextflow/main.nf -profile gcp  -params-file nextflow/params/gcp.yaml  # one-time Batch/FUSE check
```

### Publishing a container release

The [GHCR workflow](.github/workflows/publish-container.yml) builds and publishes the image whenever
a semantic version tag such as `v0.1.1` is pushed. The tag must match `Version:` in `DESCRIPTION`; the
workflow publishes both `ghcr.io/micpreuss/anchormap:0.1.1` and `ghcr.io/micpreuss/anchormap:latest`.
It uses GitHub's built-in `GITHUB_TOKEN`, so no registry secret is required.

```bash
git tag -a v0.1.1 -m "AnchorMap 0.1.1"
git push origin v0.1.1
```

After the first successful publication, open the `anchormap` package settings on GitHub and change
its visibility to **Public**. Public GHCR images can be pulled anonymously; GitHub warns that this
visibility change cannot be reversed. Do not move or reuse a released version tag; increment
`DESCRIPTION` and publish a new tag for the next release.

## Validation (cross-language parity vs the Python reference)

- **Positive control:** C5_sub0 anthro → **"Anthropometric [sharp]"** — `auc_abs=0.9164`,
  `pooled_rg=0.2473`, `vif_p=0.03489`, exactly reproduced.
- **Anthro track:** all deterministic columns match to machine precision (Δ≤6e-16); labels/ranks identical.
- **Disease track (1005 rows):** 0 deterministic mismatches; all 14 cluster auto-labels + shapes
  identical (C4→Psychiatric[sharp], C5_sub0→Cardiovascular[diffuse], …); forbidden-FP holds
  (`Quantitative` never labels C5_sub0).
- **`perm_p`/`q`** agree within Monte-Carlo error (numpy and R RNG streams differ; the gate anchors
  on the deterministic `vif_p` + label stability). Near-tie category **ranks** can reorder within MC
  noise without changing any auto-label.
- **Sensitivity (Phase 3):** the z = `h2_z_threshold` slice of `sensitivity_z_scores.tsv` is
  byte-identical to `category_anchor_scores.tsv`; the sweep is invariant to `--threads`; gate counts
  are monotonic in z (anthro 14/15 clusters label-stable across `{3,4,5}`, disease 8/15); no
  `anchor_eligible=FALSE` category labels any cluster at any z.
- **Figures (Phase 4):** the cross-cluster specificity z (`cluster_distinctive_categories.tsv`) is
  byte-identical to the Python reference on the disease track (15/15 clusters — `distinctive_category`
  and `spec_z`); all eight PNG+PDF figures render headless for the anthro + disease tracks; the AUC↔rg
  divergence is preserved (C2_sub0/C2_sub1 render blue at high AUC — strong enrichment, inverse
  direction); the positive control C5_sub0 anthro lollipop is red, ringed (`q<0.05`), and starred
  (auto-label).

## Conventions

Short pointers only — the authoritative text lives in [CLAUDE.md](CLAUDE.md).

- **Config-over-CLI:** all params live in a `--config <yaml|name>`; reuse parent configs unchanged (adjust
  paths only). CLI flags cover output (`--out-dir`/`--run-label`), input overrides (`--rds`, `--rg-long`,
  `--trait-rg`, `--ontology`), and engine controls (`--threads`, `--z-vector`); parsed via `optparse`.
- **Schemas:** the rg long-table, LDSC `--rg`, GenomicSEM `.rds`, ontology, and output column
  contracts are documented in CLAUDE.md → *Data schemas*.
- **Porting gotchas:** scipy `fisher_exact` returns the *sample* OR `(a·d)/(b·c)` (not R's
  conditional-MLE); `perm_p` is not bit-reproducible across languages; `poolr::meff` needs cleaned/PD
  matrices. See CLAUDE.md → *Gotchas*.
- **Compute (built):** the image pins everything (base `rocker/r-ver:4.6.0` by tag + CRAN via one dated
  P3M snapshot reproducing the validated package set); referenced by version tag, never `latest`; carries
  the parent's `procps` / `USER root` / `ENTRYPOINT []` fixes. GenomicSEM is omitted (the engine consumes
  `ldsc()` `.rds`, never runs it). Build context excludes the scaffold + results via `.dockerignore`.
