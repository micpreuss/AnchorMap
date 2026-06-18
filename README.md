# AnchorMap

Top-level index for the project. AnchorMap is a portable, reproducible **R (+ planned Nextflow)**
tool that generalizes the `cluster_anchoring` method: given latent **cluster factors** and their
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

**Phases 1–2 complete — the R engine is a validated drop-in for the Python reference and ingests a
GenomicSEM `.rds` directly.** Later phases are designed in the ADD but not yet built:

1. **Phase 1 — R engine port + fixture** ✅ — `R/` + [anchor_map.R](anchor_map.R), validated by
   cross-language parity on the anthro + disease tracks.
2. **Phase 2 — input generalization** ✅ — GenomicSEM `ldsc()` `.rds` ingestion
   ([R/ingest_rds.R](R/ingest_rds.R): rg = S/√(diag·diag); `rg_se` via the exact delta-method on the
   3×3 `V`-submatrix, column-major vech indexing) via `--rds`/`cfg$rds`, plus the
   `vif_correlation: auto` redundancy auto-fallback (trait×trait → cluster-profile proxy → `VIF=1`).
   Validated by [tests/test_phase2.R](tests/test_phase2.R) (delta-method numeric-diff, `.rds`↔TSV
   round-trip, fallback branches, VIF-invariance); Phase-1 oracle parity preserved byte-for-byte.
3. **Phase 3 — sensitivity + parallelism** *(designed)* — parallel z-threshold sweep and multi-CPU
   `perm_p` via `future`.
4. **Phase 4 — visualization** *(designed)* — publication-ready figures (lollipop small-multiples,
   cluster×category dot-heatmap, AUC-vs-coherence diagnostic, cross-cluster specificity heatmap +
   diagonal) via `R/plot.R` (ggplot2), config-driven and headless, from the scored TSVs.
5. **Phase 5 — Docker + Nextflow** *(designed)* — pinned `rocker/r-ver:4.4.2` image + `ANCHORMAP`
   process. The container/orchestration rows in CLAUDE.md remain *designed*.

The project is **under git** (GitHub: `micpreuss/AnchorMap`, private); the vendored
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
  - `anchormap.log` — timestamped steps ending in a `FINISHED` line (status, elapsed, output manifest).
  - `figures/` *(designed, Phase 4)* — per-track lollipop small-multiples, cluster×category
    dot-heatmap, AUC-vs-coherence diagnostic, and cross-cluster specificity heatmap + diagonal (PNG + PDF).

  See [CLAUDE.md](CLAUDE.md) for the full column contracts, units, sign, and rounding rules.

## Canonical datasets and runs

| Dataset | Canonical run | Notes |
| --- | --- | --- |
| Carey RINT-15 **anthro** track | [configs/carey_rint15_anthro.yaml](configs/carey_rint15_anthro.yaml) | **The positive control.** C5_sub0 → "Anthropometric [sharp]". Single `anthro_class` level vs the full disease universe. |
| Carey RINT-15 **disease** track | [configs/carey_rint15.yaml](configs/carey_rint15.yaml) | Comprehensive oracle: all clusters × {native, domain, icd_chapter}, primary level = `domain`. |
| Input `rg` long-table (A) | `../UKBB_CLUSTER_GWAS/scripts/FinnGen_PheWAS_RG/results/carey_rint_tuned_15clusters_neff_max_empirical_covz/rg/cluster_trait_rg_long_with_p.tsv` | 42,795 rows. Clusters scored independently → safe to subset whole clusters for fast fixtures. |
| Input trait×trait `rg` (B) | `../UKBB_CLUSTER_GWAS/data/finngen_rg/finngen_R12_FIN.ldsc.summary.tsv` | FinnGen R12 FIN LDSC `--rg` summary; 100% disease-endpoint coverage. |
| Synthetic `.rds` route smoke | [configs/synthetic_rds.yaml](configs/synthetic_rds.yaml) | End-to-end Phase-2 `.rds` ingestion on a self-contained synthetic fixture ([tests/fixtures/](tests/fixtures/), built by `make_synthetic_ldsc.R`); `vif_correlation: auto` → trait_rg @ 100% coverage. |

Configs reuse the parent `cluster_anchoring` configs verbatim — only paths are adapted (big inputs →
absolute parent-repo paths; `ontology`/`out_dir` → AnchorMap-relative).

## Workflow order

The engine is a single CLI ([anchor_map.R](anchor_map.R)) that sources the modules below; read them in
execution order to catch up on the method (they port the reference 1:1, single z).

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

Driver: [anchor_map.R](anchor_map.R) — `Rscript anchor_map.R --config <yaml> [--threads N] [--rds <file>]`
→ the two TSVs + `anchormap.log`. *(Phase 3 adds `R/sensitivity.R` for the parallel z-sweep; Phase 4
adds `R/plot.R` for the figures — new list items anchor here.)*

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

R ≥ 4.4 with `data.table`, `yaml` (required) and `poolr` (optional, `n_eff` cross-check). No
`argparse`/`testthat` needed — the CLI and tests use base R.

## Run

```bash
Rscript anchor_map.R --config configs/carey_rint15_anthro.yaml   # TSV route -> results/carey_rint15_anthro/
Rscript anchor_map.R --config configs/carey_rint15.yaml --rds <ldsc_output.rds>   # .rds route (override input)

Rscript tests/fixtures/make_synthetic_ldsc.R                     # build the synthetic .rds fixtures
Rscript anchor_map.R --config configs/synthetic_rds.yaml         # .rds-route end-to-end smoke
Rscript tests/run_tests.R                                        # Phase-1 analytic unit tests (base-R asserts)
Rscript tests/test_phase2.R                                      # Phase-2: delta-method, round-trip, fallback, VIF-invariance
bash    validation/run_oracle.sh                                 # full cross-language parity vs Python
```

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

## Conventions

Short pointers only — the authoritative text lives in [CLAUDE.md](CLAUDE.md).

- **Config-over-CLI:** all params live in a `--config <yaml>`; reuse parent configs unchanged (adjust
  paths only). CLI flags reserved for `--threads` / `--z-vector` (later phases).
- **Schemas:** the rg long-table, LDSC `--rg`, GenomicSEM `.rds`, ontology, and output column
  contracts are documented in CLAUDE.md → *Data schemas*.
- **Porting gotchas:** scipy `fisher_exact` returns the *sample* OR `(a·d)/(b·c)` (not R's
  conditional-MLE); `perm_p` is not bit-reproducible across languages; `poolr::meff` needs cleaned/PD
  matrices. See CLAUDE.md → *Gotchas*.
- **Compute (designed):** pin everything (base image by tag, CRAN via dated P3M snapshot, GenomicSEM
  by commit); reference images by version, never `latest`; carry the parent's `procps` / `USER root` /
  `ENTRYPOINT []` container fixes.
