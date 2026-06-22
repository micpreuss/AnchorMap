# AnchorMap — Analysis Design Document

> Generalizing the `cluster_anchoring` method into a portable, reproducible, Dockerized R + Nextflow tool.
> **Status:** design (pre-implementation). **Author:** preuss_micha@icloud.com. **Date:** 2026-06-18.
> **Spinoff of:** `UKBB_CLUSTER_GWAS/scripts/cluster_anchoring` (reference engine `anchor_categories.py`).

---

## 1. Overview

**Question.** Given a set of latent **cluster factors** (e.g. gPCA factors of a phenotype cluster) and
their genetic correlations (`rg`) to a panel of traits, **which ontology domain does each cluster
"anchor" to**, with what confidence, and is that anchor *sharp* (one dominant domain) or *diffuse*
(a broad smear)? The score must be **competitive** (in-category vs out-category), **size-aware**, and
**correlation-aware** (it must discount redundant, genetically correlated traits), and it must be
**robust across reliability thresholds**.

**Hypothesis.** Exploratory / methodological — there is no single prior effect to confirm. The
operating expectation is mechanistic: a cluster whose latent factor is genetically driven by a
coherent biological domain will show competitive rank-enrichment of that domain's traits in its `rg`
profile, recoverable independent of the exact heritability-reliability cut.

**Decision supported.** Replaces the manual "eyeball which category a cluster maps to" step with one
reproducible, defensible metric + auto-label per cluster, plus a sensitivity profile that tells the
analyst *how stable* that label is. It matters now because the parent project produces dozens of
clusters across multiple cohorts, and a hand-curated anchoring does not scale or reproduce.

---

## 2. Aims

- **A1 — Faithful R port.** Reproduce the reference `anchor_categories.py` algorithm in R, validated
  bit-for-bit (within tolerance) against the Python output on the documented worked example.
- **A2 — Standardized inputs.** Ingest GenomicSEM / LDSC outputs (`ldsc()` `.rds` S/V objects and LDSC
  `--rg` summaries) and the cluster×trait rg long-table with **no or minimal reformatting**.
- **A3 — Default + fallback redundancy.** Make the **trait×trait rg matrix the default** source for the
  within-category redundancy (n_eff / VIF), with **automatic fallback to the cluster-profile proxy**
  when the matrix is absent or low-coverage.
- **A4 — Parallel sensitivity + delivery.** Run a **z-threshold sensitivity sweep in parallel**
  (multi-CPU `perm_p`), emit detailed TSVs + a step log, and ship a **version-pinned, Nextflow-ready
  Docker image** carrying the project's known container fixes.
- **A5 — Visualization.** Render the anchor profile and cross-cluster specificity as
  **publication-ready figures** (lollipop small-multiples, cluster×category dot-heatmap,
  AUC-vs-coherence diagnostic, specificity heatmap + diagonal reduction), config-driven and headless,
  from the same scored TSVs — porting the reference `plot_anchors.py` / `plot_specificity*.py`
  encodings to R (ggplot2).

---

## 3. Background & rationale

The reference pipeline ([`anchor_categories.py`](../UKBB_CLUSTER_GWAS/scripts/cluster_anchoring/anchor_categories.py))
takes an `rg` long-table (one row per cluster × trait) and, per (cluster, ontology level, category),
computes a competitive Mann–Whitney **AUC** of in-category vs out-category `rg` significance, a
label-**permutation p**, a CAMERA-style **VIF-corrected z** (deflating for within-category genetic
correlation), **BH-FDR** across categories, an inverse-variance **pooled rg** (Fisher-z) + coherence,
and a Fisher **over-representation** odds ratio; it emits a ranked anchor profile, an **auto-label**,
and an **anchor shape** (`sharp`/`focal`/`diffuse`/`weak`).

The method synthesizes several established pieces:

- **Li & Ji (2005)** effective number of independent tests `n_eff` from the eigenvalues of the trait
  correlation matrix — here computed with **`poolR::meff(R, method="liji")`** (Cinar & Viechtbauer,
  *J Stat Soft* 101(1)).
- **CAMERA** competitive-test variance inflation factor (VIF) for correlated features.
- **Mann–Whitney/Wilcoxon AUC** as a size-aware competitive rank-enrichment statistic.
- **Inverse-variance Fisher-z meta-analysis** for the pooled signed magnitude.
- **Benjamini–Hochberg FDR** for multiplicity across competing domains.

**Gap this fills.** The reference engine is Python/numpy, single-track, single-z, serial, and wired to
FinnGen-specific paths. AnchorMap generalizes it to a stand-alone, reproducible tool with standardized
inputs, automatic redundancy-source fallback, a parallel reliability-threshold sweep, and a pinned
container that runs cleanly under Nextflow on Google Batch.

**Canonical worked example.** The two cheat-sheet PDFs (`cluster_anchoring_cheatsheet.pdf`,
`liji_neff_anchor_shape_cheatsheet.pdf`) carry one example end-to-end — cluster **C5_sub0**, anthropometric
track → **"Anthropometric [sharp]"** — and are the positive control for the port (see §8).

---

## 4. Data

> Read real headers; the schemas below are the *contracts* AnchorMap enforces.

### 4.1 Sources & access

| Object | Origin | Path / access |
|---|---|---|
| Reference engine (oracle) | parent project | `UKBB_CLUSTER_GWAS/scripts/cluster_anchoring/anchor_categories.py` |
| cluster×trait rg long-table | `FinnGen_PheWAS_RG` stage (GenomicSEM `ldsc()` per cluster×trait) | `.../FinnGen_PheWAS_RG/results/<run>/rg/cluster_trait_rg_long.tsv` |
| trait×trait rg matrix | FinnGen R12 LDSC `--rg` | `UKBB_CLUSTER_GWAS/data/finngen_rg/finngen_R12_FIN.ldsc.summary.tsv` |
| GenomicSEM `ldsc()` object | parent gPCA stage | `<cluster>.ldsc_output.rds` (`$S`, `$V`, `$I`) |
| ontology maps | curated | `.../cluster_anchoring/ontology/{category_ontology,lab_ontology,anthro_ontology}.tsv` |

### 4.2 Schemas (column contracts)

**Input A — cluster×trait rg long-table (TSV; primary input).** One row per (cluster, trait).

| column | type | meaning / convention |
|---|---|---|
| `cluster_label` | str | latent cluster id (`C0`, `C5_sub0`, …) |
| `trait_id` | str | trait/phenotype id (FinnGen endpoint, OMOP lab id, …) |
| `trait_category` | str | ontology join key (disease track) |
| `trait_group` | str | `disease` / `lab_value` / … (selects the universe) |
| `rg` | float | genetic correlation, signed, expected ∈ [−1, 1] (LDSC can exceed; clipped) |
| `rg_se` | float | SE of `rg`; **must be > 0** to pass the gate |
| `p` | float | p-value of `rg` (used by the Fisher ORA layer) |
| `h2_trait` | float | trait heritability |
| `h2_trait_se` | float | SE of `h2_trait`; **> 0** to pass; defines `h2_z` |
| `ldsc_converged` | bool-str | `TRUE`/`FALSE` |
| `negative_h2` | bool-str | `TRUE`/`FALSE` (dropped when `TRUE`) |
| `status` | str | `success` required |

Booleans are uppercase `TRUE`/`FALSE` strings; numerics coerced (`errors→NA`); missing = empty/`NA`.
Extra provenance columns (`cohort`, `transform`, `batch_id`, `trait_label`, `n_*`, `prevalence_*`, …)
are tolerated and passed through.

**Input B — trait×trait redundancy matrix (LDSC `--rg` summary, TSV).**

```
p1   p2   rg   se   z   p   h2_obs   h2_obs_se   h2_int   h2_int_se   gcov_int   gcov_int_se   CONVERGED
```

Long edge-list, each pair stored once → **symmetrize in code**, clip `rg` to [−1, 1], set diagonal = 1,
filter `CONVERGED==TRUE` (configurable). Only `p1, p2, rg, CONVERGED` are required.

**Input C — GenomicSEM `ldsc()` object (`.rds`; alternative to A/B).** Named list:

| element | shape | use |
|---|---|---|
| `$S` | k×k | genetic **covariance**; standardize `S_Stand = S / √(diag·diagᵀ)` (clamp negative h² → 0, diag → 1; mirror `run_cluster_gpca.R` L414–422) |
| `$V` | q×q, q = k(k+1)/2 | sampling covariance of `vech(S)`; `rg_se` via **delta-method** propagation through the standardization, using the relevant `diag(V)` entries |
| `$I` | k×k | LDSC intercepts (carried, not required by scoring) |

From `$S`/`$V` AnchorMap derives **both** the cluster×trait long-table (cluster-factor rows × trait
columns) **and** the trait×trait matrix (trait rows × trait columns) — so a single standard GenomicSEM
artifact suffices.

**Input D — ontology TSV.** Disease: `trait_category, domain, icd_chapter, kind, anchor_eligible, notes`.
Lab/anthro variants keyed on `trait_id` (`analyte_class` / `anthro_class`). `anchor_eligible=FALSE`
categories stay in the table but **can never become a label** (the forbidden-FP gate). Levels
`native` (= join key), `domain`, `icd_chapter`.

### 4.3 QC & inclusion (per-row reliability gate)

```
keep row iff:
  trait_group == <configured>            AND status == "success"
  AND rg, rg_se finite  AND rg_se > 0
  AND h2_trait_se finite AND h2_trait_se > 0
  AND ldsc_converged                      (if require_ldsc_converged)
  AND NOT negative_h2                     (if drop_negative_h2)
  AND h2_z = h2_trait / h2_trait_se > z   (z = reliability threshold, default 4)
```

Expected counts (FinnGen disease track, reference run): ~15 traits/cluster surviving at `z>4`, ~15
clusters; gate counts are logged per z value.

---

## 5. Scope

**In scope** ✅

- *Data:* ✅ Input A (long TSV) · ✅ Input B (LDSC `--rg`) · ✅ Input C (GenomicSEM `.rds` → derive A+B) · ✅ Input D (ontology, multi-track).
- *Methods:* ✅ full reference pipeline in R (gate → per-trait stats → ontology → n_eff/VIF → AUC → perm_p → pooled rg → ORA → FDR → label → shape) · ✅ trait×trait default + auto-proxy-fallback · ✅ z-threshold sensitivity sweep.
- *Compute:* ✅ version-pinned Docker (`rocker/r-ver:4.4.2`) with project fixes · ✅ multi-CPU `perm_p` + parallel z-sweep · ✅ Nextflow DSL2 process.
- *Visualization:* ✅ R plotting module (`R/plot.R`, ggplot2) porting the reference figures — lollipop small-multiples · cluster×category dot-heatmap · AUC-vs-coherence diagnostic · cross-cluster specificity heatmap + diagonal — config-driven and headless.
- *Deliverables:* ✅ detailed TSVs (primary + sensitivity) · ✅ step log with `FINISHED` statement · ✅ **publication-ready figures** (PNG + PDF per track).

**Out of scope / deferred** ❌

- ❌ Narrative `validation_report.md` / `summary.md` generators — deferred (numbers live in the TSVs + figures).
- ❌ Building the cluster×trait `rg` **upstream** (running GenomicSEM per pair from raw sumstats) — AnchorMap consumes `ldsc()` output, it does not run LDSC.
- ❌ Multi-cohort batch orchestration and R-package/nf-core packaging — future (§14).

---

## 6. Methods & analysis plan

Execution order; equations as implemented. **z** below is the **h²-reliability gate**
(`h2_z = h2_trait/h2_trait_se`), *not* a trait-relevance cut.

1. **Reliability gate** — filter to `h2_z > z` (§4.3).
   *Assumption:* trait h² is reliably non-zero so `rg` is trustworthy. *Check:* per-z gated counts logged.
   *Note:* z changes the universe N ⇒ **every z is a full independent re-run** (the basis of the sweep).

2. **Per-trait statistics**

   ```
   abs_z = |rg / rg_se|
   y     = arctanh(clip(rg, ±0.999))                 # Fisher-z transform
   v     = rg_se² / (1 − clip(rg, ±0.999)²)²          # delta-method variance of y
   ```

3. **Ontology join + in/out labeling** — threshold-free, purely by category membership at each level.

4. **Within-category redundancy** — build trait correlation matrix `R` over the in-category traits.

   ```
   n_eff = poolR::meff(R, method = "liji")            # Li & Ji (2005)
   rho_bar (ρ̄) = mean of finite off-diagonal entries of R
   ```
   *Assumption:* `R` is a valid correlation matrix (symmetric, diag 1; NaN→0 for missing pairs).
   *Check:* `n_eff ≤ N`; coverage % logged. **Source of `R` = trait×trait rg by default** (see §6 fallback).

5. **Competitive AUC** (primary ranker)

   ```
   U   = Σ rank_in − n_in(n_in+1)/2                   # Mann–Whitney U on abs_z
   AUC = U / (n_in · n_out)                            # also signed AUC on z
   ```
   *Assumption:* ranks comparable across in/out; ties handled by average ranks.

6. **Label-permutation p (`perm_p`)** — `K=2000` draws of `n_in` ranks without replacement; cached by
   `n_in`; **parallelized across CPUs**; RNG seeded (`random_seed`).

   ```
   perm_p = (1 + #{null_sum ≥ observed_sum}) / (K + 1)
   ```
   *Check:* should track the analytic VIF-z (`vif_p`); divergence flagged.

7. **CAMERA VIF deflation**

   ```
   VIF   = 1 + (n_eff − 1) · ρ̄
   var0  = (N + 1) / (12 · n_in · n_out)
   z_un  = (AUC − 0.5) / √var0
   vif_z = z_un / √VIF ;   vif_p = Φ̄(vif_z)
   ```

8. **IVW Fisher-z pooled rg + coherence**

   ```
   w = 1/v ;  ȳ = Σ(w·y)/Σ(w)
   pooled_rg = tanh(ȳ) ;  CI = tanh(ȳ ± 1.96·√(VIF/Σw))
   coherence = |mean(rg)| / mean(|rg|)                 # 1.0 = all same sign
   ```

9. **Fisher over-representation (ORA)** — threshold layer (not the ranker):
   `hit = |rg| ≥ hit_abs_rg (0.2) AND p < α`, `α = 0.05/N` (Bonferroni) or `0.05`; 2×2 → `fisher.test(alternative="greater")`.

10. **BH-FDR** of `perm_p` across **eligible** categories within (cluster, level) → `q`; rank by `(q↑, AUC↓)`.

11. **Auto-label gate:** `q < label_q_max (0.05) AND AUC ≥ label_auc_min (0.60) AND vif_z > 0 AND vif_p < 0.05 AND n ≥ min_category_n (3)` → else `ambiguous`.

12. **Anchor shape** (evaluate in order; first match wins):

    ```
    n_sig  = #{eligible domains with q<q_max & AUC≥auc_min}
    margin = AUC₁ − AUC₂
    focus  = 1 / Σ pᵢ²,  pᵢ ∝ max(AUCᵢ − 0.5, 0) over significant domains   # inverse-Simpson
    weak    : n_sig == 0
    sharp   : n_sig == 1  OR  margin ≥ shape_margin_sharp (0.10)
    diffuse : focus ≥ shape_focus_diffuse (3) AND margin < shape_margin_diffuse (0.05)
    focal   : otherwise
    ```

**Trait×trait → proxy auto-fallback (the key generalization).** Today's engine takes a manual
`vif_correlation: trait_rg | cluster_profile` flag and only *warns* below 50 % coverage. AnchorMap makes
it automatic:

```
coverage = fraction of gated traits with ≥1 finite off-diagonal trait×trait rg
if trait×trait provided AND coverage ≥ vif_coverage_min (0.5):  use trait_rg
elif n_clusters ≥ 3:                                            use cluster-profile proxy   # build_trait_profile_corr (needs ≥3 clusters/trait)
else:                                                           VIF = 1 (no deflation) + loud WARN   # anti-conservative; flagged in log
```

The choice and coverage % are recorded in the log. **VIF affects only `vif_p`/`perm_p`-adjacent
inference and the CI width — never the AUC point estimate, ranks, `pooled_rg` point estimate, or
coherence** (an invariant carried over from the reference, asserted in tests).

**Estimands / statistics.** Primary: competitive `AUC` per (cluster, level, category) and the per-cluster
`auto_label` + `anchor_shape`. Secondary: `pooled_rg` (signed magnitude + CI), `coherence`, `odds_ratio`.
**Multiplicity:** BH-FDR across eligible categories within each (cluster, level). **Uncertainty:**
`perm_p`, `vif_p`/`vif_z`, and the Fisher-z `pooled_rg` CI; sensitivity across the z-sweep.

---

## 7. Pipeline / compute architecture

### 7.1 Stage DAG

```
Input (TSV long-table  OR  GenomicSEM .rds  +  trait×trait LDSC summary  +  ontology)
   └─ io.R         read + schema-assert; if .rds → standardize S/V → long-table + trait×trait
        └─ gate.R       reliability gate (per z) + per-trait stats
             └─ redundancy.R   R-matrix (trait_rg default → proxy → VIF=1 fallback) + poolR::meff + ρ̄
                  └─ score.R    AUC · perm_p(parallel) · VIF · pooled rg · ORA           ┐
                       └─ label.R   BH-FDR · rank · auto-label · anchor shape            │ per z
   sensitivity.R: future_lapply over z ∈ {2,3,4,5,6,7}  ──────────────────────────────────┘
        └─ main.R    optparse CLI + YAML config + logging → TSVs + anchormap.log (FINISHED)
             └─ plot.R   ggplot2 figures from the scored TSVs (lollipop · dot-heatmap · AUC-vs-coherence · specificity + diagonal)
```

### 7.2 R package layout (smart & slim, commented)

`R/io.R` (3 readers + asserts) · `R/gate.R` · `R/redundancy.R` (poolR + fallback) · `R/score.R`
(AUC, parallel perm, VIF, pooled rg, ORA) · `R/label.R` (FDR, label, shape) · `R/sensitivity.R`
(z-sweep, parallel) · `R/plot.R` (ggplot2: lollipop · dot-heatmap · AUC-vs-coherence · specificity +
diagonal) · `R/main.R` (CLI/config/log). `config/*.yaml` canonical params · `tests/` fixtures.

### 7.3 Container / environment (version-pinned, Nextflow-ready)

- **Base:** `rocker/r-ver:4.4.2` — pinned (matches the project's GenomicSEM/coloc images); Debian-based,
  so `procps`/`ps` is available for Nextflow trace.
- **System deps:** `apt-get install --no-install-recommends procps libcurl4-openssl-dev libssl-dev
  libxml2-dev cmake git ca-certificates locales && rm -rf /var/lib/apt/lists/*`.
- **R deps (pinned):** dated **Posit Package Manager (P3M) snapshot** for CRAN packages
  (`poolr, data.table, Matrix, future, future.apply, optparse, yaml/jsonlite` + plotting:
  `ggplot2, patchwork, ragg, scales, ggrepel`), and
  `remotes::install_github("GenomicSEM/GenomicSEM@<commit>")` pinned to a commit/tag.
- **Project fixes carried** (from `docker/postgwas/README.md`): `USER root` (write to GCS-FUSE work dirs
  on Google Batch) + `ENTRYPOINT []` (prevent an upstream entrypoint from intercepting Nextflow's
  `.command.run` → exit 126).
- **Build-time smoke test:** `RUN` loads libraries and runs the worked-example fixture, asserting the
  C5_sub0 result — fail the build on regression.
- **Tagging:** `anchormap:0.1.2`, pushed to Artifact Registry; referenced **by version, never `latest`**.

### 7.4 Threading

Nextflow process exports `OMP_NUM_THREADS=OPENBLAS_NUM_THREADS=MKL_NUM_THREADS=VECLIB_MAXIMUM_THREADS=
NUMEXPR_NUM_THREADS=1` (avoid BLAS oversubscription), and passes `--threads ${task.cpus}` →
`data.table::setDTthreads(threads)` + `future::plan("multicore", workers = threads)`. The z-sweep
(`future_lapply` over z) and `perm_p` are the parallel hot spots.

### 7.5 Storage layout

```
results/<run_label>/
  primary/      category_anchor_scores.tsv   cluster_anchor_labels.tsv
  sensitivity/  sensitivity_z_scores.tsv     sensitivity_z_labels.tsv
  figures/      anchor_lollipop_<track>.{png,pdf}        anchor_dotheatmap_<track>.{png,pdf}
                anchor_auc_coherence_<track>.{png,pdf}   anchor_specificity_<track>.{png,pdf}
                anchor_specificity_diagonal_<track>.{png,pdf}
  logs/         anchormap.log
```

Nextflow publishes via the `output {}` block + `outputDir` (NOT `publishDir` in-process), per the
project's CLAUDE.md convention; `workDir`/`outputDir` are GCS paths under the `gcp` profile.

---

## 8. Validation & controls

**Positive control (shipped fixture + CI gate).** Worked example **C5_sub0**, anthropometric track:

| quantity | expected |
|---|---|
| `pooled_rg` | ≈ **0.25** (CI ≈ [0.19, 0.30]) |
| `AUC` (abs) | ≈ **0.92** |
| `n_eff` (BMI/WEIGHT/HEIGHT, ρ̄=0.50) | **2.00** |
| `vif_z` / `vif_p` | ≈ 1.77 / ≈ 0.038 |
| `auto_label` | **Anthropometric** |
| `anchor_shape` | **sharp** (`n_sig=1`, `margin = 0.917 − 0.722 = 0.19`) |

**Cross-language oracle (one-time gate).** R output must equal the Python `anchor_categories.py` output
on the fixture within tolerance (point estimates exact; `perm_p` within MC error at fixed seed).

**Negative control.** A domain with no expected enrichment for C5_sub0 (e.g. Neoplasm) → `AUC ≈ 0.5`/0,
`q` not significant. **Forbidden-FP test:** an `anchor_eligible=FALSE` category (Quantitative / Lab)
must **never** surface as the `auto_label`.

**Sensitivity check.** `auto_label` and `anchor_shape` should be **stable across z ∈ {3,4,5}**; the
sweep TSV reports the z values at which a label flips (instability flag).

**Schema / sanity checks.** Required columns present; boolean/NA encodings parse; trait×trait matrix is
symmetric with unit diagonal; `n_eff ≤ N`; trait×trait coverage % computed and logged; VIF ≥ 1.

---

## 9. Outputs & deliverables

| file | format | key columns / content | consumer |
|---|---|---|---|
| `primary/category_anchor_scores.tsv` | TSV | `cluster_label, level, category, eligible, n, n_eff, n_hit, rho_bar, vif, auc_abs, auc_signed, perm_p, vif_z, vif_p, pooled_rg, pooled_rg_ci_lo/hi, coherence, mean_abs_rg, mean_signed_rg, odds_ratio, fisher_p, q, rank` (at primary z=4) | analyst / downstream join |
| `primary/cluster_anchor_labels.tsv` | TSV | `cluster_label, auto_label, anchor_shape, anchor_margin, anchor_focus, n_sig_domains, top_auc, top_q, top_pooled_rg, top_coherence, profile` | reporting |
| `sensitivity/sensitivity_z_scores.tsv` | TSV | as scores + a `z_threshold` column, stacked across the sweep | sensitivity analysis |
| `sensitivity/sensitivity_z_labels.tsv` | TSV | labels + `z_threshold` + a `label_stable` flag | robustness audit |
| `figures/anchor_{lollipop,dotheatmap,auc_coherence}_<track>.{png,pdf}` | PNG + PDF | per-cluster **lollipop** (AUC x-axis, stem = signed `pooled_rg`, alpha = coherence, ring = `q<0.05`, `n` annotated); cluster×category **dot-heatmap** (size = AUC, colour = signed `pooled_rg`, black edge = `q<0.05`, ★ = auto-label); **AUC-vs-coherence** diagnostic (sign-split classes sit top-left) | reporting / publication |
| `figures/anchor_specificity{,_diagonal}_<track>.{png,pdf}` | PNG + PDF | cross-cluster **specificity heatmap** (within-category z of signed `pooled_rg` across clusters, significance-gated) + its **diagonal** reduction (single most-distinctive significant cell per cluster) | cross-cluster distinctiveness |
| `logs/anchormap.log` | text | timestamped steps: config echo, package versions + input file hashes + image tag, per-z gate counts, fallback decision + coverage %, per-z progress, figure manifest, and a final **`FINISHED`** line (status, elapsed, output manifest) | provenance / ops |

**Deferred:** `validation_report.md`, `summary.md` — see §14.

---

## 10. Reproducibility

- **Pinned environment:** `rocker/r-ver:4.4.2` base + dated P3M CRAN snapshot + GenomicSEM pinned commit;
  image referenced by version tag, never `latest`.
- **Determinism:** `random_seed` (default 1) fixes `perm_p`; same seed + same K → same p.
- **Canonical params:** a single YAML config is the source of truth (thresholds, z-vector, K, FDR α,
  AUC cut, ORA `|rg|` cut, shape thresholds, `vif_coverage_min`, threads).
- **Provenance capture:** the log records package versions, input file SHA hashes, the image tag, the
  fallback decision, and the resolved config.

---

## 11. Success criteria

- ✅ Worked example recovered (C5_sub0 → Anthropometric [sharp]; values in §8).
- ✅ R output == Python reference on the fixture (within tolerance).
- ✅ Trait×trait is the default redundancy source; auto-fallback to proxy (and to VIF=1) works and is logged.
- ✅ z-sweep runs in parallel; both sensitivity TSVs emitted with a stability flag.
- ✅ Figures render headless (lollipop, dot-heatmap, AUC-vs-coherence, specificity + diagonal) for the anthro + disease tracks, with the reference encodings (AUC and `pooled_rg` as distinct channels).
- ✅ Image builds; `ps` present; runs under Nextflow as non-root-safe `USER root` on GCS without exit-126/permission errors.
- ✅ All input schema contracts (§4.2) enforced with clear errors on violation.
- ✅ The log terminates with a `FINISHED` statement and an output manifest.

---

## 12. Phases / milestones

**Phase 1 — R engine port + fixture.**
Goal: faithful R reimplementation of the reference pipeline.
Deliverables: ✅ `gate.R`/`redundancy.R`/`score.R`/`label.R`; ✅ C5_sub0 fixture + test.
Gate: R == Python on the fixture (§8 values).

**Phase 2 — Input generalization.**
Goal: standardized ingestion + auto-fallback.
Deliverables: ✅ `io.R` (long-TSV, LDSC `--rg`, GenomicSEM `.rds` readers + schema asserts); ✅ trait×trait-default→proxy→VIF=1 fallback.
Gate: GenomicSEM `.rds`-derived rg/rg_se match per-pair `ldsc()` within tolerance; fallback unit-tested.

**Phase 3 — Sensitivity + parallelism.**
Goal: parallel z-sweep and `perm_p`.
Deliverables: ✅ `sensitivity.R` (z-vector, `future`); ✅ sensitivity TSVs + stability flag; ✅ threaded `perm_p`.
Gate: results invariant to thread count; z=4 slice == Phase-1 primary output.

**Phase 4 — Visualization. ✅ BUILT.**
Goal: publication-ready anchor + cross-cluster specificity figures from the scored TSVs.
Deliverables: ✅ `R/plot.R` + CLI `R/plot_anchors.R` + `configs/carey_rint15_plots.yaml` — lollipop
small-multiples (AUC x-axis, stem = signed `pooled_rg`, alpha = coherence, ring = `q<0.05`, ★ =
auto-label), cluster×category dot-heatmap (size = AUC, colour = signed `pooled_rg`, ★ = auto-label),
AUC-vs-coherence diagnostic, cross-cluster specificity heatmap + its diagonal reduction; PNG + PDF.
**New deps:** `ggplot2`, `patchwork`, `scales`, `ggrepel` (+ optional `ragg`; cairo fallback).
Gate: ✅ figures render headless (ragg if present else cairo) for the anthro + disease tracks; encodings
match the reference (AUC and `pooled_rg` are distinct channels — visible at C2_sub0/C2_sub1, blue at high
AUC); the cross-cluster specificity z (the only recomputation) is **byte-identical to the Python
`cluster_distinctive_categories.tsv`** on the disease track. *(Lab track activates once a
`results/carey_rint15_lab/` run exists — stubbed/commented in the plot config.)*

**Phase 5 — Docker + Nextflow.**
Goal: shippable, reproducible container + process.
Deliverables: ✅ pinned `Dockerfile` (procps, USER root, ENTRYPOINT [], P3M snapshot, smoke test); ✅ DSL2 `ANCHORMAP` process + `nextflow.config`; ✅ log with `FINISHED`.
Gate: image builds; smoke test passes; end-to-end Nextflow run on the fixture produces all outputs (primary + sensitivity TSVs, figures, log).

---

## 13. Risks & mitigations

| risk | mitigation |
|---|---|
| **GenomicSEM `V`-matrix delta-method `rg_se` mis-indexed** → wrong SEs feeding every downstream stat | unit-test derived `rg`/`rg_se` against per-pair `ldsc()` output; assert on the fixture |
| **R↔Python numeric drift** (rank-tie handling, RNG stream differs) | fix seed + tolerance; keep Python as a one-time oracle; use average-rank ties in both |
| **Proxy fallback undefined for few clusters** (proxy needs ≥3 clusters/trait) → silently no deflation | explicit `VIF=1` branch with a loud WARN in the log; documented as anti-conservative |
| **BLAS oversubscription under Nextflow** (R + OpenBLAS each spawn threads) | export all `*_NUM_THREADS=1`; drive parallelism only via `future`/`setDTthreads` from `task.cpus` |
| **Dependency drift** breaks reproducibility | dated P3M snapshot + pinned GenomicSEM commit + version-tagged image (never `latest`) |
| **Low trait×trait coverage** silently weakens VIF | compute + log coverage %; auto-fallback at `vif_coverage_min`; surface in output |

---

## 14. Future / follow-ups

- Narrative `validation_report.md` / `summary.md` generators (the numbers already live in the
  TSVs, and the figures are produced in Phase 4).
- Multi-cohort batch orchestration; packaging as a proper R package + an nf-core-style module.
- Configurable rank variable (`abs_z` vs `abs_rg`) and additional `n_eff` methods exposed
  (`poolR::meff` also offers `nyholt`, `gao`, `galwey`) for sensitivity.
- Liability-scale handling and per-category z thresholds if reliability varies by domain.

---

## 15. Appendix

**Reference code & data**
- Reference engine: `UKBB_CLUSTER_GWAS/scripts/cluster_anchoring/anchor_categories.py` (+ `docs/approach.md`, `ontology/`).
- Reference figures (ported in Phase 4): `UKBB_CLUSTER_GWAS/scripts/cluster_anchoring/{plot_anchors,plot_specificity,plot_specificity_diagonal}.py` (+ `docs/figures_guide.md`; encoding rationale in `docs/approach.md` §4.4/§4.7/§5.2).
- Container fixes: `UKBB_CLUSTER_GWAS/docker/postgwas/README.md` & `Dockerfile`; `docker/genomicsem/Dockerfile` (rocker base + `remotes` pattern); root `Dockerfile` (procps pattern).
- Nextflow conventions: `UKBB_CLUSTER_GWAS/scripts/pipeline/nextflow.config`; `scripts/genomic_sem/UKBB_Carey/main.nf` (threading); project `CLAUDE.md`.
- GenomicSEM standardization reference: `run_cluster_gpca.R` L414–422 (`S_Stand`).

**Worked example**
- `cluster_anchoring_cheatsheet.pdf` — rg long-table → AUC & anchor shape (C5_sub0).
- `liji_neff_anchor_shape_cheatsheet.pdf` — Li & Ji `n_eff` (PART A) + anchor-shape ruleset (PART B).

**Literature**
- Li & Ji (2005), *Heredity* 95:221 — doi:10.1038/sj.hdy.6800717 (effective number of tests).
- Cinar & Viechtbauer — `poolr`, *J Stat Soft* 101(1) (`meff`).
- Wu, Smyth — CAMERA competitive gene-set test (VIF).
