# AnchorMap

A portable, reproducible **R + Nextflow** tool that generalizes the `cluster_anchoring` method from
the sibling `UKBB_CLUSTER_GWAS` project: given latent **cluster factors** and their genetic
correlations (`rg`) to a trait panel, it scores — competitively, size-aware, and correlation-aware —
**which ontology domain each cluster anchors to**, how confidently, and whether the anchor is *sharp*
or *diffuse*.

**Status: Phase 1 complete (R engine), later phases designed.** The R engine in [R/](R/) +
[anchor_map.R](anchor_map.R) is a validated drop-in for the Python reference (exact deterministic
parity on the anthro + disease tracks; see [README.md](README.md)). Still designed-not-built:
GenomicSEM `.rds` ingestion + auto trait×trait→proxy fallback (Phase 2), the parallel z-sweep
(Phase 3), the **plotting/visualization module** (Phase 4), and the **Dockerfile + Nextflow process**
(Phase 5) — so the container/orchestration rows below remain *designed*. The project is **under git**
(GitHub: `micpreuss/AnchorMap`, private); the vendored `claude-science-scaffold/` subdir is gitignored
(it is its own repo). Read [ANALYSIS_DESIGN.md](ANALYSIS_DESIGN.md) as
the source of truth; the plan for the next phase goes in `.agents/plans/`.

---

## Project type

Five-axis classification (as **designed** — not yet realized in code):

- **Orchestration:** Nextflow DSL2 (a single `ANCHORMAP` process) **+** a local R engine. *Neither built yet.*
- **Compute backend:** Mixed — local R for the engine; Google Batch (spot) for production, mirroring `UKBB_CLUSTER_GWAS`.
- **Domain:** Statistical genetics (genetic-correlation / cluster anchoring; Li & Ji n_eff, CAMERA VIF, Mann–Whitney AUC, IVW Fisher-z, BH-FDR) **and** a reusable **R method/tool package** generalizing a Python reference.
- **Data locality:** External — consumes `UKBB_CLUSTER_GWAS` outputs (sibling-repo paths today; `gs://` under the `gcp` profile later). AnchorMap itself ships only small test fixtures.
- **Reproducibility:** Containers *(designed)* — pinned `rocker/r-ver:4.4.2` + dated Posit P3M CRAN snapshot + GenomicSEM pinned commit. **⚠ Currently nothing is pinned** (no Dockerfile/renv.lock/DESCRIPTION) — see Gaps.

No build system or test suite exists yet. The intended posture: a small R package (`R/*.R` modules +
a `--config <yaml>` CLI) validated by **cross-language parity against the Python reference**, not unit tests alone.

---

## Tech stack (designed)

| Layer | Tools |
|---|---|
| Orchestration | Nextflow DSL2 process `ANCHORMAP`; `output {}` block + `outputDir` (NOT `publishDir`) — inherited from parent |
| Compute | local R (engine); Google Batch spot for production; submit from `nf-head`, not laptop |
| Containers / envs | `rocker/r-ver:4.4.2` (pinned) + `procps` + `USER root` + `ENTRYPOINT []`; image `anchormap:0.1.0` → Artifact Registry `us-central1-docker.pkg.dev/lencz-lab-cogent-1/docker-images/`; **referenced by version tag, never `latest`** |
| Storage | reads sibling `UKBB_CLUSTER_GWAS` files; writes `results/<run_label>/{primary,sensitivity,figures,logs}/` |
| Languages | R ≥4.4 (`poolr`, `data.table`, `Matrix`, `future`/`future.apply`, `yaml`, `argparse`, `testthat`; plotting `ggplot2`/`patchwork`/`ragg`/`scales`); R `remotes` for GenomicSEM |
| Methods | Li & Ji (2005) n_eff via `poolr::meff(R,"liji")`; CAMERA VIF; Mann–Whitney/Wilcoxon AUC; label-permutation null; IVW Fisher-z pooled rg; BH-FDR; anchor-shape ruleset |
| External reference data | FinnGen **R12** FIN LDSC `--rg` summary (trait×trait rg); curated `category/anthro/lab` ontologies; GenomicSEM `ldsc()` S/V objects |

---

## Repository structure

```
AnchorMap/
├── CLAUDE.md                 ← this file (start here)
├── ANALYSIS_DESIGN.md        ← ★ source of truth: the 15-section ADD
├── .agents/plans/            ← build plans
│   └── anchormap-phase1-r-engine-port.md   ← ★ the next thing to build
├── .claude/                  ← injected scaffold (skills, settings); do not treat as project code
├── claude-science-scaffold/  ← vendored source of .claude (own git repo) — not part of AnchorMap
└── (planned, not yet created):
    ├── R/{io,gate,redundancy,score,label,sensitivity,plot,main}.R   ← engine modules (plot.R = Phase 4)
    ├── anchor_map.R          ← CLI entry (--config <yaml>)
    ├── config/*.yaml         ← canonical params (reuse parent configs verbatim)
    ├── tests/{testthat,fixtures}/   ← analytic + oracle-parity fixtures
    ├── validation/           ← R-vs-Python oracle comparator
    ├── docker/Dockerfile     ← pinned rocker image (Phase 5)
    └── nextflow/             ← ANCHORMAP process + nextflow.config (Phase 5)
```

The **reference engine being ported** lives in the sibling project:
[`../UKBB_CLUSTER_GWAS/scripts/cluster_anchoring/anchor_categories.py`](../UKBB_CLUSTER_GWAS/scripts/cluster_anchoring/anchor_categories.py)
— it is the bit-for-bit spec.

---

## Data flow (high-level)

```
GenomicSEM ldsc() .rds  ─┐
  (S genetic-cov, V)     ├─► io.R ─► gate.R ─► redundancy.R ─► score.R ─► label.R ─► TSVs + log
rg long-table (TSV)  ────┤        (z gate)   (n_eff/VIF)    (AUC,perm_p,  (FDR,label,
LDSC --rg summary  ──────┘                    trait_rg|proxy) pooled rg,ORA) shape)
ontology TSV  ───────────┘                                   └─ sensitivity.R sweeps z in parallel ┘
```

- **Two input routes (same engine):** the standardized `rg` long-table TSV **or** a GenomicSEM `ldsc()` `.rds`
  from which both the long-table (rg = S/√(diag·diag); rg_se via **delta-method on V**) and the trait×trait
  matrix are derived. *Phase 1 = TSV route only; the `.rds` route is Phase 2.*
- **Redundancy source:** trait×trait rg is the **default**; auto-fallback to the cluster-profile proxy when
  coverage < `vif_coverage_min` (0.5); else VIF=1 with a loud WARN. *Auto-fallback is Phase 2; Phase 1 honors
  the config `vif_correlation` flag verbatim.*
- **Must NOT cross:** `anchor_eligible=FALSE` categories (Quantitative/Lab/administrative) may be scored but
  can **never** become a cluster's `auto_label` — the forbidden-FP gate. VIF affects only `vif_p`/CI width,
  **never** the AUC, ranks, `pooled_rg` point estimate, or coherence.

---

## Canonical datasets / runs

| Dataset | Canonical run | Notes |
|---|---|---|
| Carey RINT-15 anthro track | `config/carey_rint15_anthro.yaml` (= parent `cluster_anchoring/configs/carey_rint15_anthro.yaml`) | **The positive control.** C5_sub0 → "Anthropometric [sharp]". Single category vs the full disease universe. |
| Carey RINT-15 disease track | `config/carey_rint15.yaml` | Comprehensive oracle: all clusters × {native,domain,icd_chapter} × many categories. |
| Input rg long-table | `../UKBB_CLUSTER_GWAS/scripts/FinnGen_PheWAS_RG/results/carey_rint_tuned_15clusters_neff_max_empirical_covz/rg/cluster_trait_rg_long_with_p.tsv` | 42,795 rows, 34 cols. Clusters are scored **independently** → safe to subset whole clusters for fast fixtures. |
| Trait×trait rg | `../UKBB_CLUSTER_GWAS/data/finngen_rg/finngen_R12_FIN.ldsc.summary.tsv` | LDSC `--rg` summary; 100% disease-endpoint coverage. |

---

## Conventions

### Compute / infrastructure
- Reproducibility is the headline: **pin everything** (base image by tag, CRAN via dated P3M snapshot,
  GenomicSEM by commit), reference images **by version, never `latest`**.
- Carry the parent project's container fixes: `procps` (so Nextflow `ps`/trace works), `USER root` (GCS-FUSE
  write perms on Google Batch), `ENTRYPOINT []` (don't let an upstream entrypoint intercept `.command.run`).
  See [`../UKBB_CLUSTER_GWAS/docker/postgwas/README.md`](../UKBB_CLUSTER_GWAS/docker/postgwas/README.md).
- Threading: export `OMP_NUM_THREADS=OPENBLAS_NUM_THREADS=MKL_NUM_THREADS=1` (avoid BLAS oversubscription);
  drive parallelism only via `--threads ${task.cpus}` → `setDTthreads` + `future::plan(multicore)`.
- Nextflow output via `output {}` block + `outputDir`, not in-process `publishDir` (parent convention).

### Data schemas (load-bearing — verified from real headers/code)
- **rg long-table (Input A):** required `(cluster_label, trait_id, trait_category, trait_group, rg, rg_se, p,
  h2_trait, h2_trait_se, ldsc_converged, negative_h2, status)`. One row per (cluster, trait). Booleans are
  `TRUE/FALSE` strings; numerics coerce non-numeric→NA; missing = empty/`NA`. **Ignore any precomputed `z`/`abs_rg`
  columns — the engine recomputes them.** Gate: `status==success`, `rg_se>0`, `h2_trait_se>0`, converged,
  `!negative_h2`, then `h2_z = h2_trait/h2_trait_se > z`.
- **trait×trait rg (Input B):** LDSC `--rg` summary `(p1, p2, rg, …, CONVERGED)`; engine uses only `p1,p2,rg,CONVERGED`.
  Long edge-list, each pair once → symmetrize, clip `rg`∈[−1,1], diag=1, missing pairs NaN.
- **GenomicSEM `.rds` (Input C):** list `$S` (k×k genetic **covariance**, standardize to rg), `$V` (sampling cov of
  vech(S) → `rg_se` via delta method), `$I` intercepts.
- **ontology (Input D):** disease `(trait_category, domain, icd_chapter, kind, anchor_eligible, notes)` joins on
  `trait_category`; anthro/lab `(trait_id, …_class, anchor_eligible)` join on `trait_id`. `native` level = the join key.
- **output `category_anchor_scores.tsv`:** `(cluster_label, level, category, eligible, n, n_eff, n_hit, rho_bar,
  vif, auc_abs, auc_signed, perm_p, vif_z, vif_p, pooled_rg, pooled_rg_ci_lo, pooled_rg_ci_hi, coherence,
  mean_abs_rg, mean_signed_rg, odds_ratio, fisher_p, q, rank)`.
- **output `cluster_anchor_labels.tsv`:** `(cluster_label, auto_label, anchor_shape, anchor_margin, anchor_focus,
  n_sig_domains, top_auc, top_q, top_pooled_rg, top_coherence, profile)`.
- **Units / sign / rounding:** `rg`∈[−1,1] LDSC observed scale, signed (sign must survive pooling — coherence
  depends on it). Output rounding replicates Python (`auc`→4, `pooled_rg`→4, `vif`→2, `vif_z`→3, `n_eff`→2,
  `rho_bar`→3, `odds_ratio`→3); `perm_p/vif_p/fisher_p/q` are full precision. `eligible` printed `True/False`.

### Code style
- **Config-over-CLI:** all params live in a YAML config read via `--config`; reuse the parent `cluster_anchoring`
  configs unchanged (only adjust paths). CLI flags reserved for `--threads` / `--z-vector` (later phases).
- Read real headers before parsing; long-format over wide; comment the code but keep it "smart and slim".
- Match the Python reference function-for-function in Phase 1 (the plan carries the line-referenced mapping).

### Gotchas (these save hours)
- **scipy `fisher_exact` returns the *sample* OR `(a·d)/(b·c)`** — NOT R's conditional-MLE `fisher.test$estimate`.
  Compute OR manually; only borrow `fisher.test(...)$p.value`. The single most likely silent R-port mismatch.
- **`perm_p` is not bit-reproducible across languages** (numpy PCG64 ≠ R RNG). Anchor parity on the deterministic
  `vif_p` + label stability; treat `perm_p`/`q` distributionally (MC tolerance).
- **`poolr::meff` on non-PD matrices:** pre-clean R (NaN→0, diag=1, symmetrize) and clip eigenvalues ≥0 to match
  the Python numpy implementation.
- Subset **whole clusters** for fixtures, never traits within a cluster (changes N → changes AUC).

---

## Running stages (designed; not yet runnable)

```bash
# Engine (once R/ + anchor_map.R exist):
Rscript anchor_map.R --config config/carey_rint15_anthro.yaml

# Regenerate the Python oracle to validate against (in the sibling repo):
cd "../UKBB_CLUSTER_GWAS/scripts/cluster_anchoring"
python3 anchor_categories.py --config configs/carey_rint15_anthro.yaml

# Cross-language parity check (once validation/ exists):
Rscript validation/compare_oracle.R --r-out results/.../category_anchor_scores.tsv \
  --oracle "../UKBB_CLUSTER_GWAS/scripts/cluster_anchoring/output/carey_rint15_anthro/category_anchor_scores.tsv"
```

---

## Validation / testing

No unit-test suite exists yet; the designed validation:
- **Positive control (the gate):** C5_sub0 anthro → "Anthropometric [sharp]" with the exact real-run values:
  `auc_abs=0.9164, n_eff=3.0, rho_bar=0.443, vif=1.89, vif_p=0.03489021688956177, pooled_rg=0.2473 [0.1965,0.2968],
  coherence=1.0, odds_ratio=3.523, q≈0.005497, rank=1`.
- **Negative control:** anthro C3 → `auc_abs=0.0653, q≈0.999, auto_label=ambiguous, shape=weak`; disease-track
  `Quantitative` must never label C5_sub0 (forbidden-FP).
- **Cross-language parity:** every deterministic output column matches the committed Python oracle; labels/shapes/ranks identical; `perm_p` within MC tolerance.
- **Schema/sanity:** output col order/names/rounding exact; `n_eff ≤ n`; `vif ≥ 1`; no NaN explosion.
- **Sensitivity (later):** auto-label stable across z∈{3,4,5}.

---

## Key files to know

- [ANALYSIS_DESIGN.md](ANALYSIS_DESIGN.md) — the ADD (question, schemas, methods, controls, phases). **Source of truth.**
- [.agents/plans/anchormap-phase1-r-engine-port.md](.agents/plans/anchormap-phase1-r-engine-port.md) — the next build, with the full Python→R function mapping and oracle values.
- [`../UKBB_CLUSTER_GWAS/scripts/cluster_anchoring/anchor_categories.py`](../UKBB_CLUSTER_GWAS/scripts/cluster_anchoring/anchor_categories.py) — the reference engine (bit-for-bit spec).
- Parent configs/ontologies/oracle: `../UKBB_CLUSTER_GWAS/scripts/cluster_anchoring/{configs,ontology,output,docs/approach.md}/`.
- Container fixes reference: [`../UKBB_CLUSTER_GWAS/docker/postgwas/README.md`](../UKBB_CLUSTER_GWAS/docker/postgwas/README.md); `docker/genomicsem/Dockerfile` (rocker base + `remotes` pattern).

---

## Working with this repo via Claude

- **READMEs:** use the `create-readme` skill — one top-level index that lists each subproject as a linked list
  item, plus a canonical per-subproject README (Orientation / Inputs / Outputs / Workflow / How to run / Gotchas
  / Results / Related).
- **Finishing a subproject:** use the `report-findings` skill to write a `REPORT.md` (results, findings, and a
  clearly non-binding exploratory interpretation) and link it from the subproject README's Results section.
- **Planning a stage:** `plan-analysis`. **Implementing from a plan:** `execute`. **Orientation:** `prime`. **Committing:** `commit`.

---

## Out-of-scope / external links

- **Sibling project (data source + reference engine):** `../UKBB_CLUSTER_GWAS` (read its `CLAUDE.md` for the
  GWAS pipeline, Google Batch conventions, and the known postGWAS beta-vs-OR bug).
- **Reference data:** FinnGen R12 (FIN LDSC summary); not redistributed here.
- **Methods literature:** Li & Ji (2005) `doi:10.1038/sj.hdy.6800717`; `poolr` (Cinar & Viechtbauer, *J Stat Soft* 101(1)); CAMERA (Wu & Smyth).
