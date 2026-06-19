# AnchorMap

A portable, reproducible **R + Nextflow** tool that generalizes the `cluster_anchoring` method from
the sibling `UKBB_CLUSTER_GWAS` project: given latent **cluster factors** and their genetic
correlations (`rg`) to a trait panel, it scores — competitively, size-aware, and correlation-aware —
**which ontology domain each cluster anchors to**, how confidently, and whether the anchor is *sharp*
or *diffuse*.

**Status: Phases 1–5 complete.** The R engine in [R/](R/) (the `anchor_map` CLI,
[inst/scripts/anchor_map.R](inst/scripts/anchor_map.R)) is a validated drop-in for the Python reference (exact deterministic
parity on the anthro + disease tracks; see [README.md](README.md)). **Phase 2 built**: GenomicSEM
`.rds` ingestion ([R/ingest_rds.R](R/ingest_rds.R) — vech-indexed delta-method `rg_se`, `S`→rg
standardization, factor/panel partition) via `--rds`/`cfg$rds`, plus the `vif_correlation: auto`
redundancy auto-fallback (trait_rg → cluster-profile proxy → `VIF=1`) in
[R/redundancy.R](R/redundancy.R)`:select_corr_source`; validated by the
[tests/testthat/](tests/testthat/) suite (delta-method numeric-diff, round-trip, fallback,
VIF-invariance) against synthetic fixtures, with Phase-1 oracle parity preserved byte-for-byte.
**Phase 3 built**: the parallel h²-reliability **z-sweep** ([R/sensitivity.R](R/sensitivity.R) —
`score_at_z` per-z re-run + `run_sensitivity` over `cfg$z_vector` (default `{3,4,5}`) via
`future.apply` `parallel_lapply`, with `--z-vector`/`--threads`), emitting two extra TSVs
(`sensitivity_z_scores.tsv`, `sensitivity_z_labels.tsv` with a per-cluster `label_stable` flag)
alongside the unchanged primaries; validated by the [tests/testthat/](tests/testthat/) suite
(primary-slice parity incl. `perm_p`, thread-invariance, `label_stable`, gate monotonicity). **New
deps:** `future`, `future.apply`. Determinism is engineered: each z-task re-seeds with `random_seed`
**and pins the RNG kind to Mersenne-Twister** (future.seed flips it to L'Ecuyer), so the
z = `h2_z_threshold` slice is byte-identical to the Phase-1/2 single-z primaries and the whole sweep
is thread-count- and backend-invariant; `perm_p` stays serial in `score.R` to protect that parity.
**Phase 4 built**: the **visualization module** ([R/plot.R](R/plot.R) + the `plot_anchors` CLI
[inst/scripts/plot_anchors.R](inst/scripts/plot_anchors.R), config [inst/configs/example_plots.yaml](inst/configs/example_plots.yaml))
ports the reference `plot_anchors.py`/`plot_specificity*.py` encodings to ggplot2 — per-track lollipop
small-multiples, a cluster×category dot-heatmap, an AUC-vs-coherence diagnostic, and the cross-cluster
specificity heatmap + its diagonal reduction — PNG+PDF, headless (ragg if present else cairo),
config-driven, reading only the scored TSVs. The four channels stay distinct (AUC = x/size,
signed `pooled_rg` = diverging colour, coherence = alpha, `q<q_sig` = ring/mask) so the AUC↔rg
divergence at sign-split classes survives. The only recomputation (cross-cluster specificity z) is
**byte-identical to the Python reference's `cluster_distinctive_categories.tsv`** on the disease track.
**New deps:** `ggplot2`, `patchwork`, `scales`, `ggrepel` (and optional `ragg`). **Phase 5 built**:
the pinned, self-validating **Docker image** ([docker/Dockerfile](docker/Dockerfile)) is now the tool +
the primary run interface (`docker run anchormap:0.1.0 Rscript /opt/anchormap/anchor_map.R --config <yaml>`):
`rocker/r-ver:4.6.0` (= the host R the engine was validated on) + a single dated P3M snapshot reproducing
the validated `future.apply 1.20.2` / `ggplot2 4.0.3`, the `procps` / `USER root` / `ENTRYPOINT []` fixes,
and two build-time self-tests (the synthetic-`.rds` engine run recovers C5_sub0 → anthro [sharp]; the
ggplot stack renders a figure — so a bad dep/engine/figure fails `docker build`). **GenomicSEM is omitted**
(the engine reads `ldsc()` `.rds` with base `readRDS`, never runs `ldsc()`). **Nextflow is a
container-validation harness, NOT an orchestration layer** ([nextflow/main.nf](nextflow/main.nf)): a single
`ANCHORMAP_SMOKE` process proves the image runs flawlessly under Nextflow — `test` (local) is the CI gate
for the entrypoint/procps/output-capture contract; `gcp` (Google Batch, spot) is a one-time check for the
`USER root`/GCS-FUSE write fix. **Two deliberate divergences from the ADD §7.3**: base image `4.6.0` not
`4.4.2` (to match the validated env — the 4.4.2 "parent-parity" reason was GenomicSEM image alignment,
which AnchorMap doesn't use), and Nextflow scoped to validation, not production orchestration. The project
is **under git** (GitHub: `micpreuss/AnchorMap`, private); the vendored `claude-science-scaffold/` subdir is
gitignored (it is its own repo). Read [ANALYSIS_DESIGN.md](ANALYSIS_DESIGN.md) as the source of truth (note
the two Phase-5 divergences above); the plan for the next phase goes in `.agents/plans/`.

---

## Project type

Five-axis classification:

- **Orchestration:** a local R engine (the tool) **+** a Nextflow DSL2 **container-validation harness**
  (a single `ANCHORMAP_SMOKE` process — built). Nextflow is *not* a production run path; AnchorMap is run
  directly via `docker run` / `Rscript`. The harness only proves the image obeys the Nextflow container contract.
- **Compute backend:** local R / `docker run` for the engine; the `gcp` (Google Batch, spot) profile exists
  only for the one-time container-on-Batch FUSE check, mirroring `UKBB_CLUSTER_GWAS`.
- **Domain:** Statistical genetics (genetic-correlation / cluster anchoring; Li & Ji n_eff, CAMERA VIF, Mann–Whitney AUC, IVW Fisher-z, BH-FDR) **and** a reusable **R method/tool package** generalizing a Python reference.
- **Data locality:** External — consumes `UKBB_CLUSTER_GWAS` outputs (sibling-repo paths today; `gs://` under the `gcp` profile later). AnchorMap itself ships only small test fixtures.
- **Reproducibility:** Containers **(built)** — pinned `rocker/r-ver:4.6.0` (= the validated host R) + a single
  dated Posit P3M CRAN snapshot (no GenomicSEM; see Phase-5 status). The image self-validates at build time.

The posture: a small R package (`R/*.R` modules + a `--config <yaml>` CLI) validated by **cross-language
parity against the Python reference** (analytic tests in [tests/](tests/)), shipped as a pinned image.

---

## Tech stack (designed)

| Layer | Tools |
|---|---|
| Orchestration | Nextflow DSL2 process `ANCHORMAP`; `output {}` block + `outputDir` (NOT `publishDir`) — inherited from parent |
| Compute | local R (engine); Google Batch spot for production; submit from `nf-head`, not laptop |
| Containers / envs | `rocker/r-ver:4.6.0` (pinned; = the validated host R — see the Phase-5 divergence note above) + `procps` + `USER root` + `ENTRYPOINT []`; image `anchormap:0.1.0` → Artifact Registry `us-central1-docker.pkg.dev/lencz-lab-cogent-1/docker-images/`; **referenced by version tag, never `latest`** |
| Storage | reads sibling `UKBB_CLUSTER_GWAS` files; writes `results/<run_label>/{primary,sensitivity,figures,logs}/` |
| Languages | R ≥4.4 (`data.table`, `future`/`future.apply`, `yaml`, `optparse` (CLI), `testthat`, optional `poolr`; plotting `ggplot2`/`patchwork`/`scales`/`ggrepel`/`ragg`) |
| Methods | Li & Ji (2005) n_eff via `poolr::meff(R,"liji")`; CAMERA VIF; Mann–Whitney/Wilcoxon AUC; label-permutation null; IVW Fisher-z pooled rg; BH-FDR; anchor-shape ruleset |
| External reference data | FinnGen **R12** FIN LDSC `--rg` summary (trait×trait rg); curated `category/anthro/lab` ontologies; GenomicSEM `ldsc()` S/V objects |

---

## Repository structure

Installable R package (`anchormap`); standard package layout (`R/` source, `man/` roxygen docs,
`inst/` installed assets, `tests/testthat/`):

```
AnchorMap/
├── DESCRIPTION / NAMESPACE / LICENSE        ← R package metadata (Imports incl. optparse; CLI parser)
├── CLAUDE.md                 ← this file (start here)
├── ANALYSIS_DESIGN.md        ← ★ source of truth: the 15-section ADD
├── README.md / README.production.md   ← public front page / detailed dev+provenance notes
├── .agents/plans/            ← per-phase build plans (phase1…phase5)
├── .claude/                  ← injected scaffold (skills, settings); do not treat as project code
├── claude-science-scaffold/  ← vendored source of .claude (own git repo, gitignored) — not part of AnchorMap
├── R/{io,gate,redundancy,score,label}.R   ← engine modules (Phase 1)
│   ├── ingest_rds.R          ← Phase-2 GenomicSEM .rds reader (vech delta-method, partition, standardize)
│   ├── sensitivity.R         ← Phase-3 parallel z-sweep (score_at_z, run_sensitivity, future.apply parallel_lapply)
│   ├── plot.R                ← Phase-4 ggplot2 figures (lollipop, dot-heatmap, AUC-vs-coherence, specificity + diagonal)
│   ├── run_anchormap.R       ← library entry points: run_anchormap() / run_plots() (config resolve, orchestration)
│   ├── cli.R                 ← optparse CLI: parse_*/cli_* for the anchor_map + plot_anchors entry points
│   └── anchormap-package.R   ← package-level roxygen / imports
├── man/                      ← roxygen2-generated .Rd help (exported run_*/read_*/score_at_z)
├── inst/configs/             ← shipped example configs (example_{anthro,disease,plots}.yaml + synthetic_rds{,_plots}.yaml)
├── inst/fixtures/            ← self-contained synthetic .rds + ontology (make_synthetic_ldsc.R builds them)
├── inst/ontology/            ← disease/anthro/lab ontology TSVs (Input D)
├── inst/scripts/             ← Rscript CLI shims: anchor_map.R, plot_anchors.R (call anchormap:::cli_*)
├── local/configs/            ← machine-specific Carey/FinnGen configs (absolute paths; gitignored)
├── tests/testthat/           ← testthat suite (analytic, ingest, sensitivity, cli) + helper-fixtures.R
├── validation/               ← R-vs-Python oracle comparator (compare_oracle.R, run_oracle.sh)
├── results/<run_label>/      ← engine outputs (TSVs + anchormap.log) + figures/ — generated, gitignored
├── docker/                   ← Phase-5 THE TOOL: Dockerfile (rocker/r-ver:4.6.0 + P3M pin + 2 self-tests) + bin/ shims + README
└── nextflow/                 ← Phase-5 container-validation harness: main.nf (ANCHORMAP_SMOKE) + nextflow.config + params/{test,gcp}.yaml
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
                              scored TSVs ─► plot.R (Phase 4) ─► figures/ (lollipop · dot-heatmap · AUC-vs-coherence · specificity + diagonal)
```

- **Two input routes (same engine):** the standardized `rg` long-table TSV **or** a GenomicSEM `ldsc()` `.rds`
  from which both the long-table (rg = S/√(diag·diag); rg_se via **delta-method on V**) and the trait×trait
  matrix are derived. *Both routes built: TSV (Phase 1) + the `.rds` route (`--rds`/`cfg$rds`, Phase 2,
  [R/ingest_rds.R](R/ingest_rds.R)).*
- **Redundancy source:** explicit `vif_correlation` modes (`trait_rg`/`cluster_profile`) honored verbatim
  (Phase-1 parity); `vif_correlation: auto` (Phase 2, `select_corr_source`) self-selects trait×trait rg when
  coverage ≥ `vif_coverage_min` (0.5), else the cluster-profile proxy (≥3 clusters), else VIF=1 with a loud WARN.
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

### Git / version control
- **Do NOT create new branches.** Work is already on a per-phase feature branch — commit (and merge)
  on the current branch; never branch again, even for a new phase. Confirm the branch first if unsure.
- Commits use Conventional Commits tags and carry **no** `Co-Authored-By` footer (match existing history).

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

## Running stages

**The tool is now an installable R package (`anchormap`).** Entry points: `run_anchormap()` /
`run_plots()` / `run_sensitivity()` (library) and the `anchor_map` / `plot_anchors` CLI shims (or
`Rscript inst/scripts/*.R`). There is no top-level `anchor_map.R` any more; `--config` accepts a YAML
path **or** a bare shipped-config name (resolved from `inst/configs/`). Machine-specific Carey/FinnGen
configs live in the gitignored `local/configs/`.

```bash
# Install the package (host):
R CMD INSTALL .            # or: remotes::install_github("micpreuss/AnchorMap")

# Engine (host, after install) — bare config name resolves to inst/configs, or pass a YAML path:
anchor_map --config synthetic_rds --out-dir results/synthetic_rds --threads 4
#  (without the PATH shim:  Rscript inst/scripts/anchor_map.R --config synthetic_rds --out-dir ... )
#  your real runs:  anchor_map --config local/configs/carey_rint15_anthro.yaml --out-dir results/...

# Engine — via the pinned image (THE primary, reproducible run interface; mount cwd as /work):
docker run --rm -v "$PWD:/work" -w /work anchormap:0.1.0 \
  anchor_map --config synthetic_rds --out-dir results/synthetic_rds --threads 4

# Figures (reads the scored TSVs the engine wrote; PNG+PDF into the --out-dir):
plot_anchors --config synthetic_rds_plots --out-dir results/synthetic_rds/figures

# Tests (testthat; load_all harness so internals are visible):
Rscript -e 'testthat::test_local()'

# Build the image (build IS the self-test: C5_sub0 anthro sharp + a figure render):
docker build -t anchormap:0.1.0 -f docker/Dockerfile .          # release: add --platform linux/amd64

# Nextflow container-validation harness (NOT how you run AnchorMap):
nextflow run nextflow/main.nf -profile test -params-file nextflow/params/test.yaml   # local CI gate
# nextflow run nextflow/main.nf -profile gcp -params-file nextflow/params/gcp.yaml   # one-time Batch/FUSE check (needs push + GCP creds)

# Regenerate the Python oracle to validate against (in the sibling repo):
cd "../UKBB_CLUSTER_GWAS/scripts/cluster_anchoring"
python3 anchor_categories.py --config configs/carey_rint15_anthro.yaml

# Cross-language parity check (once validation/ exists):
Rscript validation/compare_oracle.R --r-out results/.../category_anchor_scores.tsv \
  --oracle "../UKBB_CLUSTER_GWAS/scripts/cluster_anchoring/output/carey_rint15_anthro/category_anchor_scores.tsv"
```

---

## Validation / testing

The `testthat` suite lives in [tests/testthat/](tests/testthat/) (run `Rscript -e 'testthat::test_local()'`);
it ports the former analytic/Phase-2/Phase-3 scripts. The validation:
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
- [README.production.md](README.production.md) — the detailed dev/provenance README (phase history, validation, full schemas). The visible [README.md](README.md) is the de-branded public front page; **keep dev notes in the production copy.**
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
