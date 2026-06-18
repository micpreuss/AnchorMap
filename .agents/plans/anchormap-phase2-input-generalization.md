# Analysis: AnchorMap Phase 2 — Input generalization (GenomicSEM `.rds` ingestion + auto redundancy fallback)

Validate the schema contracts and method usage against the actual code/data before implementing.
Pay special attention to: the GenomicSEM `V`-matrix **vech ordering** (column-major lower-triangle),
the **delta-method** gradient of `cov2cor`, the **sign** of `rg` surviving standardization, and
**preserving Phase-1 parity** for the existing configs (the `auto` path is opt-in).

> Scope: **Phase 2 of `ANALYSIS_DESIGN.md`** only — Aims **A2** (standardized GenomicSEM/LDSC ingestion)
> and **A3** (default trait×trait → proxy → VIF=1 auto-fallback). Build on the Phase-1 R engine
> (`R/{io,gate,redundancy,score,label}.R` + `anchor_map.R`), already validated bit-for-bit against the
> Python oracle. **Out of Phase 2 (do NOT build here):** the parallel z-sweep + multi-CPU `perm_p`
> (Phase 3), the plotting module (Phase 4), the Docker image + Nextflow process (Phase 5). **Invariant to preserve:** every existing
> config (`carey_rint15.yaml`, `carey_rint15_anthro.yaml`) and the oracle parity must keep passing
> unchanged — Phase 2 only *adds* the `.rds` route and the `vif_correlation: auto` mode.

## Question & object
Make the AnchorMap engine consume a **single standard GenomicSEM `ldsc()` artifact** (`$S`, `$V`, `$I`
in a `.rds`) as an alternative to the cluster×trait long-TSV + trait×trait LDSC summary, and make the
within-category redundancy source **self-selecting** (trait×trait rg by default, automatic fallback to
the cluster-profile proxy, then to no-deflation `VIF=1`, all logged). Object in → a GenomicSEM `.rds`
(named list of genetic **covariance** `S`, sampling covariance `V`, intercepts `I`) **or** the Phase-1
TSV pair; object out → the *same* two TSVs as Phase 1 (`category_anchor_scores.tsv`,
`cluster_anchor_labels.tsv`). The inference licensed: "AnchorMap can be driven directly from a
GenomicSEM run with no hand reformatting, and it never silently runs un-deflated."

## Analyst story
As an analyst / I want AnchorMap to ingest a GenomicSEM `ldsc()` `.rds` directly (deriving both the
cluster-factor×trait rg long-table and the trait×trait redundancy matrix from `S`/`V`), and to pick the
redundancy source automatically based on coverage / so that I can anchor clusters straight from a
GenomicSEM artifact across tracks (disease, lab, anthro) without per-track manual flag-tuning, and never
get an anti-conservative VIF without a loud warning.

## Pipeline position
- **Upstream:**
  - **New (Input C):** GenomicSEM gPCA stage → `<cluster>.ldsc_output.rds` (`$S`, `$V`, `$I`); produced
    by `UKBB_CLUSTER_GWAS/scripts/genomic_sem/scripts/run_cluster_gpca.R` (`saveRDS(ldsc_output, …)` at
    L495). **NB: no real `.rds` exists in the parent repo today and GenomicSEM is not installed locally**
    (see Validation — the gate rests on analytic + numeric-diff + round-trip, per the chosen strategy).
  - **Existing (Inputs A/B):** the rg long-TSV + `finngen_R12_FIN.ldsc.summary.tsv` (Phase-1 route, unchanged).
- **Downstream:** the unchanged Phase-1 scoring chain `gate.R → redundancy.R → score.R → label.R` and the
  two output TSVs; Phase 3 (`sensitivity.R`), Phase 4 (`plot.R`) and Phase 5 (Nextflow `ANCHORMAP`) build on this.
- **Orchestration pattern to mirror:** [anchor_map.R](../../anchor_map.R) (the `run_anchormap()` driver,
  the `[load]/[gate]/[vif]/[write]` log lines, `resolve_path`/`stage_root_of`). Reference for the `.rds`
  standardization: `run_cluster_gpca.R` **L414–422** (`S_Stand`).

---

## CONTEXT REFERENCES — READ BEFORE IMPLEMENTING

### Input schemas (verified from code/headers)

**Input C — GenomicSEM `ldsc()` object (`.rds`; NEW).** `readRDS(path)` → a named list. The parent runs
`ldsc(..., stand = FALSE)` (see `run_cluster_gpca.R` L401–411, deliberately, to avoid an internal crash
on negative h²), so the `.rds` carries the **unstandardized** objects and AnchorMap must standardize
itself:

| element | shape | meaning / use |
|---|---|---|
| `$S` | k×k | genetic **covariance** (observed scale). `diag(S)` = genetic variances (≈ h² on observed scale). Row/col `dimnames` are the variable names (cluster factors **and** panel traits, joint run). |
| `$V` | q×q, `q = k(k+1)/2` | sampling covariance of **`vech(S)`** in **column-major lower-triangle order** (the order of `S[lower.tri(S, diag = TRUE)]`): `(1,1),(2,1),…,(k,1),(2,2),(3,2),…,(k,k)`. **Load-bearing — index it wrong and every SE is wrong.** |
| `$I` | k×k | LDSC intercepts (carried as provenance; not consumed by scoring). |

- **Partition of `S` (chosen contract — name regex + explicit override).** Split `dimnames(S)` into
  **cluster factors** (→ long-table rows) and **panel traits** (→ trait×trait redundancy matrix). Config:
  - `cluster_factor_pattern` (regex, **default `"^C[0-9]"`** — matches `C0`, `C5`, `C5_sub0`); names matching → factors.
  - `cluster_factors:` (optional explicit character list) **overrides** the regex when present.
  - Everything in `dimnames(S)` not a factor is a **panel trait**. Error if either set is empty.
- **`trait_category` / `trait_group` for the long-table.** The `.rds` carries only names, not ontology
  columns. Rule: `trait_group = cfg$trait_group` for every panel trait; `trait_category` is filled from an
  **optional** `rds_trait_meta` TSV (`trait_id, trait_category[, trait_group]`) when provided; if
  `ontology_key == "trait_id"` (anthro/lab tracks) `trait_category` is unused and the map is unnecessary;
  if `ontology_key == "trait_category"` and no map is given → **error** with a clear message (can't join
  the disease ontology without it).

**Output of the `.rds` reader → the Phase-1 long-table contract** (so the rest of the engine is route-agnostic).
For each (factor `f`, panel trait `t`) build one row with exactly the `.LONG_REQUIRED` columns
([io.R:45](../../R/io.R#L45)):

| column | derivation from `.rds` |
|---|---|
| `cluster_label` | factor name `f` |
| `trait_id` | panel trait name `t` |
| `trait_category` | from `rds_trait_meta` (or `NA` when `ontology_key=="trait_id"`) |
| `trait_group` | `cfg$trait_group` |
| `rg` | `S_Stand[f, t]` (standardized; see below) |
| `rg_se` | `sqrt(V_Stand_pair(f, t))` (delta method; see below) |
| `p` | `2 * pnorm(-abs(rg / rg_se))` (two-sided z; used only by the ORA layer) |
| `h2_trait` | `S[t, t]` (trait genetic variance / observed-scale h²) |
| `h2_trait_se` | `sqrt(V[idx(t,t), idx(t,t)])` (the `vech` diagonal entry for `(t,t)`) |
| `ldsc_converged` | `TRUE` if `rg, rg_se, h2_trait, h2_trait_se` all finite, else `FALSE` |
| `negative_h2` | `S[t, t] < 0` |
| `status` | `"success"` if `ldsc_converged`, else `"failed"` |

The trait×trait redundancy matrix is the **panel-trait block** of `S_Stand` (`S_Stand[panel, panel]`),
clipped to [−1,1], diag = 1 — same object `build_trait_rg_matrix` would have produced, so it drops into
the existing `reindex_corr`/`rho_bar`/`meff_liji` path unchanged.

**Standardization (mirror `run_cluster_gpca.R` L414–422):**
```r
s_diag         <- diag(S)
s_diag_clamped <- pmax(s_diag, 0)                 # clamp negative h² so sqrt is real
denom          <- outer(sqrt(s_diag_clamped), sqrt(s_diag_clamped))
denom[denom == 0] <- 1                            # avoid 0/0 for truly-zero h² traits
S_Stand        <- S / denom
diag(S_Stand)  <- 1
```
Do **not** clip `S_Stand` off-diagonals for the long-table `rg` (Phase-1 gate clips ±0.999 inside `y`/`v`
and the matrix builder clips ±1 for the redundancy matrix — keep those downstream behaviours).

**Delta-method `rg_se` (the highest-risk step).** For `r_ij = S_ij / sqrt(S_ii·S_jj)`, `r_ij` depends only
on `S_ii, S_jj, S_ij`, so its variance is exact (not approximate) from the 3×3 `V`-submatrix:
```
g = ( ∂r/∂S_ij , ∂r/∂S_ii , ∂r/∂S_jj )
  = ( 1/sqrt(S_ii·S_jj) , −r_ij/(2·S_ii) , −r_ij/(2·S_jj) )
idx3 = vech positions of (i,i),(j,j),(i,j)
Var(r_ij) = gᵀ · V[idx3, idx3] · g      ;   rg_se_ij = sqrt(Var)
```
Guard: if `S_ii ≤ 0` or `S_jj ≤ 0` → `rg_se = NA`, `ldsc_converged = FALSE` (the gate then drops the row).
The `vech` index helper: for a k×k matrix, position of `(i,j)` with `i ≥ j` is
`sum(k:(k-j+2)) + (i-j+1)` (column-major lower triangle); implement + unit-test it directly.

**Input A/B (unchanged, Phase-1):** rg long-TSV (`read_long`, [io.R:49](../../R/io.R#L49)); trait×trait
LDSC summary (`build_trait_rg_matrix`, [redundancy.R:12](../../R/redundancy.R#L12)). Schemas already verified
in the Phase-1 plan — not repeated here.

**Input D — ontology (unchanged):** disease keyed on `trait_category`, anthro/lab on `trait_id`
(`read_ontology`, [io.R:66](../../R/io.R#L66)).

### Output schema (contract — unchanged from Phase 1)
`category_anchor_scores.tsv` and `cluster_anchor_labels.tsv` exactly as `.SCORE_COLS` /`.LABEL_COLS`
([anchor_map.R:39-44](../../anchor_map.R#L39-L44)). Phase 2 adds **no output columns** — it adds an input
route and a redundancy-source decision, both recorded in `anchormap.log`.

### Reference data & methods
- **GenomicSEM `ldsc()`** — produces `S`/`V`/`I`; `V` is `vech(S)` sampling covariance, column-major lower
  triangle. Standardization reference: `run_cluster_gpca.R` L414–422. Delta method for `cov2cor` =
  what GenomicSEM does internally when `stand=TRUE` (here recomputed, since the parent ran `stand=FALSE`).
- **Coverage / fallback** — current manual logic: [anchor_map.R:73-84](../../anchor_map.R#L73-L84) (trait_rg
  branch computes `cov` and warns at <0.5) and Python `main()` L502–514. Proxy builder
  `build_trait_profile_corr` ([redundancy.R:37](../../R/redundancy.R#L37)) needs ≥3 clusters/trait
  (`min_periods=3`). ADD §6 "Trait×trait → proxy auto-fallback" is the spec.
- **Methods unchanged:** Li & Ji `n_eff` (`meff_liji`), CAMERA VIF, etc. — Phase-1 functions consume the
  `.rds`-derived objects with no change.
- **Assumptions/failure modes:** `S` may be non-PD / have negative diagonals (clamp for the denom; rows with
  negative h² get `negative_h2=TRUE` and are gate-dropped). Low trait×trait coverage → fallback. `V`
  mis-indexing is the dominant failure mode → numeric-diff test is mandatory.

### Files to read / create
- READ: [R/io.R](../../R/io.R) — extend with the `.rds` readers; reuse `resolve_path`/`stage_root_of`/`default_config`.
- READ: [R/redundancy.R](../../R/redundancy.R) — the proxy + matrix builders the fallback selects between.
- READ: [anchor_map.R](../../anchor_map.R) — the driver; refactor the `[vif]` block (L73–84) to call the new selector.
- READ: `../UKBB_CLUSTER_GWAS/scripts/genomic_sem/scripts/run_cluster_gpca.R` L395–500 — `S_Stand` + `saveRDS` shape.
- READ: [tests/run_tests.R](../../tests/run_tests.R) — current test harness (no `testthat` locally; mirror its plain-`stopifnot` style).
- CREATE: `R/ingest_rds.R` — the GenomicSEM `.rds` reader/standardizer/delta-method/partition (new module; keeps `io.R` slim).
- CREATE: `tests/fixtures/synthetic_ldsc.rds` (+ generator script) — a hand-built `S`/`V` with known ground truth.
- CREATE: `tests/test_phase2.R` — delta-method numeric-diff, partition, round-trip, fallback unit tests.
- CREATE (optional): `config/carey_rint15_rds.yaml` — example `.rds`-route config (only if a real/synthetic `.rds` is wired).

---

## METHOD / IMPLEMENTATION PLAN

### Phase 1: `.rds` ingestion (`R/ingest_rds.R`)
1. `read_ldsc_rds(path)` — `readRDS`; assert it's a list with finite-dimensioned `$S` (square, named) and
   `$V` (square, `nrow == k(k+1)/2`); warn-and-carry `$I`. Clear errors on shape mismatch.
2. `vech_index(k)` — return the `(i,j) → position` map for column-major lower triangle; unit-tested.
3. `standardize_S(S)` — the L414–422 block; returns `S_Stand` (+ keeps `s_diag` for h² extraction).
4. `rg_se_matrix(S, V)` — for every `(i,j)`, `i>j`, the 3×3 delta-method variance → a k×k `rg_se` matrix
   (NA where `S_ii≤0` or `S_jj≤0`).
5. `partition_S(names, cfg)` — factors via `cluster_factors` list else `cluster_factor_pattern` regex;
   panel = the rest; error if either empty.
6. `rds_to_long(S, S_Stand, rg_se, factors, panel, cfg, trait_meta)` — assemble the Phase-1 long-table
   contract (table above), one row per factor×panel pair.
7. `rds_to_trait_rg(S_Stand, panel)` — `S_Stand[panel, panel]`, clip ±1, diag 1 (drop-in for `build_trait_rg_matrix`).

### Phase 2: redundancy-source auto-selection (`R/redundancy.R`)
8. `trait_rg_coverage(corr, tids)` — factor out the existing coverage calc
   ([anchor_map.R:77-78](../../anchor_map.R#L77-L78)): reindex, NA the diagonal, fraction of traits with ≥1
   finite off-diagonal.
9. `identity_corr(tids)` — diag-1 / NA-off-diagonal matrix (→ `rho_bar = 0` → `VIF = 1`; `meff_liji` → `n_eff = n`).
10. `select_corr_source(g, cfg, sroot, trait_rg_override = NULL, emit)` — returns
    `list(corr, source, coverage, reason)`. Honors `cfg$vif_correlation`:
    - `"trait_rg"` (explicit) → build (or use override) trait_rg; compute coverage; warn at `< vif_coverage_min`. **Unchanged Phase-1 behaviour → parity preserved.**
    - `"cluster_profile"` (explicit) → proxy. **Unchanged.**
    - `"auto"` (**new**) →
      ```
      cov = coverage(trait_rg) if trait_rg available else 0
      if trait_rg available AND cov >= vif_coverage_min:   source = trait_rg
      elif n_clusters >= 3:                                source = cluster_profile (proxy)
      else:                                                source = identity (VIF=1) + loud WARN
      ```
    On the `.rds` route the trait_rg matrix is the `rds_to_trait_rg` block (passed as `trait_rg_override`),
    so coverage is typically 100%.

### Phase 3: driver integration (`anchor_map.R`)
11. Add `--rds <path>` (and/or `cfg$rds`) as an input route: if set, `read_ldsc_rds` → standardize → derive
    `df` (long-table) + `trait_rg_override`; else the Phase-1 `read_long` route. Both converge on the same
    `g <- apply_universe_gate(df, cfg)`.
12. Replace the inline `[vif]` block (L73–84) with `select_corr_source(...)`; log `source`, `coverage`,
    and the fallback `reason`. Everything from `score_cluster_level` onward is untouched.
13. Add Phase-2 keys to `default_config()` ([io.R:12](../../R/io.R#L12)): `vif_coverage_min = 0.5`,
    `cluster_factor_pattern = "^C[0-9]"`, `cluster_factors = NULL`, `rds = NULL`, `rds_trait_meta = NULL`.
    Keep all Phase-1 defaults; **do not change `vif_correlation`'s default** (`cluster_profile`).

### Phase 4: Validation (analytic + numeric-diff + round-trip + fallback)
See Validation strategy below.

---

## STEP-BY-STEP TASKS (execute top to bottom; each atomic + checkable)

### CREATE `R/ingest_rds.R` — vech index + standardize + delta-method
- **IMPLEMENT**: `vech_index(k)`, `standardize_S(S)`, `rg_se_matrix(S, V)`.
- **PATTERN**: `standardize_S` mirrors `run_cluster_gpca.R` L414–422 verbatim.
- **DATA/SCHEMA**: in = `$S` (k×k cov), `$V` (q×q, q=k(k+1)/2); out = `S_Stand` (k×k), `rg_se` (k×k, NA on bad diag).
- **GOTCHA**: `V` is **column-major lower-triangle `vech`** — verify the index map before trusting any SE.
  Clamp negative `diag(S)` only for the denom; never overwrite `S` itself (h² extraction needs the raw diag).
- **VALIDATE**: `Rscript -e 'source("R/ingest_rds.R"); stopifnot(vech_index(3)[["2,1"]]==2, vech_index(3)[["3,3"]]==6)'`

### CREATE `R/ingest_rds.R` — partition + long-table + trait_rg derivation
- **IMPLEMENT**: `read_ldsc_rds`, `partition_S`, `rds_to_long`, `rds_to_trait_rg`.
- **DATA/SCHEMA**: `rds_to_long` emits exactly `.LONG_REQUIRED` ([io.R:45](../../R/io.R#L45)); `rds_to_trait_rg`
  emits a trait×trait matrix matching `build_trait_rg_matrix`'s shape.
- **GOTCHA**: `trait_category` requires `rds_trait_meta` when `ontology_key=="trait_category"` — error clearly
  if absent. `negative_h2 = S[t,t] < 0`; `status="failed"` rows must survive into `df` so the gate (not the
  reader) drops them — matches the TSV route.
- **VALIDATE**: on the synthetic fixture, `df <- rds_to_long(...)`; `stopifnot(all(.LONG_REQUIRED %in% names(df)), is.numeric(df$rg), is.numeric(df$rg_se))`.

### CREATE `tests/fixtures/` synthetic `.rds` generator
- **IMPLEMENT**: build a small PD covariance `S` (e.g. 5×5: 2 factors `C0`,`C5_sub0` + 3 traits
  `BMI`,`WT`,`HT`) with known `S_Stand`; set `V` from a known per-element SE (diagonal `V` suffices for a
  first ground-truth, plus one off-diagonal block to exercise the 3×3 path); `saveRDS`.
- **GOTCHA**: name the factor rows to match `cluster_factor_pattern` (`^C[0-9]`); make one trait diag slightly
  negative in a variant fixture to exercise the `negative_h2`/`NA rg_se` guard.
- **VALIDATE**: `Rscript tests/fixtures/make_synthetic_ldsc.R && Rscript -e 'x<-readRDS("tests/fixtures/synthetic_ldsc.rds"); stopifnot(nrow(x$V)==15, all(dim(x$S)==5))'`

### CREATE `tests/test_phase2.R` — delta-method numeric-difference test (the key gate)
- **IMPLEMENT**: pick a pair `(i,j)`; compute analytic `rg_se` via `rg_se_matrix`; recompute by **finite
  differencing**: perturb `S_ii,S_jj,S_ij` by ε, recompute `r_ij`, build the numerical gradient `g_num`,
  then `sqrt(g_numᵀ V[idx3,idx3] g_num)`; assert `abs(analytic − numeric) < 1e-6`.
- **PATTERN**: plain `stopifnot` (mirror [tests/run_tests.R](../../tests/run_tests.R); no `testthat` locally).
- **GOTCHA**: this validates the Jacobian **and** the `vech` indexing together — it is the substitute for the
  unavailable real-`ldsc()` cross-check.
- **VALIDATE**: `Rscript tests/test_phase2.R` → prints `delta-method OK`.

### CREATE `tests/test_phase2.R` — round-trip equivalence (.rds route == TSV route)
- **IMPLEMENT**: from the synthetic `.rds`, derive `df` + `trait_rg`; write them to temp TSVs in the
  Input-A/B format; run the **engine** both ways (`.rds`-derived objects vs. re-read TSVs) through
  `apply_universe_gate → score_cluster_level → rank_and_label`; assert identical score rows.
- **GOTCHA**: clusters are scored independently — the two routes must agree to full deterministic precision
  (perm_p excepted: same `set.seed` → identical here since same R RNG).
- **VALIDATE**: `Rscript tests/test_phase2.R` → `round-trip identical`.

### CREATE `tests/test_phase2.R` — fallback selector unit tests (3 branches)
- **IMPLEMENT**: with `vif_correlation:auto`: (a) full-coverage trait_rg → `source=="trait_rg"`; (b) zero-coverage
  trait_rg + ≥3 clusters → `source=="cluster_profile"`; (c) zero-coverage + <3 clusters → `source=="identity"`,
  and assert resulting `VIF==1`, `rho_bar==0`, `n_eff==n` on a scored category.
- **GOTCHA**: assert the **invariant** — `auc_abs`, `auc_signed`, `pooled_rg`, `coherence`, `rank` are
  byte-identical across all three corr sources for a given category (VIF touches only `vif_p`/CI width).
- **VALIDATE**: `Rscript tests/test_phase2.R` → `fallback branches OK; VIF-invariance OK`.

### UPDATE `R/redundancy.R` — coverage helper + identity corr + `select_corr_source`
- **IMPLEMENT**: tasks 8–10 above.
- **PATTERN**: lift the coverage math from [anchor_map.R:77-78](../../anchor_map.R#L77-L78).
- **VALIDATE**: covered by the fallback unit tests.

### UPDATE `R/io.R` `default_config()` — Phase-2 keys
- **IMPLEMENT**: task 13. **GOTCHA**: additive only; Phase-1 defaults and `vif_correlation="cluster_profile"`
  default unchanged so existing configs are byte-identical.
- **VALIDATE**: `Rscript -e 'source("R/io.R"); d<-default_config(); stopifnot(d$vif_coverage_min==0.5, d$cluster_factor_pattern=="^C[0-9]")'`

### UPDATE `anchor_map.R` — `.rds` route + selector wiring + log lines
- **IMPLEMENT**: tasks 11–12; source `R/ingest_rds.R`; add `--rds`; emit `[ingest]`/`[vif source=… coverage=…%]`/fallback `reason`.
- **GOTCHA**: when `cfg$rds`/`--rds` set, skip `read_long`; else unchanged. Keep the `FINISHED` line + manifest.
- **VALIDATE**: re-run the **Phase-1 configs** (no `.rds`, explicit `vif_correlation`) and confirm the two TSVs
  are **unchanged** (parity preserved); then run the `.rds` example if wired.

---

## VALIDATION STRATEGY
No real `ldsc_output.rds` and no GenomicSEM locally → the gate rests on (chosen strategy):
- **Analytic:** `S_Stand` matches a hand-computed `cov2cor` on the synthetic fixture; `vech_index` matches
  the known column-major lower-triangle positions.
- **Numeric-difference (delta method):** analytic `rg_se` == finite-difference `rg_se` to `1e-6` — the
  substitute for the real-`ldsc(stand=TRUE)` `V_Stand` cross-check (deferred to Phase 5/container or to
  whenever a real `.rds` is supplied; documented as such).
- **Round-trip:** `.rds` route produces identical engine scores to the equivalent TSV route.
- **Fallback:** three branches selected correctly; `VIF=1` branch fires + WARNs when coverage low AND <3
  clusters; **VIF-invariance** of AUC/ranks/pooled_rg/coherence across all corr sources.
- **Parity regression (must not break):** existing `carey_rint15{,_anthro}.yaml` outputs unchanged
  (explicit `vif_correlation` modes untouched) — re-run the Phase-1 oracle comparison.
- **Schema/sanity:** `.rds`-derived `df` has all `.LONG_REQUIRED` cols; `n_eff ≤ n`; `vif ≥ 1`; trait×trait
  block symmetric, unit diagonal; no NaN explosion; negative-h² rows gate-dropped.

## VALIDATION COMMANDS (run all; zero deterministic/control failures)
```bash
# deps present locally (poolr/data.table/yaml/Matrix yes; GenomicSEM/testthat NO — plain stopifnot tests)
Rscript -e 'stopifnot(all(c("poolr","data.table","yaml","Matrix") %in% rownames(installed.packages())))'

# build the synthetic fixture, then the Phase-2 test battery
Rscript tests/fixtures/make_synthetic_ldsc.R
Rscript tests/test_phase2.R          # delta-method, round-trip, fallback, VIF-invariance

# PARITY REGRESSION — existing TSV-route outputs must be unchanged
Rscript anchor_map.R --config config/carey_rint15_anthro.yaml
Rscript anchor_map.R --config config/carey_rint15.yaml
Rscript validation/compare_oracle.R \
  --r-out  results/carey_rint15_anthro/category_anchor_scores.tsv \
  --oracle "../UKBB_CLUSTER_GWAS/scripts/cluster_anchoring/output/carey_rint15_anthro/category_anchor_scores.tsv"

# (optional) .rds route once a real/synthetic .rds config is wired
# Rscript anchor_map.R --config config/carey_rint15_rds.yaml
```

## ACCEPTANCE CRITERIA
- [ ] `.rds` reader asserts `S`/`V` shapes; `vech_index` unit-tested; partition splits factors/panel per `cluster_factor_pattern` (+ explicit override).
- [ ] Delta-method `rg_se` == numeric-difference to `1e-6`; negative-h² traits → `rg_se=NA` and gate-dropped.
- [ ] Round-trip: `.rds` route engine scores == equivalent TSV route (deterministic cols exact).
- [ ] `vif_correlation: auto` selects trait_rg ≥ `vif_coverage_min`, else proxy (≥3 clusters), else `VIF=1`+WARN — all logged with coverage %.
- [ ] **VIF-invariance** holds: AUC/ranks/pooled_rg/coherence identical across corr sources.
- [ ] **Parity regression passes:** existing configs' two TSVs unchanged; oracle comparison clean (explicit modes untouched).
- [ ] Provenance: log records input route (`.rds` vs TSV), redundancy `source`, coverage %, fallback `reason`, config + seed.

## NOTES
- **Backward-compat is a hard constraint.** `auto` is opt-in; the default stays `cluster_profile`; explicit
  `trait_rg`/`cluster_profile` keep Phase-1 behaviour byte-for-byte → the oracle parity gate is untouched.
- **Why no real-`ldsc()` parity now:** no `.rds` in the parent repo + GenomicSEM not installed locally. The
  numeric-difference test validates the Jacobian **and** the `vech` indexing simultaneously, which is the
  thing most likely to be wrong; the real-`ldsc(stand=TRUE)` `V_Stand` cross-check is deferred to the Phase-5
  container (where GenomicSEM is pinned) or to whenever a real `.rds` is provided — and is cheap to add then.
- **`S` partition contract (chosen):** name regex `cluster_factor_pattern` (default `^C[0-9]`) with an explicit
  `cluster_factors:` override. If a future `.rds` uses different factor naming, set the pattern/list in config —
  no code change.
- **Single-factor `.rds`.** If a `.rds` carries one cluster factor, the proxy fallback (`build_trait_profile_corr`,
  needs ≥3 clusters) is unavailable for that run; with a full-coverage trait_rg block from `S_Stand` this is moot,
  but the `auto` logic still correctly lands on trait_rg (or identity+WARN if somehow uncovered).
- **VIF-affects-only invariant** (ADD §6) is asserted, not assumed — it's the safety rail that lets the
  redundancy source switch automatically without perturbing the headline AUC/label.
```
