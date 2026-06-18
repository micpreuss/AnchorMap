# Analysis: AnchorMap Phase 1 — R engine port + parity fixture

Validate the schema contracts and method usage against the actual code/data before implementing.
Pay special attention to column names, units, sign, rounding, RNG, and the two cross-language gotchas
(scipy sample-OR vs R conditional-MLE; perm_p RNG is not reproducible across languages).

> Scope: **Phase 1 of `ANALYSIS_DESIGN.md`** only. Build a faithful R reimplementation of the reference
> Python engine `anchor_categories.py`, reading the **standardized long-TSV** input on a **single z**
> (read from config), and prove it reproduces the Python output. **Out of Phase 1 (do NOT build here):**
> GenomicSEM `.rds` ingestion + delta-method rg_se (Phase 2), automatic trait×trait→proxy coverage
> fallback (Phase 2), the parallel z-sweep + multi-CPU perm_p (Phase 3), the plotting module (Phase 4),
> the Docker image + Nextflow process (Phase 5). Phase 1 *honors the config's `vif_correlation` flag exactly as Python does.*

## Question & object
Given a cluster×trait genetic-correlation **long table** (one row per cluster × trait), reproduce — in R —
the competitive, correlation-aware anchoring score per (cluster, ontology level, category): Mann–Whitney
**AUC**, label-permutation **perm_p**, CAMERA **VIF-corrected z**, BH-FDR **q**, IVW Fisher-z **pooled_rg**
+ coherence, Fisher **ORA**, and the per-cluster **auto_label** + **anchor_shape**. Object in → long TSV
(+ trait×trait LDSC rg summary + ontology TSV + YAML config). Object out → two TSVs matching the Python
contract. The inference licensed: "this R engine is a trustworthy drop-in for the Python reference," which
is the precondition for every later phase.

## Analyst story
As an analyst / I want a validated R reimplementation of the cluster-anchoring scorer that reads the
existing standardized inputs and configs unchanged / so that I can build the generalized, Dockerized,
sensitivity-sweeping tool on a foundation proven equal to the Python reference.

## Pipeline position
- **Upstream:** `FinnGen_PheWAS_RG` stage → `cluster_trait_rg_long_with_p.tsv` (the rg long-table);
  `data/finngen_rg/finngen_R12_FIN.ldsc.summary.tsv` (trait×trait rg); `cluster_anchoring/ontology/*.tsv`.
- **Downstream:** AnchorMap Phase 2 (input generalization) builds on these R modules; ultimately the
  Nextflow `ANCHORMAP` process.
- **Reference engine to port (the spec):** `UKBB_CLUSTER_GWAS/scripts/cluster_anchoring/anchor_categories.py`
  (548 lines, read in full). **Orchestration to mirror later:** project R scripts use `argparse` +
  `data.table::setDTthreads` (`scripts/genomic_sem/scripts/run_cluster_gpca.R`); the Python tool's CLI is
  `--config <yaml>` only — **the R port reads the identical YAML so existing configs work unchanged.**

---

## CONTEXT REFERENCES — READ BEFORE IMPLEMENTING

### Input schemas (verified from real headers / code)

**Input A — rg long-table** `…/FinnGen_PheWAS_RG/results/carey_rint_tuned_15clusters_neff_max_empirical_covz/rg/cluster_trait_rg_long_with_p.tsv`
(42,795 rows; 34 cols). Read as character then coerce — **do not trust precomputed `z`/`abs_rg` cols; the
engine recomputes them.** Required columns (others tolerated/passed-through):

| col | type | role |
|---|---|---|
| `cluster_label` | str | group key (clusters scored independently) |
| `trait_id` | str | trait id; join key for lab/anthro ontology |
| `trait_category` | str | join key for disease ontology |
| `trait_group` | str | universe selector (`disease` / `lab_value`) |
| `rg` | float | genetic correlation (signed; clip ±0.999 for y/v, ±1 for matrix) |
| `rg_se` | float | SE; gate requires `>0` |
| `p` | float | rg p-value (ORA layer only) |
| `h2_trait` | float | trait h² |
| `h2_trait_se` | float | SE; gate requires `>0`; defines `h2_z` |
| `ldsc_converged` | `TRUE`/`FALSE` str | gate (if `require_ldsc_converged`) |
| `negative_h2` | `TRUE`/`FALSE` str | gate drop (if `drop_negative_h2`) |
| `status` | str | `success` required |

Coercion rules (from `load_long`, py L67–74): `pd.to_numeric(errors="coerce")` → R `as.numeric` (non-numeric→NA);
booleans = `toupper(x)=="TRUE"`; missing = empty/`NA`.

**Input B — trait×trait rg matrix** `…/UKBB_CLUSTER_GWAS/data/finngen_rg/finngen_R12_FIN.ldsc.summary.tsv`
(LDSC `--rg` summary). Header: `p1 p2 rg se z p h2_obs h2_obs_se h2_int h2_int_se gcov_int gcov_int_se CONVERGED`.
Engine uses only `p1, p2, rg, CONVERGED`. Build (from `load_trait_rg_matrix`, py L138–156): filter to gated
trait_ids; if `trait_rg_require_converged` keep `CONVERGED==TRUE`; `rg=clip(±1)`; pivot to symmetric matrix,
**average duplicate pairs**, `combine_first(t)` to symmetrize (file stores each pair once), reindex to sorted
trait set, set **diag=1**, missing pairs stay **NaN**.

**Input C — ontology TSV.** Disease `category_ontology.tsv` cols: `trait_category, domain, icd_chapter, kind,
anchor_eligible, notes` (join on `trait_category`). Anthro `anthro_ontology.tsv` cols: `trait_id, trait_label,
anthro_class, anchor_eligible` (join on `trait_id`; maps only `BMI_IRN/WEIGHT_IRN/HEIGHT_IRN` →
`Anthropometric`, eligible=TRUE — every other disease trait is left unmapped→NaN→skipped, with a benign
stderr WARN). `anchor_eligible` parsed `toupper=="TRUE"`; default TRUE if column absent. `native` level
aliases the join key.

**Input D — config YAML** (read identical files; keys + defaults from `load_config`, py L35–60):
`run_label, rg_long, ontology, out_dir, trait_group(=disease), require_ldsc_converged(=T), drop_negative_h2(=T),
h2_z_threshold(=4.0), levels(=[native,domain,icd_chapter]), primary_level(=domain), ontology_key(=trait_category),
min_category_n(=3), rank_variable(=abs_z), permutation_K(=2000), random_seed(=1), vif_min_rho(=0.0),
vif_correlation(=cluster_profile), trait_rg_matrix, trait_rg_require_converged(=T), hit_abs_rg(=0.2),
hit_bonferroni(=T), label_q_max(=0.05), label_auc_min(=0.60), shape_margin_sharp(=0.10),
shape_margin_diffuse(=0.05), shape_focus_diffuse(=3.0), validation{…}`.
Path resolution (py L483–491): inputs resolve relative to **stage_root** = config dir, or its parent if the
config lives in a `configs/` subdir.

### Output schema (contract — must match Python exactly)

`<out_dir>/category_anchor_scores.tsv` — one row per (cluster, level, category), sorted by `(level, cluster_label, rank)`:
```
cluster_label  level  category  eligible  n  n_eff  n_hit  rho_bar  vif  auc_abs  auc_signed
perm_p  vif_z  vif_p  pooled_rg  pooled_rg_ci_lo  pooled_rg_ci_hi  coherence
mean_abs_rg  mean_signed_rg  odds_ratio  fisher_p  q  rank
```
`<out_dir>/cluster_anchor_labels.tsv` — one row per cluster:
```
cluster_label  auto_label  anchor_shape  anchor_margin  anchor_focus  n_sig_domains
top_auc  top_q  top_pooled_rg  top_coherence  profile
```
Rounding (must replicate; R and Python 3 both use round-half-to-even, so `round()` matches): `n_eff`→2,
`rho_bar`→3, `vif`→2, `auc_abs/auc_signed`→4, `vif_z`→3, `pooled_rg`/`ci`→4, `coherence`→3, `mean_*_rg`→4,
`odds_ratio`→3, label `top_*` as shown. **Not rounded (full precision):** `perm_p, vif_p, fisher_p, q,
top_q`. `eligible` printed as Python bool `True/False`.

### Reference data & methods
- **Trait×trait rg:** FinnGen **R12** FIN LDSC summary (Finnish LD scores), `data/finngen_rg/…ldsc.summary.tsv`.
- **Ontology:** `cluster_anchoring/ontology/{category,anthro,lab}_ontology.tsv` (committed; copy into AnchorMap).
- **Method:** Li & Ji (2005) n_eff via **`poolr::meff(R, method="liji")`**; CAMERA VIF; Mann–Whitney/Wilcoxon
  AUC; IVW Fisher-z pooling; BH-FDR. Canonical impl = the Python functions below (the bit-for-bit spec):

| Python (`anchor_categories.py`) | R equivalent | gotcha |
|---|---|---|
| `apply_universe_gate` L77–102 | filter + `h2_z=h2_trait/h2_trait_se>z`; `y=atanh(clip(rg,±.999))`, `v=rg_se^2/(1-clip^2)^2` | recompute z/abs_rg/abs_z |
| `meff_li_ji` L159–173 | `poolr::meff(Rc,"liji")` on cleaned `Rc` (NaN→0, diag=1, symmetrize) | **clip eigenvalues ≥0 to match `np.clip`**; assert == manual `sum((λ≥1)+(λ-floor(λ)))` |
| `rho_bar` L176–184 | mean of finite upper-triangle off-diag | |
| `build_trait_profile_corr` L130–135 | `cor(t(pivot), use="pairwise.complete.obs")` then **mask pairs with <3 overlap → NA** | pandas `min_periods=3` |
| `auc_from_ranks` L190–192 | `U=sum(rank_in)-n_in(n_in+1)/2; AUC=U/(n_in·n_out)` | `rank()` default ties="average" == scipy `rankdata` |
| `perm_null_sums` L195–202 | `K` draws `sample(N,n_in)` of rank-sums; cache by `n_in`; seed | **RNG stream ≠ numpy → perm_p not bit-identical** |
| VIF block L242–248 | `var0=(N+1)/(12·n_in·n_out); z_un=(AUC-.5)/sqrt(var0); vif_z=z_un/sqrt(VIF); vif_p=pnorm(vif_z,lower=F)` | |
| pooled rg L250–260 | `ybar=sum(w·y)/sum(w),w=1/v; pooled_rg=tanh(ybar); ci=tanh(ybar±1.96·sqrt(VIF/sum(w)))` | |
| Fisher ORA L262–267 | `fisher.test(matrix(c(a,b,c,d),2,byrow=T),alternative="greater")$p.value`; **OR = `(a*d)/(b*c)` (sample OR), `Inf` if c==0, `0` if a==0** | **do NOT use `$estimate` (conditional MLE ≠ scipy sample OR)** |
| `bh_fdr` L285–295 | `p.adjust(perm_p, method="BH")` within (cluster, level) | verify monotone-min equivalence |
| `anchor_shape` L301–329 | n_sig/margin/focus rules; first-match weak→sharp→diffuse→focal | `focus=1/Σpᵢ²`, pᵢ∝max(AUC-.5,0)·sig |
| `rank_and_label` L332–374 | rank eligible by `(q↑,auc↓)`; label gate; top-8 profile string | non-eligible rank=NaN, never labelled |

- **Assumptions:** rg is on the LDSC observed scale; trait×trait matrix is a valid correlation source;
  ranks comparable. **Failure modes:** non-PD trait correlation matrix (negative eigenvalues → clip),
  zero-variance pooling (`v=0` when |rg|→1 — clip guards), category with `n_out<1` skipped.

### Files to read / create
- READ: `…/cluster_anchoring/anchor_categories.py` (full) — the line-referenced spec above.
- READ: `…/cluster_anchoring/configs/carey_rint15_anthro.yaml` & `carey_rint15.yaml` — configs to reuse.
- READ: `…/cluster_anchoring/output/carey_rint15_anthro/{category_anchor_scores,cluster_anchor_labels}.tsv`
  — the committed **oracle** (C5_sub0 values below).
- READ: `…/cluster_anchoring/docs/approach.md` — method rationale (not the spec; code is).
- CREATE: `AnchorMap/R/{io,gate,redundancy,score,label}.R`, `AnchorMap/anchor_map.R` (CLI), `AnchorMap/config/*.yaml`,
  `AnchorMap/tests/testthat/*`, `AnchorMap/tests/fixtures/*`, `AnchorMap/validation/{run_oracle.sh,compare_oracle.R}`,
  `AnchorMap/README.md`.

---

## METHOD / IMPLEMENTATION PLAN

### Phase 1: Inputs & harmonization (`R/io.R`)
Implement `read_long_tsv(path)`, `read_trait_rg_summary(path, traits, require_converged)`,
`read_ontology(path, key, levels)`, `read_config(path)` (yaml + `setdefault`s + stage_root resolution).
Each asserts its required columns and errors with the missing-column name. Mirror Python coercion exactly.

### Phase 2: Core analysis (`R/gate.R`, `R/redundancy.R`, `R/score.R`, `R/label.R`)
Port every function in the table above, 1:1, single-z (z = `h2_z_threshold` from config). Build correlation
matrix per the config `vif_correlation` flag (`trait_rg` → Input B builder; `cluster_profile` → proxy builder)
— **no automatic fallback in Phase 1.** Group by `cluster_label`, loop `levels`, `score_cluster_level`,
then `rank_and_label`.

### Phase 3: Outputs & integration (`anchor_map.R`)
CLI `Rscript anchor_map.R --config <yaml>` sources `R/*.R`, runs the pipeline, writes the two TSVs in the
exact contract schema + column order + rounding, and prints the same `[load]/[gate]/[vif]/[write]` progress
lines Python prints (foundation for the Phase-1 log). Write a minimal `anchormap.log` ending in `FINISHED`
(full structured log is later phases; here just prove the pattern).

### Phase 4: Validation
Unit tests per function (analytic) + real-data oracle parity on a cluster subset (integration).

---

## STEP-BY-STEP TASKS (execute top to bottom; each atomic + checkable)

### CREATE AnchorMap project skeleton
- **IMPLEMENT**: `R/` (modules), `anchor_map.R` (entry), `config/`, `tests/testthat/`, `tests/fixtures/`,
  `validation/`, `README.md`. Add `DESCRIPTION`-style header comment listing deps: R≥4.4, `poolr`,
  `data.table`, `yaml`, `argparse`, `testthat`.
- **PATTERN**: ADD §7.2 layout; mirror `run_cluster_gpca.R` for `argparse`+`setDTthreads` idiom.
- **VALIDATE**: `ls AnchorMap/R/*.R` shows io, gate, redundancy, score, label.

### CREATE `R/io.R` — readers + schema asserts
- **IMPLEMENT**: the four readers (Phase-1 §1). `read_long_tsv` uses `data.table::fread(colClasses="character")`
  then coerces the numeric/boolean cols; assert required cols present.
- **DATA/SCHEMA**: Input A/B/C/D contracts above (in: long-TSV/LDSC-summary/ontology/yaml; out: typed data.table + matrix).
- **GOTCHA**: ignore precomputed `z`,`abs_rg`; missing `n_cases`/`n_controls` empty strings → NA (not 0).
- **VALIDATE**: `Rscript -e 'source("R/io.R"); d<-read_long_tsv(LONG); stopifnot(is.numeric(d$rg), nrow(d)==42795)'`

### CREATE `R/gate.R` — reliability gate + per-trait stats
- **IMPLEMENT**: `apply_universe_gate(df,cfg)` → filter (status/converged/h2_z>z) + add `abs_rg,z,abs_z,y,v`.
- **PATTERN**: py L77–102. **GOTCHA**: clip rg to ±0.999 *inside* y/v only; `h2_z=h2_trait/h2_trait_se`.
- **VALIDATE**: on the anthro config, gated C5_sub0 BMI_IRN row has `rg=0.8085, y=atanh(0.8085)=1.1276,
  v=0.0300751^2/(1-0.8085^2)^2`; assert to 1e-6.

### CREATE `R/redundancy.R` — n_eff + ρ̄ + both correlation builders
- **IMPLEMENT**: `build_trait_rg_matrix` (Input B), `build_trait_profile_corr` (proxy, with <3-overlap masking),
  `meff_liji` (poolr::meff on cleaned matrix + clip eigenvalues ≥0; assert == manual), `rho_bar`.
- **GOTCHA**: NaN→0 + diag=1 + symmetrize before `poolr::meff`; clip negative eigenvalues to match numpy.
- **VALIDATE**: on `R=[[1,.9,.2],[.9,1,.4],[.2,.4,1]]` (cheat-sheet) → `meff_liji==2.00`, `rho_bar==0.50`,
  `VIF=1+(2-1)*0.5==1.5`; on the **real** BMI/WT/HEIGHT trait_rg submatrix → `n_eff==3.0`, `rho_bar==0.443`, `vif==1.89`.

### CREATE `R/score.R` — AUC, perm_p, VIF, pooled rg, ORA
- **IMPLEMENT**: `auc_from_ranks`, `perm_null_sums` (cache by n_in; `set.seed(random_seed)`), `score_cluster_level`.
- **PATTERN**: py L190–282. **GOTCHA**: **sample OR `(a*d)/(b*c)`** not `fisher.test$estimate`; `var0=(N+1)/(12·n_in·n_out)`;
  cache perm null by n_in; skip categories with `n_in<min_category_n` or `n_out<1`.
- **VALIDATE**: C5_sub0 anthro → `auc_abs==0.9164`, `vif_z==1.813`, **`vif_p==0.03489021688956177` (exact)**,
  `pooled_rg==0.2473`, `ci==[0.1965,0.2968]`, `coherence==1.0`, `odds_ratio==3.523`, `fisher_p==0.300211764776359`.

### CREATE `R/label.R` — BH-FDR, rank, auto-label, anchor shape
- **IMPLEMENT**: `bh_fdr`(=`p.adjust("BH")`), `rank_and_label`, `anchor_shape`.
- **PATTERN**: py L285–374. **GOTCHA**: rank only eligible cats by `(q↑,auc↓)`; label gate
  `q<0.05 & auc≥0.60 & vif_z>0 & vif_p<0.05 & n≥3`; shape first-match order weak→sharp→diffuse→focal.
- **VALIDATE**: C5_sub0 anthro → `auto_label=="Anthropometric"`, `anchor_shape=="sharp"`, `n_sig_domains==1`,
  `anchor_focus==1.0`, `anchor_margin==NA` (single category), `rank==1`.

### CREATE `anchor_map.R` — CLI + orchestration + outputs
- **IMPLEMENT**: `argparse` `--config`; source modules; run; write both TSVs (exact col order + rounding +
  `True/False` for eligible); print progress; write `anchormap.log` ending `FINISHED`.
- **GOTCHA**: write full-precision for perm_p/vif_p/fisher_p/q; round the rest; sort scores `(level,cluster,rank)`.
- **VALIDATE**: `Rscript anchor_map.R --config config/carey_rint15_anthro.yaml` produces 2 TSVs + log with `FINISHED`.

### CREATE `tests/fixtures/` — analytic + oracle-subset fixtures
- **IMPLEMENT**: (a) cheat-sheet 3×3 R-matrix + the analytic targets; (b) **cluster-subset fixture**: filter the
  real long-TSV to `{C5_sub0, C4, C5_sub1, C0, noise_re0}` (clusters scored independently ⇒ exact same per-cluster
  scores) + the trait_rg summary subset to those clusters' gated trait_ids; copy the matching Python oracle rows.
- **GOTCHA**: do NOT subset *traits within a cluster* (changes N → changes AUC); only subset whole clusters.
- **VALIDATE**: fixture long-TSV has only the 5 cluster_labels; oracle rows present for each.

### CREATE `validation/compare_oracle.R` + `run_oracle.sh`
- **IMPLEMENT**: per-column comparator: deterministic cols exact to tol (`1e-6`, or rounded-decimals);
  `perm_p`/`q` within MC tol `2*sqrt(p(1-p)/K)` AND same significance call; `auto_label`/`anchor_shape`/`rank`
  identical. `run_oracle.sh`: (1) run Python ref on the full anthro+disease configs, (2) run R, (3) compare.
- **VALIDATE**: `bash validation/run_oracle.sh` → "0 deterministic mismatches; labels identical".

### CREATE `tests/testthat/test-engine.R`
- **IMPLEMENT**: unit asserts (gate stats, meff/rho/VIF, AUC, sample-OR, BH, shape) + the C5_sub0 oracle row + a
  negative check (a `weak` cluster, e.g. C3 anthro `auto_label=="ambiguous"`, `auc_abs==0.0653`).
- **VALIDATE**: `Rscript -e 'testthat::test_dir("tests/testthat")'` → all pass.

---

## VALIDATION STRATEGY
- **Unit (analytic):** each ported function vs hand-computed values incl. the cheat-sheet R-matrix (`n_eff=2.00,
  VIF=1.5, vif_z≈1.768, vif_p≈0.0385`).
- **Positive control:** **C5_sub0 anthro → "Anthropometric [sharp]"** with the exact real-run values above
  (deterministic cols exact; perm_p≈0.0055 within MC error).
- **Negative control:** anthro track C3 → `auc_abs=0.0653`, `q≈0.999`, `auto_label=ambiguous`, shape `weak`;
  **forbidden-FP:** in the *disease* track, `Quantitative` (anchor_eligible=FALSE) must never be C5_sub0's label.
- **Integration parity:** R vs committed Python oracle over the full anthro track AND a disease-track cluster
  subset — every deterministic column matches; labels/shapes/ranks identical.
- **Sensitivity (smoke):** run with `permutation_K: 200` → labels unchanged, perm_p within wider MC band.
- **Schema/sanity:** output col order/names exact; `n_eff≤n`; `vif≥1`; no NaN explosion; `eligible∈{True,False}`.

## VALIDATION COMMANDS (run all; zero deterministic/control failures)
```bash
# 0. deps (local R or the Phase-5 container once built)
Rscript -e 'stopifnot(all(c("poolr","data.table","yaml","argparse","testthat") %in% rownames(installed.packages())))'

# 1. (re)generate the Python oracle in the reference repo
cd "…/UKBB_CLUSTER_GWAS/scripts/cluster_anchoring"
python3 anchor_categories.py --config configs/carey_rint15_anthro.yaml
python3 anchor_categories.py --config configs/carey_rint15.yaml

# 2. run the R port on the same configs
cd "…/AnchorMap"
Rscript anchor_map.R --config config/carey_rint15_anthro.yaml
Rscript anchor_map.R --config config/carey_rint15.yaml

# 3. cross-language parity (deterministic exact; perm_p/q tolerant; labels identical)
Rscript validation/compare_oracle.R \
  --r-out   results/carey_rint15_anthro/category_anchor_scores.tsv \
  --oracle  "…/cluster_anchoring/output/carey_rint15_anthro/category_anchor_scores.tsv"

# 4. unit + fixture tests
Rscript -e 'testthat::test_dir("tests/testthat")'
```

## ACCEPTANCE CRITERIA
- [ ] All inputs assert their schema; both output TSVs match the Python column contract (order, names, rounding, `True/False`).
- [ ] **Positive control recovered:** C5_sub0 → Anthropometric [sharp], deterministic values exact (`vif_p=0.03489…`, `pooled_rg=0.2473`, `auc=0.9164`).
- [ ] **Negative control null:** C3 anthro `ambiguous/weak`; disease-track `Quantitative` never labels C5_sub0.
- [ ] R vs Python: 0 deterministic-column mismatches across the full anthro track + disease cluster subset; perm_p within MC tol; labels/shapes/ranks identical.
- [ ] `poolr::meff` n_eff == manual eigen implementation == Python on the fixtures.
- [ ] Smoke run (`K=200`) passes; full run commands documented; `anchormap.log` ends with `FINISHED`.
- [ ] Existing YAML configs run unchanged (only paths adjusted to AnchorMap layout); provenance (config, seed, input paths) logged.

## NOTES
- **perm_p is the one quantity that cannot be bit-identical across languages** (numpy PCG64 ≠ R RNG). Anchor the
  gate on the deterministic `vif_p` (full precision, exact) + label stability; treat perm_p distributionally.
  If exact perm_p parity is ever required, add an optional record/replay of the Python permutation index draws
  (deferred — not needed for Phase 1).
- **scipy `fisher_exact` returns the sample OR `(a·d)/(b·c)`**, not R's conditional-MLE `fisher.test$estimate` —
  compute OR manually; only borrow `$p.value`. This is the single most likely silent mismatch.
- Clusters are scored independently → safe to subset whole clusters for fast fixtures; never subset traits within a cluster.
- Phase 1 honors `vif_correlation` from config verbatim (anthro/disease use `trait_rg`; lab uses `cluster_profile`).
  The automatic coverage-based fallback and GenomicSEM `.rds` ingestion are **Phase 2**.
- Keep `poolr::meff` as the documented n_eff path (ADD decision) but guard it with the Python-matching matrix
  cleaning + eigenvalue clip so the parity gate holds on non-PD matrices.
