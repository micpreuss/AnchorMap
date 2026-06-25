# Analysis: Phase 6 — AUC confidence intervals + shape confidence score

Validate the schema contracts and method usage against the actual code before implementing.
Pay special attention to column **order/append position**, the `auc_abs` vs `auc_signed` distinction,
the VIF inflation convention, and the RNG-reproducibility discipline (Mersenne-Twister pin).

## Question & object

**Scientific question.** Two point estimates currently drive every AnchorMap call but carry no
uncertainty: (1) the competitive **AUC** per (cluster, level, category) — the Mann–Whitney
probability-of-superiority that an in-category trait outranks an out-category trait
([R/score.R:43-44](../../R/score.R#L43-L44)); and (2) the **anchor_focus** inverse-Simpson index that,
with `n_sig`/`margin`, decides the discrete **anchor_shape** ∈ {weak, sharp, focal, diffuse}
([R/label.R:20-38](../../R/label.R#L20-L38)). Phase 6 adds:

- **A — a confidence interval for each AUC** (`auc_abs_ci_lo/hi`, plus the SE `auc_abs_se`), so a high
  AUC built on 3 correlated traits is visibly less certain than one on 30 independent traits.
- **B — a confidence score for the shape call** (`shape_confidence`, `anchor_focus_ci_lo/hi`,
  `shape_posterior`, `shape_jackknife_stable`), answering: *if the AUCs wobble within their intervals,
  how often does the sharp/focal/diffuse/weak verdict flip?*

**Object.** The two scored TSVs the engine already writes — `category_anchor_scores.tsv` (one row per
cluster×level×category) and `cluster_anchor_labels.tsv` (one row per cluster). Phase 6 **appends**
columns; it changes **no** existing column value.

**Inference licensed.** These quantify **sampling uncertainty in `rg` given `rg_se` / VIF** —
"how stable is this anchor under rg noise" — *not* ontology error, model misspecification, or
selection. They are estimation intervals; they do **not** replace the `vif_p`/`q` testing machinery.

**Classification.** Extension of a stage (the scoring + labelling engine). **Complexity: Medium.**
Main risk is **method/edge-case** (degenerate variances at perfect separation / tiny n) and
**reproducibility** (the Part-B Monte-Carlo must stay thread- and order-invariant and must not perturb
the byte-for-byte `perm_p` parity that Phases 1–3 guarantee).

## Analyst story

As an analyst,
I want each AUC to carry a redundancy-aware confidence interval and each anchor-shape call to carry a
support score,
So that I can tell a sharp, well-supported anchor from one that is an artefact of a handful of
correlated traits, and report intervals rather than bare point estimates.

## Pipeline position

- **Upstream:** unchanged — gate → redundancy → `score_cluster_level()` → `rank_and_label()`
  ([R/score.R](../../R/score.R), [R/label.R](../../R/label.R)), driven by `run_anchormap()` and the
  parallel z-sweep `run_sensitivity()`/`score_at_z()` ([R/run_anchormap.R](../../R/run_anchormap.R),
  [R/sensitivity.R](../../R/sensitivity.R)).
- **Downstream:** the scored TSVs → `plot.R` figures; the oracle comparator
  ([validation/compare_oracle.R](../../validation/compare_oracle.R)); the sensitivity TSVs.
- **Orchestration pattern to mirror:** the existing `pooled_rg_ci_lo/hi` CI (already emitted from
  `score_cluster_level` at [R/score.R:60-65](../../R/score.R#L60-L65)) is the template for Part A;
  the MC-reproducibility pattern at [R/sensitivity.R:55-60](../../R/sensitivity.R#L55-L60)
  (`set.seed(seed, kind="Mersenne-Twister", normal.kind="Inversion", sample.kind="Rejection")`) is the
  template for Part B.

## CONTEXT REFERENCES — READ BEFORE IMPLEMENTING

### Input schemas (already in-memory; no new files read)

Part A works **inside** `score_cluster_level(gc, level, corr, cfg)`, which already holds everything
needed — no new inputs:
- `ranks_abs <- rank(rv, ties.method="average")` — pooled midranks on the rank variable `rv`
  (= `cfg$rank_variable`, e.g. `abs_z`) ([R/score.R:25-26](../../R/score.R#L25-L26)). **The DeLong
  placement values derive from these midranks — no rescan of the data.**
- `inmask`, `n_in`, `n_out`, `auc_abs` ([R/score.R:39-43](../../R/score.R#L39-L43)).
- `vif` ([R/score.R:54](../../R/score.R#L54)) — the redundancy inflation factor, `vif ≥ 1` guaranteed
  (`vif_min_rho ≥ 0`, `m_eff ≥ 1`, `rho_bar ≥ 0`).

Part B works **inside** `rank_and_label()` / `anchor_shape()`, which hold the per-cluster
primary-level eligible rows with `auc_abs`, the new `auc_abs_se`, and `q`
([R/label.R:20-38](../../R/label.R#L20-L38), [R/label.R:57-74](../../R/label.R#L57-L74)).

### Output schema (contract for downstream)

**`category_anchor_scores.tsv`** — `.SCORE_COLS` at
[R/run_anchormap.R:4-7](../../R/run_anchormap.R#L4-L7). **APPEND three columns at the END** (preserving
all existing positions):

| col | meaning | precision |
|---|---|---|
| `auc_abs_se` | VIF-inflated DeLong SE of `auc_abs` | full precision |
| `auc_abs_ci_lo` | lower 95% CI bound (logit-back-transformed) | round 4 (matches `auc_abs`) |
| `auc_abs_ci_hi` | upper 95% CI bound | round 4 |

Invariant: `0 ≤ ci_lo ≤ auc_abs ≤ ci_hi ≤ 1`.

**`cluster_anchor_labels.tsv`** — `.LABEL_COLS` at
[R/run_anchormap.R:8-9](../../R/run_anchormap.R#L8-L9). **APPEND at the END:**

| col | meaning | precision |
|---|---|---|
| `shape_confidence` | fraction of MC draws whose shape == the point shape | round 3 |
| `anchor_focus_ci_lo` | 2.5% quantile of `focus*` over significant draws (NA if none) | round 2 (matches `anchor_focus`) |
| `anchor_focus_ci_hi` | 97.5% quantile of `focus*` | round 2 |
| `shape_posterior` | compact string, e.g. `"sharp=0.81;focal=0.15;diffuse=0.04;weak=0.00"` | string |
| `shape_jackknife_stable` | `True/False`: shape unchanged when any single significant domain is dropped | bool repr |

**Sensitivity contracts** `.SENS_SCORE_COLS` / `.SENS_LABEL_COLS`
([R/run_anchormap.R:11-12](../../R/run_anchormap.R#L11-L12)) extend automatically (`= .SCORE_COLS + …`),
so the new columns flow through the z-sweep TSVs with no extra wiring beyond the constant update.

### Reference data & methods

- **AUC variance — DeLong, DeLong & Clarke-Pearson (1988), *Biometrics* 44:837.** Nonparametric AUC
  variance from per-observation *placement values* (structural components). Fast form (Sun & Xu 2014,
  *IEEE SPL* 21:1389) computes them from midranks in O(N log N): for in-group element *i* with pooled
  midrank `R_i` and within-in-group midrank `Q_i`, the placement value is
  `V10_i = (R_i − Q_i)/n_out`; for out-group *j*, `V01_j = 1 − (R_j − Q_j)/n_in`. Then
  `S10 = var(V10)`, `S01 = var(V01)` (sample variances), and
  **`Var(AUC) = S10/n_in + S01/n_out`**. Canonical impl: `pROC::ci.auc(method="delong")`
  (Robin et al. 2011, *BMC Bioinformatics* 12:77) — used only as a **cross-check oracle in tests**, not
  a runtime dependency.
  - *Assumption:* independent observations. **Violated here** (correlated traits) → inflate by VIF
    (below).
  - *Failure mode:* perfect separation (AUC=0 or 1) → `Var=0` → zero-width CI. Handled by clamp +
    Hanley–McNeil fallback (see edge cases).
- **VIF inflation.** DeLong assumes independence; AnchorMap's whole `vif`/`n_eff` machinery exists
  because in-category traits are correlated. Inflate consistently with the existing test:
  `vif_z` divides z by `√vif` ([R/score.R:57](../../R/score.R#L57)), so the CI uses
  **`Var_adj = vif · Var(AUC)`**. Deterministic; never shrinks the interval (`vif ≥ 1`).
- **Logit transform (Newcombe 2006, *Stat Med* 25:559).** Keep the CI in (0,1) — essential because the
  positive control sits at `auc_abs = 0.9164`. `se_logit = √Var_adj / (a(1−a))` (delta method);
  bounds `plogis(qlogis(a) ± 1.96·se_logit)`.
- **Null vs alternative variance (the crux).** The existing `var0 = (N+1)/(12·n_in·n_out)`
  ([R/score.R:55](../../R/score.R#L55)) is the variance of U **under H₀: AUC=0.5** — correct for the
  `vif_p` *test*, **wrong for a CI around the observed AUC**. Part A must use the DeLong
  (alternative-hypothesis) variance. Do **not** reuse `var0` for the CI.
- **Hanley–McNeil (1982), *Radiology* 143:29 — variance fallback.** Closed form, strictly positive for
  clamped `a∈(0,1)`: `Q1 = a/(2−a)`, `Q2 = 2a²/(1+a)`,
  `Var_HM = [a(1−a) + (n_in−1)(Q1−a²) + (n_out−1)(Q2−a²)] / (n_in·n_out)`. Used **only** when DeLong
  `Var_adj` is 0 / non-finite (perfect separation, degenerate groups), then VIF-inflated and logit-
  transformed identically.
- **Inverse-Simpson / Hill number ²D (Hill 1973; Jost 2006, *Oikos* 113:363).** `anchor_focus` is the
  effective number of anchored domains; its uncertainty is propagated by MC, **not** the classic
  multinomial diversity-index variance (Nayak 1985), because the weights are `max(auc−0.5,0)` not
  counts and the shape verdict is a discrete threshold function of several quantities.
- **Bootstrap support analogy (Felsenstein 1985, *Evolution* 39:783).** `shape_confidence` = fraction
  of MC replicates recovering the point shape = a posterior support / selection-stability for a
  discrete label.

### Files to read / create

- READ: [R/score.R](../../R/score.R) (1-92) — Why: add Part A inside `score_cluster_level`; mirror the
  `pooled_rg_ci` block (L60-65) and the `vif` inflation idiom.
- READ: [R/label.R](../../R/label.R) (1-91) — Why: add Part B; `anchor_shape()` is the ruleset the MC
  must re-evaluate per draw, `rank_and_label()` is where per-cluster MC + jackknife live.
- READ: [R/run_anchormap.R](../../R/run_anchormap.R) (4-12, 118-142) — Why: update the four column
  contracts and the write/rounding/`eligible`-bool-repr path.
- READ: [R/sensitivity.R](../../R/sensitivity.R) (55-77) — Why: the RNG pin + serial-within-z pattern
  Part B must follow; confirm the new label MC does not disturb the per-z `perm_p` stream (it runs
  *after* all `perm_p` draws and re-seeds independently).
- READ: [validation/compare_oracle.R](../../validation/compare_oracle.R) (16-31) — Why: confirm the
  comparator is **column-name-driven** (iterates `names(tol)`), so appended R-only columns are ignored
  → **no comparator change needed**; the new columns simply have no Python counterpart.
- READ: [tests/testthat/helper-fixtures.R](../../tests/testthat/helper-fixtures.R),
  [tests/testthat/test-sensitivity.R](../../tests/testthat/test-sensitivity.R),
  [tests/testthat/test-validation.R](../../tests/testthat/test-validation.R) — Why: mirror the
  fixture/`eqnum`/parity test idioms; locate the column-contract and positive-control assertions to
  extend.
- CREATE: `R/uncertainty.R` — the DeLong/Hanley–McNeil variance + logit CI (`auc_delong_var`,
  `auc_ci_logit`, `auc_hanley_var`) and the shape-MC + jackknife helpers (`shape_confidence_mc`,
  `shape_jackknife`). Keeps `score.R`/`label.R` slim; they call into it.
- CREATE: `tests/testthat/test-uncertainty.R` — Part A/B unit + edge-case + reproducibility tests.
- UPDATE: `inst/configs/example_disease.yaml`, `example_anthro.yaml`, `synthetic_rds.yaml` — add the
  Phase-6 config block (defaults below).
- UPDATE: `man/` (roxygen) for any newly `@export`ed function.

## METHOD / IMPLEMENTATION PLAN

### Phase 1: Config + module scaffold

Add a Phase-6 config block (all with safe defaults so existing configs keep working):

```yaml
# ---- Phase 6: uncertainty ----
emit_uncertainty: true        # master toggle; false -> legacy column contracts (parity mode)
auc_ci_method: delong         # delong (default) | hanley
auc_ci_level: 0.95
shape_confidence_B: 2000      # MC draws for shape support (mirrors permutation_K scale)
shape_confidence_min_sig_frac: 0.05   # min fraction of draws with >=1 sig domain to report a focus CI
```

`run_anchormap()` must read these with defaults (treat missing as the values above) so the shipped and
machine-specific configs that predate Phase 6 still run. When `emit_uncertainty: false`, the engine
emits the legacy `.SCORE_COLS`/`.LABEL_COLS` (11/24 cols) byte-for-byte — the Phase-1 parity escape
hatch and the Docker self-test's stable contract.

### Phase 2: Part A — deterministic AUC CI (in `score_cluster_level`)

Per category, immediately after `auc_abs` is computed ([R/score.R:43](../../R/score.R#L43)):

```r
# DeLong placement values from midranks already in hand (Sun & Xu 2014 fast form)
R_in  <- ranks_abs[inmask];  R_out <- ranks_abs[!inmask]
tx    <- rank(rv[inmask],  ties.method = "average")   # within in-group midrank
ty    <- rank(rv[!inmask], ties.method = "average")   # within out-group midrank
V10   <- (R_in  - tx) / n_out
V01   <- 1 - (R_out - ty) / n_in
S10   <- if (n_in  > 1) stats::var(V10) else 0
S01   <- if (n_out > 1) stats::var(V01) else 0
var_auc <- S10 / n_in + S01 / n_out
var_adj <- vif * var_auc                              # redundancy inflation (consistent w/ vif_z)
if (!is.finite(var_adj) || var_adj <= 0)              # perfect-separation / degenerate fallback
  var_adj <- vif * auc_hanley_var(auc_abs, n_in, n_out)
# logit CI with boundary continuity clamp
eps <- 1 / (2 * n_in * n_out)
a   <- min(max(auc_abs, eps), 1 - eps)
se_logit <- sqrt(var_adj) / (a * (1 - a))
zc  <- stats::qnorm(1 - (1 - auc_ci_level) / 2)
ci_lo <- stats::plogis(stats::qlogis(a) - zc * se_logit)
ci_hi <- stats::plogis(stats::qlogis(a) + zc * se_logit)
auc_abs_se <- sqrt(var_adj)
```

Add `auc_abs_se` (full precision), `auc_abs_ci_lo = round(ci_lo,4)`, `auc_abs_ci_hi = round(ci_hi,4)`
to the returned `data.frame` ([R/score.R:77-89](../../R/score.R#L77-L89)). **Do not touch** `auc_abs`,
`auc_signed`, `vif_z`, `vif_p`, `pooled_rg*`, or any other field.

### Phase 3: Part B — shape confidence (in `rank_and_label` / `anchor_shape`)

Factor `anchor_shape`'s verdict into a pure function `decide_shape(auc, sig, focus, margin, cfg)` so
both the point call and each MC draw use **one** ruleset (no divergence). Then, per cluster at the
primary level (inside the labelling loop, [R/label.R:69-74](../../R/label.R#L69-L74)):

```r
shape_confidence_mc <- function(sub, cfg) {           # sub = primary-level eligible rows, rank-sorted
  a  <- sub$auc_abs; se <- sub$auc_abs_se; q <- sub$q
  point <- decide_shape(...)                          # the deterministic verdict on point estimates
  if (!isTRUE(cfg$emit_uncertainty)) return(point only)
  set.seed(cfg$random_seed + match(cl, sorted_clusters),   # order- & thread-invariant per-cluster seed
           kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
  B <- cfg$shape_confidence_B
  draws <- replicate(B, {
    lo <- qlogis(pmin(pmax(a, eps), 1 - eps))
    astar <- plogis(lo + rnorm(length(a), 0, se / (a * (1 - a))))   # logit-normal; q held FIXED
    sig <- (q < cfg$label_q_max) & (astar >= cfg$label_auc_min)
    margin <- if (length(astar) > 1) sort(astar, TRUE)[1] - sort(astar, TRUE)[2] else NA
    w <- pmax(astar - 0.5, 0) * sig
    focus <- if (sum(w) > 0) { p <- w / sum(w); 1 / sum(p^2) } else NA
    list(shape = decide_shape(astar, sig, focus, margin, cfg), focus = focus)
  })
  list(shape_confidence = mean(shapes == point$anchor_shape),
       shape_posterior  = compact_table(shapes),                    # "sharp=..;focal=..;.."
       anchor_focus_ci  = quantile(focus[is.finite], c(.025,.975)) or NA if < min_sig_frac)
}
```

- **`q` is held fixed** per draw (only AUC is perturbed): perturbing significance would require
  re-permuting (a separate, expensive noise source). Document this as a scope choice — the score is
  "shape support conditional on the significance calls".
- AUC draws are **independent** across categories (a deliberate approximation; the true joint shares
  the complement set). Document as a limitation.
- **Jackknife** (deterministic, no RNG): `shape_jackknife()` drops each significant domain in turn,
  recomputes `decide_shape`; `shape_jackknife_stable = all(verdicts == point)`. Flags single-domain
  dependence (the sharp/diffuse distinction).

### Phase 4: Wiring + outputs

- Append the new columns to `.SCORE_COLS`, `.LABEL_COLS`, `.SENS_SCORE_COLS`, `.SENS_LABEL_COLS`
  ([R/run_anchormap.R:4-12](../../R/run_anchormap.R#L4-L12)) — **append only**, never reorder.
- Gate the appended columns behind `emit_uncertainty`: when false, project to the legacy contracts.
- The write path ([R/run_anchormap.R:118-142](../../R/run_anchormap.R#L118-L142)) needs the
  `True/False` bool-repr applied to `shape_jackknife_stable` (mirror the `eligible`/`label_stable`
  `ifelse(...,"True","False")` idiom); `shape_posterior` and `anchor_focus_ci_*` write as-is with
  `na = ""`.
- The z-sweep needs **no** logic change: Part A lives in `score_cluster_level` and Part B in
  `rank_and_label`, both already re-run per z by `score_at_z`. Confirm the per-cluster MC reseed keeps
  thread- and order-invariance (it re-seeds independently of the `perm_p` stream).

### Phase 5: Figure integration — surface the uncertainty on the lollipop small-multiples

Extend `fig_lollipops()` ([R/plot.R:108-152](../../R/plot.R#L108-L152)) so the new columns are
**visible**, not just tabular. Two additions, both **strictly opt-in on column presence** so the figure
still renders byte-for-byte on pre-Phase-6 TSVs (the `load_track` required-column list must **not** gain
the new columns):

1. **AUC 95% CI band on each stem.** A thin black horizontal bar with short vertical end-caps spanning
   `[auc_abs_ci_lo, auc_abs_ci_hi]` at the stem's height `y` — a horizontal error bar overlaying the
   coloured stem so the uncertainty around the AUC tip reads at a glance. The stem stays the wider
   (`linewidth = 1.1`) coloured layer; the CI bar is thin (`linewidth ≈ 0.4`) and black so the signed-rg
   colour still shows around it.

   ```r
   has_ci <- all(c("auc_abs_ci_lo","auc_abs_ci_hi") %in% names(s))   # gate on presence
   # ... build the base plot `g` (vline + coloured geom_segment), then BEFORE the tip points:
   if (has_ci)
     g <- g + geom_errorbar(
       aes(y = y, xmin = pmax(auc_abs_ci_lo, 0.5), xmax = pmin(auc_abs_ci_hi, 1.0)),
       orientation = "y", width = 0.32, colour = "black", linewidth = 0.4)
   # ... then the coloured tip point (size 3) on top, sig ring, auto-label star, labels (unchanged).
   ```

   - **Clamp** `xmin`/`xmax` to the visible `[0.5, 1.0]` x-axis (a low-AUC domain's `ci_lo` can fall
     below 0.5) so the bar is not silently dropped by the scale's `oob` censor — consistent with the
     existing deliberate axis crop.
   - Use `geom_errorbar(..., orientation = "y")`, **not** the soft-deprecated `geom_errorbarh`
     (ggplot2 4.0.x — the pinned image's version).
   - Layer order: `geom_segment` (coloured stem) → `geom_errorbar` (black CI) → coloured tip point →
     sig ring → auto-label star, so the estimate point stays crisp on top of the CI bar.

2. **`shape_confidence` in the panel title.** Append the numeric support after the shape category in the
   per-panel title's bracket — `"C5_sub0  -  Anthropometric [sharp, conf=0.87]"` — only when the
   `shape_confidence` column is present and finite:

   ```r
   has_conf <- "shape_confidence" %in% names(labels)
   shp <- lb[["anchor_shape"]]
   if (has_conf && is.finite(suppressWarnings(as.numeric(lb[["shape_confidence"]]))))
     shp <- sprintf("%s, conf=%.2f", shp, as.numeric(lb[["shape_confidence"]]))
   ttl <- sprintf("%s  -  %s [%s]", cl, lb[["auto_label"]], shp)
   ```

3. **Legend note (conditional).** When `has_ci`, append `"; black bar = 95% AUC CI"` to the
   `plot_annotation` title so the new glyph is documented; leave the annotation unchanged otherwise.

Scope: **lollipops only** for this phase (it has the per-domain stem the CI band naturally rides). The
dot-heatmap / scatter already encode AUC as a size/position channel where a CI band has no clean home;
leave them untouched (note as a deferred follow-up if a per-cell uncertainty glyph is later wanted).

## STEP-BY-STEP TASKS (execute top to bottom; each atomic + checkable)

### CREATE `R/uncertainty.R`
- **IMPLEMENT**: `auc_delong_var(ranks_abs, inmask, rv, n_in, n_out)`, `auc_hanley_var(a, n_in, n_out)`,
  `auc_ci_logit(auc, var_adj, n_in, n_out, level)` → `c(se, lo, hi)`; `decide_shape(auc, sig, focus,
  margin, cfg)` (extracted ruleset); `shape_confidence_mc(sub, cl, sorted_clusters, cfg)`;
  `shape_jackknife(sub, cfg)`.
- **PATTERN**: CI idiom from [R/score.R:60-65](../../R/score.R#L60-L65); RNG pin from
  [R/sensitivity.R:55-60](../../R/sensitivity.R#L55-L60).
- **GOTCHA**: sample variance needs group size > 1; `var0` (null) ≠ DeLong (alternative) variance.
- **VALIDATE**: `Rscript -e 'devtools::load_all("."); print(anchormap:::auc_hanley_var(0.9164,3,30))'`

### UPDATE `R/score.R` — emit Part A columns
- **IMPLEMENT**: insert the DeLong-CI block after L43; add `auc_abs_se/ci_lo/ci_hi` to the row frame.
- **DATA/SCHEMA**: produces 3 new score columns; asserts `vif ≥ 1` already holds.
- **GOTCHA**: clamp `a` with `eps = 1/(2 n_in n_out)`; HM fallback when `var_adj ≤ 0`.
- **VALIDATE**: rerun the synthetic engine (command below); check `ci_lo ≤ auc_abs ≤ ci_hi ∈ [0,1]`.

### UPDATE `R/label.R` — extract `decide_shape`, add MC + jackknife
- **IMPLEMENT**: refactor `anchor_shape` to call `decide_shape`; call `shape_confidence_mc` +
  `shape_jackknife` in `rank_and_label`'s cluster loop; add the 5 label columns.
- **GOTCHA**: weak clusters (`sum(w)==0`) → `focus=NA`, focus CI=NA, but `shape_confidence` still valid;
  per-cluster reseed must be order-independent (`random_seed + match(cl, sorted)`).
- **VALIDATE**: positive control C5_sub0 → `shape_confidence ≈ 1.0`, shape `sharp`.

### UPDATE `R/run_anchormap.R` — contracts + write path + `emit_uncertainty`
- **IMPLEMENT**: append to the 4 `.*_COLS`; bool-repr `shape_jackknife_stable`; legacy projection when
  `emit_uncertainty: false`.
- **GOTCHA**: append only; verify primary == sweep[z==primary] still holds for the new columns.

### UPDATE configs — add the Phase-6 block
- **IMPLEMENT**: add the YAML block to `example_disease.yaml`, `example_anthro.yaml`,
  `synthetic_rds.yaml`; read with defaults in `load_config`/`run_anchormap`.

### UPDATE `R/plot.R` — lollipop AUC CI band + shape_confidence in title
- **IMPLEMENT**: in `fig_lollipops`, add the presence-gated `geom_errorbar(orientation="y")` 95% AUC CI
  band (clamped to `[0.5,1]`), append `conf=%.2f` to the panel-title shape bracket, and the conditional
  legend note (Phase-5 of the implementation plan).
- **PATTERN**: [R/plot.R:108-152](../../R/plot.R#L108-L152) (`fig_lollipops`); mirror the existing
  layer/`plot_annotation` idioms.
- **DATA/SCHEMA**: consumes the new `auc_abs_ci_lo/hi` (scores) + `shape_confidence` (labels)
  **optionally** — do **not** add them to `load_track`'s `required_scores`/required-labels lists.
- **GOTCHA**: gate every new layer/title-fragment on `%in% names(...)` so pre-Phase-6 TSVs render
  byte-identically; use `geom_errorbar(orientation="y")`, not deprecated `geom_errorbarh`; clamp the CI
  to the visible axis so it is not censored away.
- **VALIDATE**: render lollipops twice — once on a current TSV (no band, no `conf=`), once on a
  Phase-6 TSV (band + `conf=` present); both succeed headless.

### CREATE `tests/testthat/test-uncertainty.R`
- **IMPLEMENT**: the unit/edge/reproducibility/parity tests in the strategy below.
- **PATTERN**: [tests/testthat/test-sensitivity.R](../../tests/testthat/test-sensitivity.R) idioms
  (`fx`, `eqnum`, `same`, thread-invariance).

### RUN validation
- **VALIDATE**: the commands block below; zero failures.

## VALIDATION STRATEGY

### Part A — AUC CI
- **Cross-check oracle:** for clean (non-degenerate) categories, `auc_ci` (VIF=1) must match
  `pROC::ci.auc(method="delong")` within 1e-6 (skip the test if `pROC` absent). This pins the DeLong
  math independently.
- **Containment + shape:** `ci_lo ≤ auc_abs ≤ ci_hi`, both ∈ [0,1], for **every** scored row.
- **VIF monotonicity:** with `vif > 1` the interval is **wider** than with `vif = 1` (same point AUC).
- **Positive control (n=3):** C5_sub0 anthro `auc_abs = 0.9164`, `n_in = 3` → a **wide** interval
  (assert width > 0.1 and `ci_hi` near, not at, 1.0) — small-n honesty, not a bug.
- **Invariance:** all pre-existing columns (`auc_abs`, `vif_z`, `vif_p`, `pooled_rg*`, `perm_p`, `q`,
  ranks, labels, shapes, `anchor_focus`) byte-identical to the Phase-5 output (the VIF-style invariant).

### Part B — shape confidence
- **Positive control:** C5_sub0 → `shape_confidence` ≈ 1.0, `sharp`; `shape_jackknife_stable = True`.
- **Negative control:** C3 anthro (weak, `auc_abs = 0.0653`) → shape `weak` with high
  `shape_confidence`; `anchor_focus_ci_*` = NA (no significant domains).
- **Reproducibility:** `shape_confidence` identical across `threads ∈ {1,4}` and independent of cluster
  iteration order (per-cluster MT reseed) — the Phase-3 thread-invariance test, extended.
- **`perm_p` non-disturbance:** the primary-slice `perm_p` column is byte-identical before/after Phase 6
  (the label-stage MC re-seeds independently and runs after all `perm_p` draws).
- **Soundness:** `shape_posterior` probabilities sum to 1; `shape_confidence` = the point shape's
  posterior entry; `anchor_focus_ci_lo ≤ anchor_focus ≤ anchor_focus_ci_hi` when non-NA.

### EDGE CASES TO HANDLE BEFOREHAND (the executor must implement guards for each)

**Part A (AUC CI):**
1. **Perfect separation `auc=0|1`** → DeLong `Var=0`, `logit(0|1)=±Inf`, zero-width CI. **Guard:** clamp
   `a∈[eps,1−eps]`, `eps=1/(2 n_in n_out)`; **Hanley–McNeil variance fallback** when `var_adj ≤ 0`/
   non-finite (strictly positive at clamped `a`).
2. **`n_out < 2` or `n_in < 2`** → placement-value sample variance undefined. **Guard:** the
   corresponding `S10`/`S01 = 0`; HM fallback supplies a finite floor. (`min_category_n` makes `n_in≥3`
   typical, but `n_out` can be 1–2.)
3. **Constant rank variable (all ties)** → `auc=0.5`, all `V=` equal, `Var=0`. **Guard:** HM fallback at
   `a=0.5` → finite small interval.
4. **VIF very large (tiny n, high ρ)** → raw bounds exceed [0,1]. **Guard:** logit transform keeps
   bounds in (0,1); the wide width is correct, not clamped away.
5. **NA/Inf in `se`/`auc`** → **Guard:** emit NA CI, never crash; row still written.
6. **Append position** → new columns must be **last**; existing column indices unchanged (downstream
   positional readers, the oracle comparator, the Docker self-test).

**Part B (shape confidence):**
7. **Weak cluster `sum(w)=0`** → `focus=NA`, division by zero. **Guard:** branch on `sum(w)>0`; focus
   CI = NA; `shape_confidence` still computed.
8. **Single significant domain `n_sig=1`** → `focus=1`, sharp by rule. MC mostly sharp; no degenerate
   quantile.
9. **One category at level (`length(auc)==1`)** → `margin=NA`. **Guard:** propagate NA margin into
   `decide_shape` per draw (sharp via `n_sig==1`).
10. **Few significant draws** (< `shape_confidence_min_sig_frac·B`) → unstable focus quantiles.
    **Guard:** report `anchor_focus_ci_*` = NA.
11. **Tie at the margin sort (`astar[1]==astar[2]`)** → `margin*=0`; fine, no special case but must not
    error.
12. **`se=NA` for a category** → **Guard:** that category's AUC is held fixed (no perturbation) in draws.
13. **RNG kind drift** under the parallel sweep (`future.seed=TRUE` flips to L'Ecuyer) → **Guard:** the
    per-cluster `set.seed(..., kind="Mersenne-Twister", normal.kind="Inversion")` pin, exactly as
    [R/sensitivity.R:56-60](../../R/sensitivity.R#L56-L60).

## VALIDATION COMMANDS (run all; zero failures)

```bash
# 1. Unit + edge + reproducibility suite
Rscript -e 'testthat::test_local()'

# 2. Smoke run on the synthetic .rds fixture (writes the augmented TSVs)
Rscript inst/scripts/anchor_map.R --config synthetic_rds --out-dir results/synthetic_rds --threads 4

# 3. Schema/containment check on the augmented scores
Rscript -e '
  d <- data.table::fread("results/synthetic_rds/category_anchor_scores.tsv")
  stopifnot(all(c("auc_abs_se","auc_abs_ci_lo","auc_abs_ci_hi") %in% names(d)))
  stopifnot(all(d$auc_abs_ci_lo <= d$auc_abs & d$auc_abs <= d$auc_abs_ci_hi))
  stopifnot(all(d$auc_abs_ci_lo >= 0 & d$auc_abs_ci_hi <= 1))
  l <- data.table::fread("results/synthetic_rds/cluster_anchor_labels.tsv")
  stopifnot(all(c("shape_confidence","anchor_focus_ci_lo","anchor_focus_ci_hi",
                  "shape_posterior","shape_jackknife_stable") %in% names(l)))
  stopifnot(all(l$shape_confidence >= 0 & l$shape_confidence <= 1))
  cat("schema + containment OK\n")'

# 4. Cross-language parity unchanged on the disease track (extra R columns are ignored by name-driven
#    comparator) — every Phase-1 deterministic column still matches the Python oracle.
Rscript validation/compare_oracle.R \
  --r-out results/<disease_run>/category_anchor_scores.tsv \
  --oracle "../UKBB_CLUSTER_GWAS/scripts/cluster_anchoring/output/carey_rint15/category_anchor_scores.tsv"

# 5. Legacy contract escape hatch (emit_uncertainty:false reproduces the 24/11-col Phase-5 contract)
Rscript inst/scripts/anchor_map.R --config synthetic_rds --out-dir results/synthetic_legacy --threads 4  # with emit_uncertainty:false
```

## ACCEPTANCE CRITERIA

- [ ] `category_anchor_scores.tsv` gains `auc_abs_se/ci_lo/ci_hi` (appended); `cluster_anchor_labels.tsv`
      gains `shape_confidence`, `anchor_focus_ci_lo/hi`, `shape_posterior`, `shape_jackknife_stable`.
- [ ] **Every existing column value is byte-identical** to Phase 5 (AUC, ranks, `pooled_rg*`, `vif*`,
      `perm_p`, `q`, labels, shapes, `anchor_focus`) — Phase 6 is purely additive.
- [ ] `ci_lo ≤ auc_abs ≤ ci_hi`, all ∈ [0,1], for every row; VIF widens the interval.
- [ ] DeLong CI matches `pROC::ci.auc(method="delong")` (VIF=1) within 1e-6 on clean categories.
- [ ] Positive control C5_sub0: wide AUC CI (width > 0.1), `shape_confidence ≈ 1.0`, sharp,
      `shape_jackknife_stable = True`. Negative control C3: weak, focus CI = NA.
- [ ] All 13 edge cases guarded with a test exercising each degenerate path.
- [ ] `shape_confidence` thread- and order-invariant; cross-language parity (deterministic cols)
      unchanged; `emit_uncertainty:false` reproduces the Phase-5 contracts byte-for-byte.
- [ ] New deps confined to **tests** (`pROC` as an optional cross-check, `Suggests:` not `Imports:`);
      the runtime stays dependency-clean (pure base-R math).
- [ ] `ANALYSIS_DESIGN.md` Phase-6 section + `README.production.md` updated; configs carry the new block
      with defaults.

## NOTES

- **No new runtime dependency.** DeLong/Hanley–McNeil/logit are base-R arithmetic; `pROC` is a
  **test-only** cross-check (`Suggests`). This protects the pinned Docker image and the reproducibility
  headline.
- **Determinism split.** Part A is fully deterministic (no RNG) → parity-safe, always on. Part B's MC is
  reproducible-within-seed but Monte-Carlo across languages — treat `shape_confidence`/`shape_posterior`
  like `perm_p`/`q`: distributional (MC tolerance), never a byte-parity gate. The jackknife is
  deterministic.
- **Two complementary confidence axes.** `shape_confidence` = within-z statistical robustness (rg noise);
  the existing `label_stable` ([R/sensitivity.R](../../R/sensitivity.R)) = across-z threshold robustness.
  Report both; do not merge.
- **Scope boundary.** These intervals capture rg sampling noise given `rg_se`/VIF only — not ontology
  error, misspecification, or selection. State this in the README so readers don't over-read a tight CI
  as biological certainty.
- **Why MC, not a closed form, for shape.** The verdict is a discrete function crossing `n_sig`,
  `margin`, and `focus` thresholds simultaneously; a delta-method CI on `focus` alone cannot represent
  the threshold flips. MC propagation does, and naturally yields the full `shape_posterior`.
- **Optional future extension** (defer): a `bootstrap` `auc_ci_method` (resample traits within
  in/out groups) for the small-n regime — assumption-free but stochastic and degenerate at n_in=3, so
  not the default; and perturbing significance (re-permutation) inside the shape MC.
```
