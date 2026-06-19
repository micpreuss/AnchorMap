# Analysis: AnchorMap Phase 4 — Visualization (R/ggplot2 port of the reference figures)

Validate the schema contracts and method usage against the actual code/data before implementing.
Pay special attention to column names, the `eligible` string encoding, the `level` value per track,
and the ordering algorithms (hierarchical leaf order + natural cluster order).

## Question & object
**Question.** Given the already-scored anchoring tables (`category_anchor_scores.tsv`,
`cluster_anchor_labels.tsv`), render **publication-ready figures** that make each cluster's anchor
profile and cross-cluster distinctiveness legible — *without* recomputing any statistic. The figures
must keep **AUC** (sign-blind magnitude / ranking) and **pooled_rg** (signed direction) as *distinct
visual channels*, because they diverge exactly at sign-split classes (the load-bearing insight).

**Object.** Two TSVs per track (scores: one row per `(cluster, level, category)`; labels: one row per
cluster), read into long `data.table`s; derived **pivot matrices** (cluster × category of `pooled_rg`,
`q`) for ordering and the specificity z-transform. Output objects are **figures** (PNG + PDF) plus one
small `cluster_distinctive_categories.tsv` side-table.

**Inference licensed.** Visual reading of (a) *what each cluster is* (anchoring: lollipop / dot-heatmap
/ AUC-vs-coherence) and (b) *what makes each cluster different from its siblings* (specificity heatmap
+ diagonal). No new estimands — this is a faithful re-encoding of Phase 1–3 outputs.

## Analyst story
As an analyst / I want **headless, config-driven figures rendered from the scored AnchorMap TSVs that
port the reference `plot_anchors.py` / `plot_specificity*.py` encodings to R** / So that **I can read
each cluster's anchor (shape, direction, significance) and its cross-cluster specificity from
publication-ready plots, with the AUC-vs-pooled_rg divergence preserved**.

## Pipeline position
- **Upstream:** `anchor_map.R` (Phases 1–3) → `results/<run>/category_anchor_scores.tsv` +
  `cluster_anchor_labels.tsv` (and the sensitivity TSVs, not consumed here).
- **Downstream:** human reporting / `report-findings`; no programmatic consumer. The distinctive TSV is
  a convenience side-table.
- **Orchestration pattern to mirror:** the **reference is three standalone Python scripts** that each
  read the scored tables + a shared plot config and render to `figures/`. We unify them into **one R
  module `R/plot.R` + one CLI `plot_anchors.R --config <plots.yaml>`** (mirrors the engine's
  `anchor_map.R` entry: [anchor_map.R:18-41](anchor_map.R#L18-L41) script-dir + tiny arg parser, and
  [R/io.R:34-51](R/io.R#L34-L51) `load_config`/`stage_root_of`/`resolve_path`).

---

## CONTEXT REFERENCES — READ BEFORE IMPLEMENTING

### Input schemas (verified from real headers + the writer)
- **`results/<run>/category_anchor_scores.tsv`** — columns (exact, from
  [anchor_map.R:44-47](anchor_map.R#L44-L47)):
  `cluster_label, level, category, eligible, n, n_eff, n_hit, rho_bar, vif, auc_abs, auc_signed,
  perm_p, vif_z, vif_p, pooled_rg, pooled_rg_ci_lo, pooled_rg_ci_hi, coherence, mean_abs_rg,
  mean_signed_rg, odds_ratio, fisher_p, q, rank`.
  - **`eligible` is the string `"True"`/`"False"`** (pandas bool repr, written at
    [anchor_map.R:126](anchor_map.R#L126)). Python `read_csv` re-infers it to bool; **R `fread` reads it
    as character → you MUST convert `eligible == "True"`** before filtering. *This is the #1 silent port
    bug for the plots.*
  - `level` partitions the file: disease track uses **`domain`** (also has `native`, `icd_chapter`
    rows — filter them out); anthro track uses **`anthro_class`** (single category `Anthropometric`);
    lab track would use `analyte_class`. Verified: `cut -f2` on the real files.
  - `auc_abs` ∈ [0,1] is the **x-axis / size** channel (sign-blind). `pooled_rg` ∈ [−1,1] is the
    **colour** channel (signed). `coherence` ∈ [0,1] (NA→1.0, clip [0,1]) is **alpha**. `q` < 0.05 is
    the significance ring/mask. These four are the encoding; do not collapse AUC and pooled_rg.
- **`results/<run>/cluster_anchor_labels.tsv`** — columns
  `cluster_label, auto_label, anchor_shape, anchor_margin, anchor_focus, n_sig_domains, top_auc, top_q,
  top_pooled_rg, top_coherence, profile`. Used for panel titles (`cluster — auto_label [shape]`) and the
  ★ auto-label marker (match `auto_label == category`).

### Output schema (contracts produced)
- `results/<run>/figures/anchor_lollipops_<track>.{png,pdf}` — one per track.
- `results/<run>/figures/anchor_dotheatmap.{png,pdf}` — one combined (all tracks side by side).
- `results/<run>/figures/anchor_auc_coherence.{png,pdf}` — one combined diagnostic scatter.
- `results/<run>/figures/anchor_specificity_<track>.{png,pdf}` — one per track.
- `results/<run>/figures/anchor_specificity_diagonal_<track>.{png,pdf}` — one per track (skipped if no
  significant distinctive cell).
- `results/<run>/figures/cluster_distinctive_categories.tsv` — columns
  `track, cluster_label, distinctive_category, spec_z, pooled_rg, runner_up` (from
  `plot_specificity.py:distinctive_table`, lines 70-85; `spec_z` round 2, `pooled_rg` round 3).

### Reference data & methods
- **Reference figures (bit-for-bit encoding spec, NOT bit-for-bit pixels):**
  - [plot_anchors.py](../UKBB_CLUSTER_GWAS/scripts/cluster_anchoring/plot_anchors.py) — lollipop
    small-multiples (L79-128), dot-heatmap (L132-176), AUC-vs-coherence scatter (L180-208), shared
    colorbar+legend (L212-228).
  - [plot_specificity.py](../UKBB_CLUSTER_GWAS/scripts/cluster_anchoring/plot_specificity.py) —
    within-category z of `pooled_rg` (L60-67), distinctive table (L70-85), masked heatmap (L88-119).
  - [plot_specificity_diagonal.py](../UKBB_CLUSTER_GWAS/scripts/cluster_anchoring/plot_specificity_diagonal.py)
    — boxed assignment (L37-50), greedy diagonal column order (L53-59), reduced heatmap (L62-99).
  - Reference rendered outputs to eyeball against:
    `../UKBB_CLUSTER_GWAS/scripts/cluster_anchoring/output/carey_rint15/figures/*.png`.
  - Encoding rationale: `../UKBB_CLUSTER_GWAS/scripts/cluster_anchoring/docs/figures_guide.md` (§1-§4)
    and `docs/approach.md` §4.4/§4.7/§5.2.
- **Methods to port (no scipy/matplotlib in R):**
  - **Diverging colour, centered at 0:** matplotlib `TwoSlopeNorm(-cap,0,+cap)` + cmap `RdBu_r` (rg) /
    `PuOr_r` (specificity). R: `ggplot2::scale_*_gradient2(low=<blue>, mid="white", high=<red>,
    midpoint=0, limits=c(-cap,cap), oob=scales::squish)`. RdBu_r endpoints ≈ low `#2166AC` / high
    `#B2182B`; PuOr_r ≈ low `#542788` (purple) / high `#B35806` (orange). `rg_cap` default **0.55**;
    specificity `cap` default **2.5**.
  - **Category (column) order = hierarchical-clustering leaf order** (`plot_anchors.leaf_order`,
    L44-49): on the `category × cluster` (lollipop/dotheatmap) or `cluster × category` (specificity, on
    `Z.T`) pivot of `pooled_rg`/`Z`, `NaN→0`, `linkage(method="average", metric="euclidean")`,
    `leaves_list`. R equivalent: `stats::hclust(dist(M, method="euclidean"), method="average")$order`
    over `M` with `NA→0`; if `nrow(M) < 3` keep input order.
  - **Cluster (row) order = natural order** (`natural_order`, L52-62): regex
    `^C(\d+)(?:_sub(\d+))?$` then `^noise_re(\d+)$` then everything else; sub index −1 sorts the bare
    `C5` before `C5_sub0`. Port the comparator exactly (the C0…, sub-sorted, `noise_re*`-last order is
    used across all figures).
  - **Specificity z** (`specificity`, L60-67): `M = pivot(cluster × category, pooled_rg)`;
    `Z = (M − colMean) / colSD(ddof=0)` (population SD, **divide by n not n−1** — use
    `sd_pop = sqrt(mean((x-mean)^2))`); `mask = (q < q_sig) & (|M| ≥ rg_floor) & (n_present_per_col ≥
    min_clusters)`; `n_present` counts non-NA cluster entries per category column. Defaults
    `q_sig=0.05, rg_floor=0.10, min_clusters=5`.
  - **Size ∝ AUC** (dot-heatmap/scatter): matplotlib point area
    `s = (clip(auc−0.5,0,0.5)/0.5)*200 + 12` (+120 for scatter). In ggplot map a derived
    `auc_size = clip(auc_abs−0.5,0,0.5)` via `scale_size_area()` or `scale_radius()` — **monotone in AUC
    is the contract, exact px is not.** Legend breaks at AUC = 0.6 / 0.75 / 0.9.
  - **Lollipop** (L87-115): per cluster, top-`k` rows by `rank`; horizontal segment from x=0.5 to
    `auc_abs`; colour=`pooled_rg`, alpha=`0.30+0.70*coherence`, black ring if `q<q_sig`, open ★ if
    `category==auto_label`; category text above stem, `n<n>` right of dot; vline at 0.5; xlim [0.5,1.0].
- **Method assumptions / failure modes:** figures are deterministic given the TSVs (no RNG). Failure
  modes: a track with <3 categories (skip hclust ordering, keep input order — the anthro track has a
  **single** `Anthropometric` category, so leaf_order must no-op gracefully); a cluster with zero
  eligible rows at the level (drop panel); specificity with no significant cell (diagonal returns
  nothing → skip with a log line, mirroring `plot_diagonal` L67-68); all-NA colour column.

### Files to read / create
- READ: [anchor_map.R](anchor_map.R) (L18-41 entry/argparse; L43-52 column contracts) — Why: mirror the
  CLI + script-dir sourcing pattern.
- READ: [R/io.R](R/io.R) (L34-51 `load_config`/`stage_root_of`/`resolve_path`; L77-84 ontology reader) —
  Why: reuse config loader + path resolution; the plot config is plain YAML like the engine configs.
- READ: [R/sensitivity.R](R/sensitivity.R) — Why: mirror the module header/comment density + the
  `parallel`-free helper style for the new `R/plot.R`.
- READ: `../UKBB_CLUSTER_GWAS/scripts/cluster_anchoring/configs/carey_rint15_plots.yaml` — Why: the
  canonical plot-config schema (`out_dir, top_k, lollipop_ncols, scatter, rg_cap, row_order_track,
  tracks:[{name, level, scores, labels}]`).
- CREATE: `R/plot.R` — all figure builders + helpers (leaf_order, natural_order, diverging scale,
  specificity z, distinctive table), pure functions returning ggplot/patchwork objects + a `save_fig`.
- CREATE: `plot_anchors.R` — CLI entry: `--config <plots.yaml> [--q-sig --rg-floor --min-clusters]`;
  loads tracks, renders all five figure families + the distinctive TSV; logs `[write]` lines.
- CREATE: `configs/carey_rint15_plots.yaml` — AnchorMap plot config pointing at `results/…` (disease +
  anthro tracks that exist locally; lab track stub commented out — no lab result is committed).

---

## METHOD / IMPLEMENTATION PLAN

### Phase 1: Inputs & harmonization
- Reuse `load_config`/`stage_root_of`/`resolve_path` from `io.R`. Add plot-config defaults
  (`top_k=8, lollipop_ncols=3, scatter=TRUE, rg_cap=0.55, q_sig=0.05`).
- `load_track(t, stage_root)`: `fread` the scores TSV; **`eligible <- eligible == "True"`**; filter
  `level == t$level & eligible`; `coherence <- pmin(pmax(fifelse(is.na(coherence),1,coherence),0),1)`.
  Read labels TSV. Return `list(name, level, s, labels)`.
- Build per-track pivots with `data.table::dcast` (category × cluster of `pooled_rg`) for ordering.

### Phase 2: Core analysis (ordering + transforms; no new statistics)
- `natural_order(labels)` — exact regex comparator port (C…, sub-sorted, noise_re* last, others last).
- `leaf_order(M)` — `NA→0`; if `nrow<3` return `rownames`; else
  `rownames(M)[hclust(dist(M,"euclidean"),"average")$order]`.
- `specificity(s, q_sig, rg_floor, min_clusters)` — pivots M (pooled_rg), Q (q); population-SD column z;
  significance mask incl. `n_present_per_col ≥ min_clusters`. Return `list(M, Z, mask)`.
- `distinctive_table(M, Z, mask, track)` and `diagonal_column_order(assignment)` — direct ports.

### Phase 3: Figures & integration
- `make_rg_scale()/make_spec_scale()` — `scale_*_gradient2` diverging, oob squish, fixed limits.
- `fig_lollipops(track, row_order, cfg)` — one ggplot per cluster (`geom_segment` + `geom_point` +
  `geom_text` for category & `n`; ★ via `geom_point(shape=8/"★")`; ring via `geom_point` with
  `stroke`/black `colour` when `q<q_sig`); assemble with `patchwork::wrap_plots(ncol=lollipop_ncols)` +
  a `plot_annotation` suptitle. Panel title `cl — auto_label [shape]`.
- `fig_dotheatmap(tracks, row_order, cfg)` — one `geom_point` panel per track
  (`aes(x=category, y=cluster, size=auc_size, colour=pooled_rg, alpha=coherence)`, black `stroke` when
  `q<q_sig`, ★ overlay at `(auto_label, cluster)`); `patchwork` side-by-side with shared rg scale; y
  reversed (C0 top).
- `fig_scatter(tracks, cfg)` — `aes(x=auc_abs, y=coherence, size=auc_size, colour=pooled_rg)`; shade
  `coherence<0.5` band; dashed line at 0.6; text-label points with `auc_abs≥0.65 & coherence≤0.6`
  (`ggrepel` if available, else `geom_text`); facet/patchwork per track.
- `fig_specificity(spec, row_order, name, cap)` — `geom_tile(aes(fill=Z))` with mask→grey
  (`fill=NA`/`na.value="grey90"`), `geom_rect`/`geom_tile(colour="black")` box on each cluster's max-|z|
  significant cell; PuOr scale; columns in `leaf_order(Z.T)`, rows in natural order.
- `fig_diagonal(spec, row_order, name, cap)` — reduce to one boxed cell per cluster, greedy diagonal
  column order; skip (log `[skip]`) if no significant cell.
- `save_fig(plot, path, w, h)` — write **PNG + PDF**. Use `ragg::agg_png` if available, else
  `grDevices::png(type="cairo")` (or `ggsave(device="png")`); PDF via `ggsave(device=cairo_pdf)`.
  **Headless:** never require X11; set `options(bitmapType="cairo")` guard.
- CLI `plot_anchors.R`: parse args, load tracks, compute `row_order = natural_order(∪ clusters)`, render
  all figures + write `cluster_distinctive_categories.tsv`, emit `[write] <path>` per artifact.

### Phase 4: Validation
- Run the CLI on `configs/carey_rint15_plots.yaml`; assert all expected files exist and are non-empty;
  eyeball against the reference PNGs for the same run.

---

## STEP-BY-STEP TASKS (execute top to bottom; each atomic + checkable)

### CREATE configs/carey_rint15_plots.yaml
- **IMPLEMENT**: `out_dir: results/carey_rint15/figures`; `top_k: 8`; `lollipop_ncols: 3`;
  `scatter: true`; `rg_cap: 0.55`; `row_order_track: disease`; `tracks:` → disease
  (`level: domain`, `scores: results/carey_rint15/category_anchor_scores.tsv`, labels alongside) and
  anthro (`level: anthro_class`, `scores: results/carey_rint15_anthro/…`). Lab track commented
  (no committed result).
- **PATTERN**: reference `configs/carey_rint15_plots.yaml`.
- **DATA/SCHEMA**: paths relative to stage root (AnchorMap repo root); `level` values verified above.
- **GOTCHA**: AnchorMap writes flat `results/<run>/…`, not the reference `output/<run>/figures/` — point
  paths at the real local TSVs.
- **VALIDATE**: `Rscript -e 'yaml::read_yaml("configs/carey_rint15_plots.yaml")'`.

### CREATE R/plot.R (helpers + figure builders)
- **IMPLEMENT**: `natural_order`, `leaf_order`, `load_track`, `specificity`, `distinctive_table`,
  `diagonal_column_order`, `make_rg_scale`, `make_spec_scale`, `fig_lollipops`, `fig_dotheatmap`,
  `fig_scatter`, `fig_specificity`, `fig_diagonal`, `save_fig`. Pure functions; `suppressPackageStartupMessages`
  for `ggplot2`, `patchwork`, `scales`, `data.table`.
- **PATTERN**: module header + comment density of [R/sensitivity.R](R/sensitivity.R); config/path helpers
  from [R/io.R](R/io.R).
- **DATA/SCHEMA**: consumes the score/label contracts above; the only mutation is
  `eligible == "True"` + coherence clamp.
- **GOTCHA**: population SD (ddof=0) in `specificity`; `leaf_order` no-op for <3 rows (anthro single
  category); `pooled_rg` colour limits fixed to ±`rg_cap` with `oob=scales::squish` so out-of-range
  doesn't drop to grey.
- **VALIDATE**: `Rscript -e 'source("R/plot.R"); stopifnot(identical(natural_order(c("noise_re0","C5_sub1","C0","C5_sub0","C2")), c("C0","C2","C5_sub0","C5_sub1","noise_re0")))'`.

### CREATE plot_anchors.R (CLI)
- **IMPLEMENT**: script-dir sourcing of `R/io.R` + `R/plot.R`; tiny arg parser
  (`--config`, `--q-sig`, `--rg-floor`, `--min-clusters`); load all tracks; `row_order`; render lollipop
  per track, one dot-heatmap, one scatter (if `cfg$scatter`), specificity + diagonal per track; write the
  distinctive TSV; print `[write]` per file.
- **PATTERN**: [anchor_map.R:18-41,164-167](anchor_map.R#L18-L41).
- **DATA/SCHEMA**: out files per the Output schema above.
- **GOTCHA**: create `out_dir` recursively; skip diagonal gracefully when empty; PDF + PNG for every
  figure except the scatter (reference writes scatter PNG only — match it, or add PDF and note the
  deviation).
- **VALIDATE**: `Rscript plot_anchors.R --config configs/carey_rint15_plots.yaml`.

### ADD ragg/ggrepel to the designed dependency list (docs only this phase)
- **IMPLEMENT**: note in `ANALYSIS_DESIGN.md` §7.3 R-deps already lists `ggplot2/patchwork/ragg/scales`;
  add `ggrepel` (optional, for scatter labels) so Phase 5's Dockerfile pins it.
- **GOTCHA**: `ragg` is **not** installed in the local R right now (verified) — code must fall back to
  `cairo` so local validation works without it; the container pins `ragg`.
- **VALIDATE**: `Rscript -e 'cat(requireNamespace("ragg",quietly=TRUE))'` (FALSE locally is fine).

### UPDATE CLAUDE.md + README (wire the stage in)
- **IMPLEMENT**: flip Phase 4 from "designed-not-built" to built in [CLAUDE.md](CLAUDE.md) status
  paragraph + repo-structure (`R/plot.R`, `plot_anchors.R`, `configs/*_plots.yaml`); add a "Running
  stages" line. Use `create-readme` for the per-stage README.
- **VALIDATE**: `grep -n "plot.R" CLAUDE.md`.

## VALIDATION STRATEGY
No unit-test suite for figures; validate by:
- **Smoke run** of `plot_anchors.R` on the committed `carey_rint15` + `carey_rint15_anthro` results.
- **Positive controls (visual, against the reference PNGs):**
  - **Anthro lollipop**: `C5_sub0` panel shows the single `Anthropometric` lollipop **red** (positive
    `pooled_rg≈0.25`), **black ring** (`q<0.05`), **open ★** (auto-label) — its adiposity core.
  - **Disease dot-heatmap**: `C0` row's largest dot at `Musculoskeletal`/`Psychiatric` with a ★ on the
    auto-label cell; `C2_sub0`/`C2_sub1` show **blue** (negative `pooled_rg`) dots at **high AUC** — the
    AUC↔pooled_rg divergence that proves the two channels are not redundant.
  - **AUC-vs-coherence scatter**: at least one high-AUC (≥0.65) low-coherence (≤0.6) point is labelled
    (sign-split class) — the diagnostic's whole purpose.
  - **Specificity diagonal (disease)**: boxes run down a diagonal; each cluster's most-distinctive domain
    is a single boxed cell.
- **Negative / robustness control:** a diffuse cluster (`C5_sub0` on the **disease** track) stands out on
  *no* specificity cell (grey row) even though it anchors on the anthro track — the anchoring-vs-
  specificity distinction.
- **Schema/sanity checks:** every expected PNG **and** PDF exists and is >1 KB; `cluster_distinctive_categories.tsv`
  has the 6-column contract; column ordering for category axes is the hierarchical leaf order; cluster
  rows are natural order; no figure errors on the single-category anthro track.
- **Ordering port check:** `natural_order` and `leaf_order` unit assertions (above) — the only nontrivial
  ported logic.

## VALIDATION COMMANDS (run all; zero schema/control failures)
```bash
# config parses
Rscript -e 'yaml::read_yaml("configs/carey_rint15_plots.yaml")'

# ordering ports
Rscript -e 'source("R/plot.R"); stopifnot(identical(natural_order(c("noise_re0","C5_sub1","C0","C5_sub0","C2")), c("C0","C2","C5_sub0","C5_sub1","noise_re0"))); cat("natural_order OK\n")'

# smoke render
Rscript plot_anchors.R --config configs/carey_rint15_plots.yaml

# all expected artifacts exist and are non-empty
ls -l results/carey_rint15/figures/{anchor_lollipops_disease,anchor_lollipops_anthro,anchor_dotheatmap,anchor_auc_coherence,anchor_specificity_disease,anchor_specificity_anthro,anchor_specificity_diagonal_disease}.png \
      results/carey_rint15/figures/cluster_distinctive_categories.tsv
Rscript -e 'f<-Sys.glob("results/carey_rint15/figures/*.png"); stopifnot(length(f)>=5, all(file.size(f)>1024)); cat(length(f),"PNGs OK\n")'

# distinctive TSV contract
Rscript -e 'd<-data.table::fread("results/carey_rint15/figures/cluster_distinctive_categories.tsv"); stopifnot(all(c("track","cluster_label","distinctive_category","spec_z","pooled_rg","runner_up") %in% names(d))); cat("distinctive TSV OK\n")'

# eyeball vs the reference renders
echo "compare: results/carey_rint15/figures/*.png  vs  ../UKBB_CLUSTER_GWAS/scripts/cluster_anchoring/output/carey_rint15/figures/*.png"
```

## ACCEPTANCE CRITERIA
- [ ] `eligible` parsed as `=="True"` (not string-truthiness); only `level`-matching eligible rows plotted.
- [ ] Five figure families render headless (PNG+PDF; scatter per reference) for the disease + anthro
      tracks without X11.
- [ ] AUC and pooled_rg are **distinct channels** — the divergence is visible (C2_sub* blue at high AUC;
      a labelled sign-split point in the scatter).
- [ ] Positive controls recovered (C5_sub0 anthro red+ring+★; C0 disease star on auto-label).
- [ ] Negative/robustness control holds (C5_sub0 grey on disease specificity).
- [ ] Category axes in hierarchical leaf order; cluster rows in natural order; single-category anthro
      track renders without error.
- [ ] `cluster_distinctive_categories.tsv` matches the 6-column contract; rounding (spec_z 2, pooled_rg 3).
- [ ] Stage wired into CLAUDE.md / ADD / README; `plot_anchors.R` usage documented.
- [ ] Provenance: config echoed; `[write]` manifest printed (mirror the engine's log style).

## NOTES
- **Not bit-for-bit.** The Phase-4 gate (ADD §12) is *encoding fidelity + headless render*, not pixel
  parity — matplotlib→ggplot point-area math differs, so size/alpha are matched **monotonically**, not
  exactly. Keep AUC (x-position/size), pooled_rg (diverging colour), coherence (alpha), q (ring/mask) as
  the four channels; that is the contract.
- **One CLI vs three scripts.** The reference splits into `plot_anchors.py` + `plot_specificity.py` +
  `plot_specificity_diagonal.py`. Unifying into `plot_anchors.R` is intentional (slimmer, single config
  read); the `R/plot.R` functions stay separable so a future split is trivial. If you prefer parity with
  the reference invocation surface, expose `--only lollipop|dotheatmap|scatter|specificity|diagonal`.
- **`ragg` optional locally.** Verified absent in the current R; the `save_fig` cairo fallback keeps
  local validation green. Phase 5 pins `ragg` in the Dockerfile for crisp raster output.
- **Lab track deferred** in the local config (no committed `results/carey_rint15_lab/`). The config keeps
  a commented lab stub so it activates the moment a lab run exists — matches the reference 3-track config.
- **scales/colours:** use explicit RdBu_r (rg) and PuOr_r (specificity) endpoints so the two heatmap
  families stay visually distinct (the reference deliberately uses different cmaps — `plot_specificity.py`
  L33 comment).
- Do **not** branch (CLAUDE.md git convention) — commit on the current branch via `commit`.
```
