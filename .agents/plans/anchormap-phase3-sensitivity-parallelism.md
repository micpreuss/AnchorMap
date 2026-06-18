# Analysis: AnchorMap Phase 3 — Reliability-threshold sensitivity sweep + parallelism

Validate the contracts and the RNG/parallel behaviour against the actual code before implementing.
Pay special attention to: **per-z RNG seeding** (the single load-bearing correctness point — it must make
results *invariant to thread count* AND make the primary-z slice reproduce the committed Phase-1/2 output
byte-for-byte), the **gate-per-z re-run** semantics (z changes the universe N), and **not perturbing
`perm_p`** (parallelize over z only, never the permutation draws).

> Scope: **Phase 3 of `ANALYSIS_DESIGN.md`** only — Aim **A4 (first half)**: a parallel z-threshold
> sensitivity sweep emitting two extra TSVs + a stability flag, with thread-count-invariant results.
> Build on the Phase-1/2 R engine (`R/{io,gate,redundancy,score,label,ingest_rds}.R` + `anchor_map.R`),
> validated bit-for-bit against the Python oracle and (Phase 2) the synthetic `.rds` battery.
> **Out of Phase 3 (do NOT build here):** the plotting module (Phase 4), the Docker image + Nextflow
> `ANCHORMAP` process + the `results/<run_label>/{primary,sensitivity,figures,logs}/` subdir reorg
> (Phase 5). **Invariant to preserve:** the two *primary* TSVs (`category_anchor_scores.tsv`,
> `cluster_anchor_labels.tsv`) must stay byte-identical to today for the existing configs and the
> `.rds` route — Phase 3 only *adds* two `sensitivity_z_*.tsv` files and parallelizes the sweep.

## Question & object
**Scientific question.** Is each cluster's `auto_label` / `anchor_shape` an artefact of the (arbitrary)
h²-reliability cut `z`, or is it *robust* across a bracket of reasonable cuts? The h²-reliability gate
`h2_z = h2_trait / h2_trait_se > z` defines which traits are trustworthy enough to enter the competitive
universe; **z changes the universe N**, so every z is a full independent re-run of the whole engine
(gate → redundancy → AUC/perm/VIF/pooled-rg/ORA → FDR → label/shape). The object operated on is the
same cluster×trait `rg` long-table; the new representation is a **stacked, z-indexed** version of the two
output tables plus a per-cluster **`label_stable`** flag.

**Inference licensed.** "Cluster C's anchor is/ isn't stable to the reliability threshold" — the sweep
TSV says *at which z values* a label flips, turning the single-z auto-label into an auditable robustness
claim. (ADD §1 "robust across reliability thresholds"; §6 step 1 note; §8 "Sensitivity check".)

## Analyst story
As an analyst / I want AnchorMap to re-score every cluster across a vector of h²-reliability thresholds
**in parallel** and flag whether each cluster's auto-label is stable / So that I can report a defensible,
reproducible anchor *with* its robustness, instead of trusting a single arbitrary z — and get identical
numbers no matter how many CPUs the run used.

## Pipeline position
- **Upstream:** the Phase-1/2 ingestion + gate chain — `read_long` / `read_rds_route` → `df`
  ([anchor_map.R:67-81](../../anchor_map.R#L67-L81)). Unchanged.
- **Core (re-used per z):** `apply_universe_gate` → `attach_ontology` → `select_corr_source` →
  `score_cluster_level` (×clusters×levels) → `rank_and_label`
  ([anchor_map.R:83-104](../../anchor_map.R#L83-L104)). Phase 3 wraps this in a per-z function and maps it
  over a z-vector.
- **Downstream:** Phase 4 (`R/plot.R`) consumes `category_anchor_scores.tsv` (primary) and may overlay the
  `sensitivity_z_*` tables; Phase 5 (Nextflow `ANCHORMAP`) wires `--threads ${task.cpus}` + the
  `output {}`/`outputDir` publish. Neither built here.
- **Orchestration pattern to mirror:** [anchor_map.R](../../anchor_map.R) `run_anchormap()` (the
  `emit()`/`[start]/[config]/[gate]/[vif]/[write]/FINISHED` log lines; `setDTthreads`; the
  `.SCORE_COLS`/`.LABEL_COLS` write contract). Test harness: [tests/run_tests.R](../../tests/run_tests.R)
  + [tests/test_phase2.R](../../tests/test_phase2.R) (plain `stopifnot`, **no `testthat`** locally).

---

## CONTEXT REFERENCES — READ BEFORE IMPLEMENTING

### Input schemas (unchanged from Phase 1/2)
- rg long-table (Input A) — `read_long` ([io.R:55](../../R/io.R#L55)); GenomicSEM `.rds` (Input C) —
  `read_rds_route` ([R/ingest_rds.R](../../R/ingest_rds.R)). Phase 3 adds **no new input file**; it adds
  one config key (`z_vector`) and one CLI flag (`--z-vector`).
- The gate's z is `cfg$h2_z_threshold` ([gate.R:18](../../R/gate.R#L18)); Phase 3 makes that an
  **overridable argument** so the sweep can vary it without mutating `cfg`.

### Output schema (contract for downstream)
Two **new** TSVs written flat into `out_dir` (alongside — not replacing — the two primary TSVs):

- `sensitivity_z_scores.tsv` — columns `= c(.SCORE_COLS, "z_threshold")` i.e. the full
  [anchor_map.R:42-45](../../anchor_map.R#L42-L45) score contract with `z_threshold` **appended last**;
  one block of score rows per swept z, stacked (sorted `z_threshold ↑, level, cluster_label, rank`).
  `eligible` printed `True/False` exactly as the primary table.
- `sensitivity_z_labels.tsv` — columns `= c(.LABEL_COLS, "z_threshold", "label_stable")`
  ([anchor_map.R:46-47](../../anchor_map.R#L46-L47)) with `z_threshold` and `label_stable` appended;
  one labels block per swept z. **`label_stable`** is per *cluster* (broadcast onto each of that cluster's
  z-rows): `True` iff `auto_label` is identical across **every** swept z, else `False`; printed
  `True/False` (pandas-style repr, matching `eligible`).

**Primary tables unchanged.** `category_anchor_scores.tsv` / `cluster_anchor_labels.tsv` keep their exact
paths, columns, ordering, rounding, and values — they are the **primary-z slice** of the sweep (see the
parity construction below). The `results/<run_label>/{primary,sensitivity,figures,logs}/` subdir layout in
ADD §7.5 is **deferred to Phase 5** (moving the primary TSVs now would break the oracle-comparison path and
the Phase-1/2 committed results).

### Reference data & methods
- **No Python oracle for the sweep** — `anchor_categories.py:main()` (L479-543) is single-z, serial,
  single-track; it has **no** `z_vector`/parallel/`sensitivity` code (verified by grep). Phase 3 is a
  clean AnchorMap generalization, so the parity anchor is **internal**: the primary-z slice must equal the
  existing committed primary output (which itself already passed the cross-language oracle in Phase 1).
- **Parallel backend (decision — declared dependency, not availability-contingent).** Use
  **`future` + `future.apply`** as named by ADD §7.4 / CLAUDE.md (`future_lapply` over z;
  `future::plan("multicore", workers = threads)`). `future.apply` is the right tool on merit — backend-
  agnostic via `plan()` (multicore/multisession/cluster, incl. the Phase-5 container), parallel-safe
  L'Ecuyer RNG (`future.seed = TRUE`), and a clean `lapply` drop-in — strictly better than hand-rolling
  `parallel::mclapply` (unix-fork-only, manual seeding). It is now **installed locally**
  (`future` 1.x, `future.apply` 1.20.2; verified `future_lapply` runs under `multicore` + reseed-
  determinism holds) and will be pinned in the Phase-5 Docker image (dated P3M snapshot). **Do not gate the
  code on what is installed** — `future.apply` is a hard dependency; `parallel::mclapply` is *not* used.
- **RNG / determinism (the crux).** `perm_p` is a label-permutation MC p-value seeded by
  `cfg$random_seed` ([score.R:13-16](../../R/score.R#L13-L16), [anchor_map.R:59](../../anchor_map.R#L59)).
  - **Each z-task calls `set.seed(cfg$random_seed)` at its top**, then runs the clusters×levels loop in the
    *same order* as today. ⇒ (a) the z = `h2_z_threshold` task reproduces the current serial run's RNG
    stream **byte-for-byte** (so primary == sensitivity[z==primary] == committed Phase-1/2 output); (b) every
    z-task is self-contained and deterministic ⇒ **invariant to thread count and to backend**. `future.seed
    = TRUE` is still passed (silences the parallel-RNG warning) but the explicit `set.seed` dominates.
  - **Do NOT parallelize `perm_p`.** Parallelizing the permutation draws would reorder RNG consumption and
    break the z==primary parity. Phase 3 parallelizes the **outer z axis only**; `score.R` is untouched.
    (ADD §7.4 also lists "threaded perm_p" — deliberately deferred; noted in §NOTES.)
- **Gate-per-z = full re-run.** `select_corr_source` depends on the gated trait set (coverage; the
  `cluster_profile` proxy is rebuilt from `g`), so it is re-run **inside** each z-task on that z's `g`.
  *Optimization:* for the TSV route, build the trait×trait matrix **once** over the union of all-z gated
  `trait_id`s and pass it as `trait_rg_override` to every z-task (avoids re-reading the big LDSC summary
  per z; `trait_rg_coverage`/`reindex_corr` still subset per-z, so the result is identical). For the `.rds`
  route the override is already the full panel block (z-independent). The `cluster_profile` proxy stays
  per-z (cheap, in-memory).
- **Method assumptions / failure modes.** Higher z ⇒ stricter gate ⇒ smaller N ⇒ fewer categories clear
  `min_category_n` (some clusters may drop to `ambiguous`/`weak` at high z — that is a *finding*, not a
  bug). `n_eff ≤ n` and `vif ≥ 1` must hold at every z. If a z gates a cluster below `min_category_n*2`,
  `score_cluster_level` returns `list()` ([score.R:23](../../R/score.R#L23)) → that cluster simply has no
  rows at that z → it must still appear in `sensitivity_z_labels` as `ambiguous/weak` (mirror
  `rank_and_label`'s empty-sub branch, [label.R:61-67](../../R/label.R#L61-L67)).

### Files to read / create
- READ: [anchor_map.R](../../anchor_map.R) — the driver to refactor (extract the per-z core; keep the
  log/`FINISHED`/write contract). Lines 83-104 are the loop to lift.
- READ: [R/gate.R](../../R/gate.R) — add an optional `z` arg to `apply_universe_gate`.
- READ: [R/redundancy.R](../../R/redundancy.R) `select_corr_source` (L119) — re-used per z unchanged.
- READ: [R/score.R](../../R/score.R) / [R/label.R](../../R/label.R) — **unchanged** (touching them risks
  the parity/RNG invariant).
- READ: [tests/run_tests.R](../../tests/run_tests.R), [tests/test_phase2.R](../../tests/test_phase2.R) —
  mirror the plain-`stopifnot` harness; reuse the synthetic `.rds` fixtures + the Carey configs.
- CREATE: `R/sensitivity.R` — `parallel_lapply()`, `score_at_z()` (per-z core), `run_sensitivity()`
  (sweep + stack + `label_stable`). Keeps the driver slim.
- CREATE: `tests/test_phase3.R` — primary-slice parity, thread-invariance, `label_stable`, gate-monotonicity.
- UPDATE: [R/io.R](../../R/io.R) `default_config()` — add `z_vector` (default `c(3,4,5)`).
- UPDATE: [anchor_map.R](../../anchor_map.R) — `--z-vector`, source `R/sensitivity.R`, call the sweep,
  write the two sensitivity TSVs, new log lines, `.SENS_*` col contracts.
- (no config edits required — the default `z_vector` applies to the existing Carey configs; optionally add
  an explicit `z_vector: [3, 4, 5]` comment to document it.)

---

## METHOD / IMPLEMENTATION PLAN

### Phase 1: per-z core extraction (`R/sensitivity.R` + `R/gate.R`)
1. `apply_universe_gate(df, cfg, z = NULL)` — default `z <- if (is.null(z)) as.numeric(cfg$h2_z_threshold)
   else as.numeric(z)`; use `z` in the `h2_z > z` filter ([gate.R:18](../../R/gate.R#L18)). **Zero
   behaviour change when `z` omitted** ⇒ Phase-1/2 callers unaffected.
2. `score_at_z(df, ont, cfg, sroot, z, trait_rg_override, emit)` — the lifted core, deterministic:
   ```
   set.seed(as.integer(cfg$random_seed))      # <- per-z reseed: parity + thread-invariance
   g   <- apply_universe_gate(df, cfg, z)
   g   <- attach_ontology(g, ont, cfg$ontology_key, cfg$levels)
   sel <- select_corr_source(g, cfg, sroot, trait_rg_override, emit_quiet)
   rows <- (for cl in clusters) (for lvl in levels) score_cluster_level(gc, lvl, sel$corr, cfg)
   if (!length(rows)) return(list(ranked = <empty>, labels = <all-ambiguous>, n_gated = nrow(g), ...))
   rank_and_label(do.call(rbind, rows), cfg)  ->  ranked, labels
   return list(ranked, labels, z = z, n_gated = nrow(g), n_clusters, source = sel$source, coverage)
   ```
   Build the all-ambiguous labels for clusters with no rows by reusing `rank_and_label`'s empty branch (it
   already emits `ambiguous/weak` for clusters absent from `prim`, [label.R:59-67](../../R/label.R#L59-L67))
   — i.e. feed it the full cluster set so missing clusters fall through to that branch. **Verify**: a
   cluster gated to zero categories still yields one `ambiguous/weak` label row.
3. `parallel_lapply(X, FUN, threads)` — **always via `future.apply`** (declared dependency): `workers <-
   max(1, min(threads, length(X)))`; `future::plan(if (workers == 1) "sequential" else "multicore", workers
   = workers)`; `on.exit(future::plan("sequential"), add = TRUE)`; `future.apply::future_lapply(X, FUN,
   future.seed = TRUE)`. (`future.seed = TRUE` only silences the parallel-RNG warning — the per-task
   `set.seed(random_seed)` in `score_at_z` dominates, so results are identical for any `workers`.) Set
   `setDTthreads(1)` inside the region (avoid data.table × workers oversubscription). The `workers == 1`
   → `"sequential"` branch is plan selection, **not** an availability fallback.

### Phase 2: the sweep (`run_sensitivity`)
4. `run_sensitivity(df, ont, cfg, sroot, z_vector, threads, trait_rg_override, emit)`:
   - `z_primary <- as.numeric(cfg$h2_z_threshold)`; `zs <- sort(unique(c(as.numeric(z_vector),
     z_primary)))` (**primary z is always in the sweep** ⇒ the parity slice always exists).
   - `res <- parallel_lapply(zs, function(z) score_at_z(df, ont, cfg, sroot, z, trait_rg_override,
     emit_quiet), threads)`.
   - Stack: `scores_stacked <- rbind over z of (res$ranked with z_threshold = z)`;
     `labels_stacked <- rbind over z of (res$labels with z_threshold = z)`.
   - `label_stable`: per `cluster_label`, `TRUE` iff `length(unique(auto_label)) == 1` across `zs`; join
     back onto `labels_stacked`.
   - Order: scores by `(z_threshold, level, cluster_label, rank)`; labels by `(z_threshold, cluster_label)`.
   - Return both stacked frames **and** the `res` element whose `z == z_primary` (for the driver to write
     the primary TSVs — guaranteeing primary == sweep[z==primary] from the *same* computation).

### Phase 3: driver integration (`anchor_map.R`)
5. Source `R/sensitivity.R`; add `--z-vector "3,4,5"` parsing (comma/space split → numeric) overriding
   `cfg$z_vector`. CLAUDE.md reserves `--z-vector` for exactly this.
6. Replace the inline single-z block ([anchor_map.R:83-104](../../anchor_map.R#L83-L104)) with: read
   ontology once; (TSV route) build `trait_rg_override` once over the union of all-z gated `trait_id`s when
   `vif_correlation ∈ {trait_rg, auto}` (else `NULL`); call `run_sensitivity(...)`.
7. Write the **primary** TSVs from the returned `z==primary` slice exactly as today (same `.SCORE_COLS`
   ordering/rounding, `eligible → True/False`, same paths) ⇒ byte-identical. Then write
   `sensitivity_z_scores.tsv` / `sensitivity_z_labels.tsv` (the `.SENS_*` contracts).
8. Log: keep `[start]/[config]/[gate]/[vif]/[write]/FINISHED`; **add** `[sweep] z ∈ {…} on N workers`, one
   `[sweep z=k] gated=… clusters=… source=… (coverage …%)` per z, `[stable] X/Y clusters label-stable`,
   and the two extra `[write]` lines. The `FINISHED` manifest lists all **four** TSVs.

### Phase 4: validation — see VALIDATION STRATEGY.

---

## STEP-BY-STEP TASKS (execute top to bottom; each atomic + checkable)

### UPDATE `R/gate.R` — overridable gate threshold
- **IMPLEMENT**: add `z = NULL` arg to `apply_universe_gate`; resolve to `cfg$h2_z_threshold` when `NULL`;
  use it in the `h2_z >` filter.
- **PATTERN**: minimal diff at [gate.R:18](../../R/gate.R#L18).
- **GOTCHA**: keep the default-path behaviour identical (Phase-1/2 callers pass no `z`).
- **VALIDATE**: `Rscript -e 'source("R/io.R");source("R/gate.R"); d<-read_long("/abs/cluster_trait_rg_long_with_p.tsv"); c4<-nrow(apply_universe_gate(d,modifyList(default_config(),list()) )); c5<-nrow(apply_universe_gate(d,default_config(),5)); stopifnot(c5<=c4)'`
  (higher z ⇒ ≤ rows).

### CREATE `R/sensitivity.R` — `parallel_lapply` + `score_at_z` + `run_sensitivity`
- **IMPLEMENT**: tasks 2-4 above. `score_at_z` re-seeds with `cfg$random_seed`; loops in the **same
  cluster/level order** as [anchor_map.R:97-100](../../anchor_map.R#L97-L100).
- **PATTERN**: lift [anchor_map.R:83-104](../../anchor_map.R#L83-L104) verbatim into `score_at_z`.
- **DATA/SCHEMA**: `run_sensitivity` returns stacked frames with `z_threshold` appended + the primary slice.
- **GOTCHA**: clusters with no scored category at a given z **must** still appear in that z's labels block as
  `ambiguous/weak` (drive `rank_and_label` with the full cluster universe). `setDTthreads(1)` inside workers.
- **VALIDATE**: covered by `tests/test_phase3.R` (next).

### CREATE `tests/test_phase3.R` — primary-slice parity (the key gate)
- **IMPLEMENT**: run `run_sensitivity` on the **anthro** Carey config (or synthetic `.rds`) with
  `z_vector=c(3,4,5)`, `threads=1`; assert the `z==4` slice of `scores_stacked` (drop `z_threshold`) equals
  the standalone Phase-1/2 primary output **column-for-column including `perm_p`** (`all.equal`, tol 0).
- **PATTERN**: plain `stopifnot`/`isTRUE(all.equal(...))` ([tests/test_phase2.R](../../tests/test_phase2.R)).
- **GOTCHA**: this only holds because `score_at_z` re-seeds with `random_seed` and preserves loop order —
  if it fails, the RNG-consumption order drifted.
- **VALIDATE**: `Rscript tests/test_phase3.R` → prints `primary-slice parity OK`.

### CREATE `tests/test_phase3.R` — thread-invariance
- **IMPLEMENT**: run `run_sensitivity` with `threads=1` and `threads=4`; assert `scores_stacked` and
  `labels_stacked` are **identical** (`all.equal`, tol 0), incl. `perm_p`.
- **GOTCHA**: invariance comes from the per-task `set.seed(random_seed)` (not from `future.seed`), so it
  holds for any `workers` under the `multicore` plan. `future.apply` is a hard dependency — the test
  `library(future.apply)`s directly (no backend branching).
- **VALIDATE**: `Rscript tests/test_phase3.R` → `thread-invariance OK`.

### CREATE `tests/test_phase3.R` — `label_stable` + gate monotonicity + sanity
- **IMPLEMENT**: (a) on the anthro config, assert C5_sub0 `auto_label=="Anthropometric"` and
  `anchor_shape=="sharp"` at **every** z∈{3,4,5} ⇒ its `label_stable=="True"`; (b) construct a frame where a
  cluster's label differs across z ⇒ `label_stable=="False"`; (c) `n_gated` strictly non-increasing in z;
  (d) every score row satisfies `n_eff ≤ n` and `vif ≥ 1` at every z.
- **GOTCHA**: `label_stable` is per-cluster broadcast — assert it's constant within a cluster across its z-rows.
- **VALIDATE**: `Rscript tests/test_phase3.R` → `label_stable OK; monotonic gate OK; sanity OK`.

### UPDATE `R/io.R` `default_config()` — `z_vector`
- **IMPLEMENT**: add `z_vector = c(3, 4, 5)` to the defaults list ([io.R:12-27](../../R/io.R#L12-L27));
  coerce `cfg$z_vector` to numeric in `load_config` (like `levels`).
- **GOTCHA**: additive only; **do not** touch any Phase-1/2 default. Existing configs gain the sweep via the
  default without edits.
- **VALIDATE**: `Rscript -e 'source("R/io.R"); stopifnot(identical(as.numeric(default_config()$z_vector), c(3,4,5)))'`

### UPDATE `anchor_map.R` — `--z-vector`, sweep call, two sensitivity writes, log
- **IMPLEMENT**: tasks 5-8. Add `.SENS_SCORE_COLS <- c(.SCORE_COLS, "z_threshold")` and
  `.SENS_LABEL_COLS <- c(.LABEL_COLS, "z_threshold", "label_stable")`; print `label_stable` as `True/False`.
- **GOTCHA**: build `trait_rg_override` once (TSV route) over the **union** of all-z gated `trait_id`s;
  keep the primary write path/columns/rounding identical; `FINISHED` manifest lists 4 TSVs.
- **VALIDATE**: full validation commands below (primary unchanged + sweep emitted).

---

## VALIDATION STRATEGY
No Python oracle for the sweep (the reference is single-z). Validate by:
- **Primary-slice parity (internal oracle):** the `z == h2_z_threshold` slice of `sensitivity_z_scores.tsv`
  equals the committed `category_anchor_scores.tsv` **byte-for-byte, incl. `perm_p`** — re-using Phase 1/2's
  already-oracle-validated numbers. This is the headline gate.
- **Thread-invariance:** `--threads 1` vs `--threads 8` produce identical sensitivity TSVs (diff = empty).
- **Positive control (stability):** anthro C5_sub0 → `Anthropometric [sharp]` at z∈{3,4,5} ⇒
  `label_stable=True`; the z=4 row reproduces ADD §8 values (`auc_abs=0.9164`, `pooled_rg=0.2473`,
  `vif_p≈0.03489`, `q≈0.005497`, `rank=1`).
- **Negative / forbidden-FP control:** on the **disease** track, no `anchor_eligible=FALSE` category
  (Quantitative/Lab) becomes any cluster's `auto_label` at **any** swept z; a weak/ambiguous cluster
  (e.g. anthro C3, ADD §8) stays `ambiguous/weak` across z (its own `label_stable` should be `True`).
- **Gate monotonicity:** per-z gated row count is non-increasing in z (logged + asserted).
- **Schema/sanity:** `sensitivity_z_scores.tsv` has `c(.SCORE_COLS,"z_threshold")`; `sensitivity_z_labels`
  has `c(.LABEL_COLS,"z_threshold","label_stable")`; every cluster present at every z; `n_eff ≤ n`,
  `vif ≥ 1`; no NaN explosion; the two **primary** TSVs unchanged (oracle comparison still clean).

## VALIDATION COMMANDS (run all; zero parity/control failures)
```bash
# 0. deps — future + future.apply are REQUIRED (installed locally: future 1.x, future.apply 1.20.2)
Rscript -e 'stopifnot(all(c("future","future.apply") %in% rownames(installed.packages()))); cat("future.apply present\n")'

# 1. unit + Phase-3 battery (parity slice, thread-invariance, label_stable, monotonic gate, sanity)
Rscript tests/run_tests.R
Rscript tests/test_phase2.R
Rscript tests/test_phase3.R

# 2. PRIMARY PARITY — the two primary TSVs must be byte-identical to the committed run
cp results/carey_rint15_anthro/category_anchor_scores.tsv /tmp/anthro_scores.before.tsv
Rscript anchor_map.R --config configs/carey_rint15_anthro.yaml --threads 4
diff -q /tmp/anthro_scores.before.tsv results/carey_rint15_anthro/category_anchor_scores.tsv   # -> identical
Rscript validation/compare_oracle.R \
  --r-out  results/carey_rint15_anthro/category_anchor_scores.tsv \
  --oracle "../UKBB_CLUSTER_GWAS/scripts/cluster_anchoring/output/carey_rint15_anthro/category_anchor_scores.tsv"

# 3. SWEEP emitted + primary-slice parity + z values present
Rscript -e 'library(data.table); s<-fread("results/carey_rint15_anthro/sensitivity_z_scores.tsv"); \
  stopifnot("z_threshold" %in% names(s), all(c(3,4,5) %in% unique(s$z_threshold))); \
  p<-fread("results/carey_rint15_anthro/category_anchor_scores.tsv"); \
  slice<-s[z_threshold==4][, !"z_threshold"]; setcolorder(slice, names(p)); \
  stopifnot(isTRUE(all.equal(as.data.frame(slice), as.data.frame(p), check.attributes=FALSE))); \
  cat("sweep+parity OK\n")'

# 4. label_stable + forbidden-FP across z (disease track)
Rscript anchor_map.R --config configs/carey_rint15.yaml --threads 4
Rscript -e 'library(data.table); l<-fread("results/carey_rint15_anthro/sensitivity_z_labels.tsv"); \
  c5<-l[cluster_label=="C5_sub0"]; stopifnot(all(c5$auto_label=="Anthropometric"), all(c5$label_stable=="True")); \
  d<-fread("results/carey_rint15/sensitivity_z_labels.tsv"); \
  stopifnot(!any(d$auto_label %in% c("Quantitative","Lab"))); cat("label_stable + forbidden-FP OK\n")'

# 5. THREAD-INVARIANCE
Rscript anchor_map.R --config configs/carey_rint15_anthro.yaml --threads 1
cp results/carey_rint15_anthro/sensitivity_z_scores.tsv /tmp/sens.t1.tsv
Rscript anchor_map.R --config configs/carey_rint15_anthro.yaml --threads 8
diff -q /tmp/sens.t1.tsv results/carey_rint15_anthro/sensitivity_z_scores.tsv   # -> identical
```

## ACCEPTANCE CRITERIA
- [ ] **Primary-slice parity:** `sensitivity_z_scores.tsv[z==h2_z_threshold]` (minus `z_threshold`) ==
      committed `category_anchor_scores.tsv` byte-for-byte (incl. `perm_p`); the two primary TSVs unchanged;
      oracle comparison still clean.
- [ ] **Thread-invariant:** identical sensitivity TSVs for `--threads 1` vs `--threads 8` (same
      `future.apply` backend, varying only `workers`).
- [ ] **Sweep emitted:** both `sensitivity_z_*.tsv` carry every cluster at every swept z; `z_vector` default
      `{3,4,5}` plus the primary z; `--z-vector` overrides it.
- [ ] **`label_stable`** correct: `True` iff `auto_label` constant across all swept z (per cluster); printed
      `True/False`. Anthro C5_sub0 stable; forbidden-FP categories never label at any z.
- [ ] Gate counts non-increasing in z; `n_eff ≤ n`, `vif ≥ 1` at every z; no NaN explosion.
- [ ] Provenance: log records `z_vector`, worker count/backend, per-z gate counts + redundancy source, the
      `X/Y label-stable` summary, all four output files, and a `FINISHED` line.

## NOTES
- **Parity is engineered, not hoped-for.** The primary TSVs are *defined* as the sweep's primary-z slice,
  and each z-task re-seeds with `random_seed` while preserving the cluster/level loop order ⇒ the primary
  slice is bit-identical to Phase-1/2 by construction, and the whole sweep is invariant to thread count.
- **Why `perm_p` is left serial.** Parallelizing the permutation draws would reorder RNG consumption and
  break the z==primary parity. ADD §7.4's "threaded perm_p" is therefore **deferred**; if ever added it must
  use per-(cluster,level,n_in) deterministic sub-seeds that don't perturb the serial stream — out of scope
  here, noted as future work.
- **Backend is `future.apply`, full stop.** Chosen on merit (backend-agnostic `plan()`, parallel-safe RNG,
  clean `lapply` drop-in), not on local availability — `future`/`future.apply` are now installed locally
  (and will be pinned in the Phase-5 P3M snapshot). The code `library(future.apply)`s directly; no
  availability branching, no `parallel::mclapply`. `workers` varies (1 → `sequential` plan, else
  `multicore`) but results are worker-count-invariant because each z-task re-seeds with `random_seed`.
- **Subdir layout deferred.** ADD §7.5's `primary/ sensitivity/ figures/ logs/` reorg lands in Phase 5 with
  the Nextflow `output {}`/`outputDir` publish; Phase 3 keeps the flat `out_dir` to preserve Phase-1/2
  parity and the oracle path.
- **`z_vector` default `{3,4,5}`** matches the ADD §8 stability window and the cheapest defensible bracket
  (3 re-runs); the broader `{2,3,4,5,6,7}` (ADD §7.1) is available via config/`--z-vector` when wanted.
```
