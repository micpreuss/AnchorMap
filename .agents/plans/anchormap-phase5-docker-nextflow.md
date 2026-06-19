# Analysis: AnchorMap Phase 5 — pinned Docker image (the tool) + Nextflow container-validation harness

Validate the schema contracts and method usage against the actual code/data before implementing.
Pay special attention to: the exact runtime R dependency set (NOT what the ADD lists), how the engine
resolves config paths (stage-root-relative), and the three parent container fixes + *why* each exists.

> **Scope decision (locked with the user, 2026-06-18):** **AnchorMap proper = the R engine + the pinned
> Docker image.** That image is the shippable, reproducible tool and the primary run interface
> (`docker run anchormap:0.1.0 Rscript anchor_map.R --config <yaml>`). **Nextflow is NOT an orchestration
> layer here** — there is no DAG to schedule and no scale that needs a batch executor (the engine is one
> fast single-process job; its z-sweep/`perm_p` parallelism is internal, via `future`/`setDTthreads`).
> Nextflow's *only* job in Phase 5 is a **thin harness that proves the container runs flawlessly under
> Nextflow** — i.e. it exercises the three Nextflow-specific container failure modes that a plain
> `docker run` cannot catch. This is a deliberate, recorded divergence from ADD §5/§7.1/§12 (which framed
> Nextflow as a production `ANCHORMAP` process); see NOTES + the offer to update the ADD.
>
> The engine is byte-for-byte validated (Phases 1–4). Phase 5 adds **zero** lines to `R/`.

## Question & object
Two distinct artifacts, two distinct questions:
1. **The tool (deliverable):** *Can the validated engine run reproducibly — pinned deps, no host R, no
   host paths — as a container?* Object = the image `anchormap:0.1.0` carrying the engine + figures +
   a build-time self-test. Inference it licenses: any AnchorMap run is reproducible from a version tag + a config.
2. **The harness (a test, not a deliverable to *use*):** *Does that image execute flawlessly under
   Nextflow's execution model?* Object = a one-process `ANCHORMAP_SMOKE` workflow. It validates the
   container ⇄ Nextflow contract, nothing more.

## Analyst story
As an analyst / I want AnchorMap shipped as one pinned Docker image I run directly with a config, plus a
minimal Nextflow smoke test proving that image obeys the lab's container contract (no exit-126, GCS-FUSE
write perms, working `ps`/trace) / So that AnchorMap is reproducible and *Batch-ready* without paying for
a production orchestration layer the workload doesn't need.

## Pipeline position
- **The tool's run interface (primary):** `docker run … anchormap:0.1.0 Rscript /opt/anchormap/anchor_map.R --config <yaml>`
  and the figures CLI `Rscript /opt/anchormap/R/plot_anchors.R --config <plots.yaml>`. Inputs are the same
  as Phases 1–4 (Input A/B/D TSVs or Input C `.rds`) from the sibling `UKBB_CLUSTER_GWAS`.
- **The harness:** runs the image's **self-contained synthetic fixture** (no external data) and asserts the
  engine's outputs land where Nextflow captures them. Local `test` profile = the CI gate; `gcp` profile =
  a one-time Google Batch check (the only place the FUSE/USER-root fix is genuinely exercised).
- **Patterns to mirror:**
  - Container fixes + rationale: [../UKBB_CLUSTER_GWAS/docker/postgwas/Dockerfile](../UKBB_CLUSTER_GWAS/docker/postgwas/Dockerfile)
    + [README.md](../UKBB_CLUSTER_GWAS/docker/postgwas/README.md) (exit-126 / GCS-FUSE perms).
  - Rocker base + apt + R installs: [../UKBB_CLUSTER_GWAS/docker/genomicsem/Dockerfile](../UKBB_CLUSTER_GWAS/docker/genomicsem/Dockerfile).
  - Process threading exports + `output {}` block: [../UKBB_CLUSTER_GWAS/scripts/pipeline/main.nf](../UKBB_CLUSTER_GWAS/scripts/pipeline/main.nf)
    (L143–157 the `*_NUM_THREADS=1` exports + `--threads ${task.cpus}`; L332–349 `output {}` with `path/mode 'copy'`).
  - Profiles + `google.batch.spot`: [../UKBB_CLUSTER_GWAS/scripts/pipeline/nextflow.config](../UKBB_CLUSTER_GWAS/scripts/pipeline/nextflow.config).

## CONTEXT REFERENCES — READ BEFORE IMPLEMENTING

### ⚠ Decision-grade findings (verified from the code; override the ADD where noted)

1. **GenomicSEM is NOT a runtime dependency.** [R/ingest_rds.R](R/ingest_rds.R) reads the `.rds` with base
   `readRDS()` (L39); `grep -rn GenomicSEM R/` returns only comments. The engine *consumes* the `ldsc()`
   artifact, it never *runs* `ldsc()`. → **Drop GenomicSEM from the image** (vs ADD §7.3): removes a heavy
   GitHub install + dep tree + a pinning surface. Note the divergence in `docker/README.md`.

2. **Exact runtime R deps** (`grep -rEn "library\(|::|requireNamespace" R/ anchor_map.R`):
   - Engine: `data.table`, `yaml`, `poolr` (+ `Matrix`, poolr dep), `future`, `future.apply`; base `stats/utils/grDevices`.
   - Figures: `ggplot2`, `patchwork`, `scales`, `ggrepel`, and **optional** `ragg` ([R/plot.R](R/plot.R) L337:
     `requireNamespace("ragg")` else base `png(type="cairo")`; PDF via `grDevices::cairo_pdf` L344).
   - **No `argparse`/`optparse`/`jsonlite`** — both CLIs hand-roll `parse_args` (anchor_map.R L28). Don't install them.

3. **Config paths resolve stage-root-relative** ([R/io.R](R/io.R) L44–51): `stage_root_of()` = config dir
   (or its parent if under `configs/`); `resolve_path()` keeps absolute paths verbatim, else joins to the
   stage root; `out_dir` resolves the same way. → The harness sets the synthetic config's input paths to the
   **absolute baked-in fixture paths** and `out_dir` to an **absolute work-dir path**, so the engine writes
   straight into Nextflow's work dir (this is what makes the `gcp` profile a real FUSE-write test).

4. **The self-contained smoke already works.** `Rscript anchor_map.R --config configs/synthetic_rds.yaml`
   runs the whole engine from `tests/fixtures/synthetic_ldsc_panel.rds` + `synthetic_panel_ontology.tsv`
   (**no external paths**), recovering `C5_sub0 → anthro [sharp]` + all 4 TSVs + `FINISHED`
   (confirmed in `results/synthetic_rds/anchormap.log`). It is both the **build-time** self-test and the
   **Nextflow** smoke payload. The C5_sub0 *real-data* oracle (ADD §8) stays a host check (needs FinnGen inputs).

### Where each container fix is actually testable (drives the two profiles)
| Fix | Failure mode it prevents | Caught by plain `docker run`? | Caught by **local** NF? | Needs **Batch**? |
|---|---|---|---|---|
| `ENTRYPOINT []` | upstream entrypoint intercepts `.command.run` → exit 126 | ❌ | ✅ | — |
| `procps`/`ps` | Nextflow trace/metrics poll `ps` | ❌ | ✅ (`-with-trace`) | — |
| `USER root` | can't write GCS-FUSE-mounted work dir | ❌ | ❌ (local bind mount as root won't reproduce) | ✅ **only here** |

→ Local `test` profile validates entrypoint + procps + output-capture; the FUSE/USER-root fix is only
genuinely validated by the one-time `gcp` Batch run. (rocker has no problematic entrypoint, so `ENTRYPOINT []`
is belt-and-suspenders — the test confirms it stays clean.)

### Input/Output schemas (unchanged from Phases 1–4 — the engine asserts them)
- Inputs A/B/C/D and the 4 output TSVs + figures exactly as in Phases 1–4 (anchor_map.R `.SCORE_COLS`
  L44–47, `.LABEL_COLS` L48–49; plot_anchors.R outputs L7–11). Phase 5 does **not** touch these. The harness
  only checks the synthetic-route outputs *exist + parse + end in FINISHED* — it does not re-validate values
  (that's Phases 1–3's job).

### Files to read / create
- READ: postgwas/Dockerfile + README, genomicsem/Dockerfile, pipeline/main.nf (threading + `output {}`),
  pipeline/nextflow.config (profiles + `google`), [anchor_map.R](anchor_map.R), [R/io.R](R/io.R),
  [R/plot_anchors.R](R/plot_anchors.R), [configs/synthetic_rds.yaml](configs/synthetic_rds.yaml).
- CREATE: `docker/Dockerfile` — pinned rocker image (engine + figures baked in; build-time engine **and** figure smoke; procps, USER root, ENTRYPOINT []).
- CREATE: `docker/README.md` — build/push + the three-fixes rationale + the GenomicSEM-omission note (mirror postgwas/README.md).
- CREATE: `.dockerignore` (repo root) — exclude `claude-science-scaffold/`, `.git/`, `.agents/`, `results/`, `validation/`, R history.
- CREATE: `configs/synthetic_rds_plots.yaml` — a tiny single-track plots config over the synthetic run's TSVs (for the build-time figure smoke).
- CREATE: `nextflow/main.nf` — one `ANCHORMAP_SMOKE` process + a trivial workflow + a minimal `output {}` (capture proof only).
- CREATE: `nextflow/nextflow.config` — params + `test` (local) + `gcp` (Batch) profiles + `google` block + resource defaults.
- CREATE: `nextflow/params/test.yaml` (local image tag) and `nextflow/params/gcp.yaml` (Artifact Registry tag + bucket/slug/user).

## METHOD / IMPLEMENTATION PLAN

### Phase 1 — The pinned image (the deliverable): `docker/Dockerfile`
The image **is** the tool: engine + figures baked in, self-validating at build time.
```dockerfile
FROM rocker/r-ver:4.4.2

# system deps: procps (NF ps/trace) + build libs + figure-render libs (cairo/freetype/png/jpeg/tiff)
RUN apt-get update && apt-get install -y --no-install-recommends \
      procps libcurl4-openssl-dev libssl-dev libxml2-dev \
      libcairo2-dev libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev \
      cmake git ca-certificates locales \
    && rm -rf /var/lib/apt/lists/*

# pin CRAN to a dated P3M snapshot (confirm <codename> via /etc/os-release; pick a date with ggplot2 >= 4.0)
RUN echo 'options(repos=c(P3M="https://packagemanager.posit.co/cran/__linux__/<codename>/<YYYY-MM-DD>"))' \
      >> "${R_HOME}/etc/Rprofile.site" \
 && R -e "install.packages(c('data.table','yaml','poolr','Matrix','future','future.apply', \
                             'ggplot2','patchwork','scales','ggrepel','ragg'))" \
 && R -e "if(!all(sapply(c('data.table','yaml','poolr','future','future.apply','ggplot2', \
                           'patchwork','scales','ggrepel','ragg'),requireNamespace))) quit(status=1)"

# bake the tool in (pinned with the image tag) — NO GenomicSEM, NO argparse/optparse
WORKDIR /opt/anchormap
COPY R/ ./R/
COPY anchor_map.R ./anchor_map.R
COPY configs/ ./configs/
COPY ontology/ ./ontology/
COPY tests/fixtures/synthetic_ldsc_panel.rds tests/fixtures/synthetic_panel_ontology.tsv ./tests/fixtures/

# build-time self-test #1 (engine): synthetic .rds run must recover C5_sub0 -> anthro [sharp] + FINISHED
RUN Rscript anchor_map.R --config configs/synthetic_rds.yaml --threads 2 \
 && grep -Eq 'C5_sub0.*anthro.*sharp' results/synthetic_rds/cluster_anchor_labels.tsv \
 && grep -q '^FINISHED ok' results/synthetic_rds/anchormap.log

# build-time self-test #2 (figures): the ggplot/ragg/cairo stack must render a PNG from those TSVs
RUN Rscript R/plot_anchors.R --config configs/synthetic_rds_plots.yaml \
 && ls results/synthetic_rds/figures/anchor_lollipops_*.png \
 && rm -rf results/synthetic_rds

# project container fixes (postgwas/README.md): write perms on GCS-FUSE + don't intercept .command.run
USER root
ENTRYPOINT []
```
Notes: **`ragg` is installed** (not left optional) so PNG rendering is headless-robust without X11; PDF uses
`cairo_pdf` ⇒ `libcairo2-dev`. Verify the TSV separator (`head -3 cluster_anchor_labels.tsv`) before fixing
the `grep`. `--platform linux/amd64` for Google Batch when building on Apple Silicon.

### Phase 2 — The Nextflow container-validation harness: `nextflow/main.nf`
One process. No fan-out, no figures pipeline, no production output tree, no real-data params.
```groovy
nextflow.enable.dsl = 2
params.anchormap_container = null   // set per-profile (local tag vs Artifact Registry tag)
params.smoke_cpus = 2

process ANCHORMAP_SMOKE {
    container params.anchormap_container
    cpus   params.smoke_cpus
    memory '4 GB'

    output:
    path "smoke_out/*.tsv",       emit: tsv
    path "smoke_out/anchormap.log", emit: log

    script:
    """
    export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
           VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1

    # point the baked synthetic config at absolute in-image fixtures + an absolute work-dir out_dir,
    # so the ENGINE writes straight into Nextflow's work dir (this is the real FUSE-write test on Batch)
    Rscript -e '
      cfg <- yaml::read_yaml("/opt/anchormap/configs/synthetic_rds.yaml")
      cfg\$rds      <- "/opt/anchormap/tests/fixtures/synthetic_ldsc_panel.rds"
      cfg\$ontology <- "/opt/anchormap/tests/fixtures/synthetic_panel_ontology.tsv"
      cfg\$out_dir  <- file.path(getwd(), "smoke_out")
      yaml::write_yaml(cfg, "run_config.yaml")'

    Rscript /opt/anchormap/anchor_map.R --config run_config.yaml --threads ${task.cpus}

    grep -q '^FINISHED ok' smoke_out/anchormap.log   # fail the task if the engine didn't finish
    """
}

workflow {
    main: ANCHORMAP_SMOKE()
    publish:
      smoke_tsv = ANCHORMAP_SMOKE.out.tsv
      smoke_log = ANCHORMAP_SMOKE.out.log
}

output {
    smoke_tsv { path 'smoke'; mode 'copy' }
    smoke_log { path 'smoke'; mode 'copy' }
}
```
Why this proves the contract: the task runs *through* Nextflow's `.command.run` wrapper (entrypoint OK),
`ps`-based trace works (procps), the engine writes into the work dir and Nextflow captures it (output
capture; on Batch that work dir is the GCS-FUSE mount → USER-root/FUSE proven). `output {}` + `outputDir`,
never `publishDir` (parent convention).

### Phase 3 — Profiles + params: `nextflow/nextflow.config`
- `params`: `anchormap_container`, `smoke_cpus`, plus `gcp_bucket`/`project_slug`/`user` for the gcp work/outputDir.
- `process { executor='local'; cpus=1; memory='4 GB'; time='30m' }` defaults.
- `profiles`:
  - **`test`** (the CI gate): `process.executor='local'; docker.enabled=true; workDir='work'; outputDir='results/anchormap_smoke'`.
  - **`gcp`** (one-time Batch/FUSE check): `process { executor='google-batch'; errorStrategy='retry'; maxRetries=3; cache='lenient' }; docker.enabled=true; workDir="${params.gcp_bucket}/AnchorMap/work/${params.user}"; outputDir="${params.gcp_bucket}/AnchorMap/results/${params.user}"`.
- `google { project='lencz-lab-cogent-1'; location='us-central1'; batch { spot=true; network=…; subnetwork=… } }` — copy verbatim from the parent nextflow.config.
- `params/test.yaml`: `anchormap_container: anchormap:0.1.0` (local). `params/gcp.yaml`: the Artifact Registry
  tag `us-central1-docker.pkg.dev/lencz-lab-cogent-1/docker-images/anchormap:0.1.0` + bucket/slug/user.

### Phase 4 — Docs + wiring
- `docker/README.md` (mirror postgwas/README.md): Overview; pinning rationale; the three fixes + *why*; the
  GenomicSEM-omission note; build+push to Artifact Registry `anchormap:0.1.0` (`--platform linux/amd64`,
  by version tag never `latest`).
- Update [CLAUDE.md](CLAUDE.md) + [README.md](README.md): flip Phase-5 rows to "built"; add `docker/` +
  `nextflow/` to the tree; document **primary run = `docker run …`**, and **Nextflow = container-validation
  harness only** (test profile = CI gate; gcp profile = one-time Batch smoke). Use the `create-readme` skill for the README.

## STEP-BY-STEP TASKS (execute top to bottom; each atomic + checkable)

### CREATE `.dockerignore` (repo root)
- **IMPLEMENT**: exclude `claude-science-scaffold/`, `.git/`, `.agents/`, `results/`, `validation/`, `.DS_Store`, `*.Rhistory`; keep `R/`, `anchor_map.R`, `configs/`, `ontology/`, `tests/fixtures/synthetic_*`.
- **GOTCHA**: `claude-science-scaffold/` is a separate gitignored repo — never COPY it in.
- **VALIDATE**: `du -sh` of non-ignored files is small (≪ the scaffold).

### CREATE `configs/synthetic_rds_plots.yaml`
- **IMPLEMENT**: a single-track plots config over the synthetic run: `out_dir: results/synthetic_rds/figures`;
  one `tracks:` entry `{name: synthetic, level: anthro_class, scores: results/synthetic_rds/category_anchor_scores.tsv, labels: results/synthetic_rds/cluster_anchor_labels.tsv}`.
- **PATTERN**: [configs/carey_rint15_plots.yaml](configs/carey_rint15_plots.yaml).
- **GOTCHA**: synthetic run has 3 clusters → specificity panels skip (`spec_min_clusters` default 5); lollipop/dot-heatmap still render — that's enough to validate the figure stack. plot_anchors.R L94–97 skips gracefully; don't treat the skip as failure.
- **VALIDATE**: `Rscript R/plot_anchors.R --config configs/synthetic_rds_plots.yaml` writes `results/synthetic_rds/figures/anchor_lollipops_synthetic.png` (after a synthetic engine run).

### CREATE `docker/Dockerfile` + `docker/README.md`
- **IMPLEMENT**: the Phase-1 skeleton — rocker:4.4.2, sysdeps incl. procps+cairo/freetype, P3M-pinned installs
  (engine + figure deps + ragg; **no GenomicSEM/argparse/optparse**), bake engine+figures, **two** build-time
  self-tests (engine recovers C5_sub0 anthro sharp; figures render a PNG), `USER root`, `ENTRYPOINT []`.
  README mirrors postgwas/README.md + the GenomicSEM-omission note + the build/push block.
- **PATTERN**: genomicsem/Dockerfile (base+installs), postgwas/Dockerfile + README (fixes + rationale).
- **GOTCHA**: confirm `/etc/os-release` codename for the P3M URL; pick a snapshot with `ggplot2 ≥ 4.0` (host validated on 4.0.3 — see NOTES); `--platform linux/amd64`.
- **VALIDATE**: `docker build --platform linux/amd64 -t anchormap:0.1.0 -f docker/Dockerfile .` (the build *is* both self-tests).

### CREATE `nextflow/main.nf`, `nextflow/nextflow.config`, `nextflow/params/{test,gcp}.yaml`
- **IMPLEMENT**: the single `ANCHORMAP_SMOKE` process (threading exports + config-rewrite snippet + baked-in
  `Rscript` + `FINISHED` assert), trivial workflow, minimal `output {}`; `test`/`gcp` profiles + `google` block; the two params files.
- **PATTERN**: pipeline/main.nf L143–157 + L332–349; pipeline/nextflow.config.
- **DATA/SCHEMA**: process emits `smoke_out/*.tsv` + `smoke_out/anchormap.log`; no external inputs staged (payload is the baked-in synthetic fixture).
- **GOTCHA**: escape `$` as `\$` inside the Groovy double-quoted script for the inline R (`cfg\$rds`); `out_dir` must be **absolute** (`file.path(getwd(), …)`) so it lands in the work dir, not in-image.
- **VALIDATE**: `nextflow run nextflow/main.nf -profile test -params-file nextflow/params/test.yaml -preview`.

### RUN the local CI-gate smoke
- **IMPLEMENT**: `-profile test` end-to-end; assert published `results/anchormap_smoke/smoke/{category_anchor_scores.tsv,anchormap.log}` exist and the log ends `FINISHED ok`.
- **VALIDATE**: see Validation Commands #4–5.

### (one-time, manual) RUN the `gcp` Batch smoke — the only real FUSE/USER-root check
- **IMPLEMENT**: push the image to Artifact Registry; `nextflow run … -profile gcp -params-file nextflow/params/gcp.yaml`; confirm the task succeeds (no exit-126, no permission-denied on the GCS-FUSE work dir) and outputs publish to the gs:// outputDir.
- **GOTCHA**: requires GCP creds + a pushed image + network/subnetwork; this *costs* (spot VM) — it's a documented one-time gate, not a CI step.
- **VALIDATE**: Batch task state `SUCCEEDED`; `gsutil ls <outputDir>/smoke/` shows the TSVs + log.

### UPDATE [CLAUDE.md](CLAUDE.md) + [README.md](README.md)
- **IMPLEMENT**: flip Phase-5 rows to "built"; add `docker/`+`nextflow/` to the tree; document the primary
  run (`docker run`) and the Nextflow harness's *test-only* role (test=CI gate, gcp=one-time Batch smoke).
- **GOTCHA**: don't branch (CLAUDE.md convention); Conventional Commits, no `Co-Authored-By`.
- **VALIDATE**: links resolve; "Running stages" shows the container + the two Nextflow profiles.

## VALIDATION STRATEGY
- **Build = self-test (primary).** The Dockerfile's two build-time runs (engine recovers `C5_sub0 → anthro
  [sharp]` + `FINISHED`; figures render a PNG) fail the build on any engine/dep/figure-stack regression.
- **Determinism control.** In-container, run the synthetic engine `--threads 1` vs `--threads 4`;
  `cluster_anchor_labels.tsv` byte-identical (Phase-3 thread-invariance must survive containerisation).
- **Container-fixes check.** `docker run` asserts `ps` present, `whoami`=root, empty entrypoint honoured.
- **Local Nextflow contract (CI gate).** `-profile test` proves entrypoint + procps + output-capture; the
  engine's `FINISHED` log is captured under `outputDir`.
- **Batch contract (one-time).** `-profile gcp` is the only run that exercises USER-root/GCS-FUSE writes.
- **Positive control (host, real data).** Outside the container, C5_sub0 anthro §8 values remain ground
  truth; the containerised run on the same inputs reproduces them.

## VALIDATION COMMANDS (run all; zero failures)
```bash
# 1. Build (build-time self-tests = engine C5_sub0 anthro sharp + a rendered figure PNG)
docker build --platform linux/amd64 -t anchormap:0.1.0 -f docker/Dockerfile .

# 2. Container fixes
docker run --rm anchormap:0.1.0 which ps         # procps
docker run --rm anchormap:0.1.0 whoami           # root
docker run --rm --entrypoint= anchormap:0.1.0 bash -lc 'echo entrypoint-empty-ok'

# 3. In-container engine determinism (thread-invariance survives containerisation)
docker run --rm anchormap:0.1.0 bash -lc \
 'cd /opt/anchormap && Rscript anchor_map.R --config configs/synthetic_rds.yaml --threads 1 && cp results/synthetic_rds/cluster_anchor_labels.tsv /tmp/t1 && \
  Rscript anchor_map.R --config configs/synthetic_rds.yaml --threads 4 && diff -q /tmp/t1 results/synthetic_rds/cluster_anchor_labels.tsv && echo DETERMINISM_OK'

# 4. Nextflow DAG dry-run
nextflow run nextflow/main.nf -profile test -params-file nextflow/params/test.yaml -preview

# 5. Local CI-gate smoke → captured outputs + FINISHED
nextflow run nextflow/main.nf -profile test -params-file nextflow/params/test.yaml -with-trace
ls results/anchormap_smoke/smoke/category_anchor_scores.tsv results/anchormap_smoke/smoke/anchormap.log
grep -q '^FINISHED ok' results/anchormap_smoke/smoke/anchormap.log && echo SMOKE_OK

# 6. (one-time, manual) Batch/FUSE check — requires push + GCP creds
#   nextflow run nextflow/main.nf -profile gcp -params-file nextflow/params/gcp.yaml

# 7. (host, optional) real-data positive control parity
Rscript anchor_map.R --config configs/carey_rint15_anthro.yaml
grep -E 'C5_sub0.*Anthropometric.*sharp' results/carey_rint15_anthro/cluster_anchor_labels.tsv
```

## ACCEPTANCE CRITERIA
- [ ] Image builds `--platform linux/amd64`; **both** build-time self-tests pass (engine C5_sub0 anthro sharp + figure PNG).
- [ ] Container fixes verified: `ps` present, `whoami`=root, empty `ENTRYPOINT` honoured.
- [ ] CRAN deps pinned via dated P3M snapshot; **no GenomicSEM / argparse / optparse**; `ragg` present.
- [ ] In-container engine thread-invariant (Phase-3 determinism preserved).
- [ ] `nextflow … -preview` resolves; `-profile test` runs the one-process harness on the synthetic `.rds` (no host paths) and captures outputs + `FINISHED`.
- [ ] `gcp` profile authored + the one-time Batch smoke documented (the FUSE/USER-root validation path exists).
- [ ] Engine outputs unchanged from Phases 1–4 (the harness checks existence/parse/FINISHED, not values).
- [ ] `output {}` + `outputDir` used (not `publishDir`); image referenced by version tag, never `latest`.
- [ ] **No `R/` engine lines changed.**
- [ ] Docs flipped to "built": primary run = `docker run`; Nextflow = validation harness (test gate + one-time gcp). ADD divergence noted (NOTES) + ADD update offered.

## NOTES
- **Why Nextflow is a test, not orchestration (recorded rationale).** AnchorMap is a single fast
  single-process job with internal `future`/`setDTthreads` parallelism — there is no DAG to schedule and no
  fan-out/scale that a batch executor accelerates. A production Nextflow layer would be ceremony around one
  task. The value Nextflow *does* add is catching the three Nextflow-specific container failure modes that a
  plain `docker run` cannot (exit-126 / FUSE perms / `ps`). So the harness is scoped to exactly that.
- **Deliberate divergence from the ADD.** ADD §5 ("Compute: ✅ Nextflow DSL2 process"), §7.1 (engine+figures
  in the DSL2 flow), and §12 Phase 5 ("DSL2 `ANCHORMAP` process") imply Nextflow as a production run path.
  The user has scoped it down to a container-validation harness (figures + real runs happen via `docker run`,
  not Nextflow). **Offer to update ADD §5/§7.1/§12 Phase 5** to "container is the tool; Nextflow validates the
  container contract; full orchestration deferred to §14." The ADD already defers multi-cohort batch
  orchestration (§5) and nf-core packaging (§14), so this is consistent with its own future-work framing.
- **GenomicSEM omission** (finding #1): justified — engine reads the `.rds` with base `readRDS`. Reintroduce a
  pinned `remotes::install_github("GenomicSEM/GenomicSEM@<commit>")` only if a future stage runs `ldsc()`
  inside AnchorMap (explicitly out of scope, ADD §5).
- **ggplot2 4.0 pin risk (the single most likely silent break).** The figures module was authored against
  `ggplot2 4.0.3` (host). A P3M snapshot pinning `ggplot2 3.5.x` could alter/break rendering — hence the
  build-time **figure** self-test (self-test #2) renders a PNG so a bad snapshot fails the build, not a user.
- **Engine untouched.** The only Nextflow glue is the in-script synthetic-config rewrite (absolute fixture
  paths + absolute work-dir `out_dir`); it lives in the process script, not in `R/`. If you find yourself
  editing `R/*.R`, stop.
- **Scale/locality.** Workload is tiny (synthetic 0.2s; real run seconds–minutes). The `gcp` profile + spot
  VM is for *contract parity* with the parent (and the one-time FUSE check), not throughput; default few-GB
  memory/disk suffice.
```
