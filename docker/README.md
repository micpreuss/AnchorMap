# AnchorMap image

## Overview

**Purpose:** Ships the validated AnchorMap R engine (Phases 1–4: gate → redundancy → score → label →
sensitivity z-sweep → ggplot2 figures) as a single pinned, version-constrained container. **This image
*is* the tool** — the primary way to run AnchorMap:

The image installs the `anchormap` R package and puts `anchor_map` / `plot_anchors` shims on `PATH`:

```bash
# engine (rg long-TSV or GenomicSEM .rds route, per your config)
docker run --rm -v "$PWD:/work" -w /work \
    ghcr.io/micpreuss/anchormap:0.1.2 \
    anchor_map --config configs/your_run.yaml --out-dir results/your_run --threads 4

# figures (reads the scored TSVs the engine wrote)
docker run --rm -v "$PWD:/work" -w /work \
    ghcr.io/micpreuss/anchormap:0.1.2 \
    plot_anchors --config configs/your_run_plots.yaml --out-dir results/your_run/figures
```

Nextflow is **not** how you run AnchorMap — there is no DAG to orchestrate. The [`nextflow/`](../nextflow/)
harness exists only to *validate that this image runs flawlessly under Nextflow* (the lab's container
contract: no exit-126, GCS-FUSE write perms, working `ps`/trace). See [`nextflow/README` / main.nf](../nextflow/main.nf).

## Pinning rationale

- **Base `rocker/r-ver:4.6.0`** — pinned to **the same R the engine was validated on** (the host's
  R 4.6.0). This intentionally departs from the ADD §7.3 `rocker/r-ver:4.4.2` "parent-parity" pin: that
  parity reason was GenomicSEM/coloc image alignment, which doesn't apply to AnchorMap (it omits
  GenomicSEM), and a 4.4.2 base would ship CRAN packages *older* than the validated set — which is what
  surfaced a `future.apply` z-sweep regression during the build. Matching the validated R removes the
  skew and lets a single snapshot reproduce the validated `future.apply`/`ggplot2` versions.
  Debian/Ubuntu-based ⇒ `procps`/`ps` available for Nextflow trace.
- **CRAN via one dated Posit P3M snapshot** (`--build-arg P3M_SNAPSHOT=YYYY-MM-DD`, default `2026-06-01`) —
  every CRAN package (including the bootstrap `remotes`, installed after the snapshot repo is wired in)
  resolves to the version current at that date. On the 4.6.0 base this yields the validated
  `future.apply 1.20.x` and `ggplot2 4.0.x` from a single snapshot (no split pin). The codename is
  detected from the base image's `/etc/os-release` so the binary URL is always correct; the verify step
  prints the resolved versions.
- **Reproducibility scope.** R and every CRAN package are *version-pinned* (base tag + snapshot date), so
  the R layer is deterministic. The base **tag** and its apt packages are not digest-pinned by default,
  so the OS layer can drift as rocker rebuilds `4.6.0` upstream. For a fully byte-reproducible release,
  pin the resolved base digest: `--build-arg BASE_IMAGE=rocker/r-ver:4.6.0@sha256:<digest>` (resolve it
  once with `docker buildx imagetools inspect rocker/r-ver:4.6.0`).
- Version tags such as `0.1.2` are the run interface; GHCR also carries `latest` as a convenient pointer
  to the newest release.

## Container fixes (carried from the parent)

Same three fixes the parent project's `postgwas` image needs to run under Nextflow on Google Batch
(see [`../../UKBB_CLUSTER_GWAS/docker/postgwas/README.md`](../../UKBB_CLUSTER_GWAS/docker/postgwas/README.md)):

1. **`procps`** (apt) — Nextflow's trace/metrics poll `ps`; without it trace is empty.
2. **`USER root`** — write access to GCS-FUSE-mounted work dirs on Google Batch (non-root users hit
   permission-denied on file operations).
3. **`ENTRYPOINT []`** — prevents any upstream image entrypoint from intercepting Nextflow's bash
   wrapper (`.command.run`), which otherwise exits 126 before the script runs. (rocker has no
   problematic entrypoint; this is belt-and-suspenders, asserted by the Nextflow smoke.)

## What is deliberately NOT in the image

- **GenomicSEM.** The ADD §7.3 lists a pinned `remotes::install_github("GenomicSEM/GenomicSEM@<commit>")`,
  but the engine **never loads it**: [`R/ingest_rds.R`](../R/ingest_rds.R) reads the `ldsc()` `.rds`
  (a plain named list `$S`/`$V`/`$I`) with base `readRDS()`. AnchorMap *consumes* the artifact; it does
  not *run* `ldsc()` (explicitly out of scope, ADD §5). Omitting it removes a heavy GitHub install + its
  dependency tree + a pinning surface. Reintroduce it only if a future stage runs `ldsc()` inside AnchorMap.
- **`argparse` / `optparse`.** Both CLIs hand-roll their argument parsing (`R/cli.R`).

## Build-time self-tests (the build fails on regression)

The build runs the self-contained synthetic fixture (`inst/fixtures/synthetic_ldsc_panel.rds`, no
external data) twice:

1. **Engine** — must recover `C5_sub0 → anthro [sharp]` and a `FINISHED ok` log.
2. **Figures** — the `ggplot2`/`ragg`/`cairo` stack must render `anchor_lollipops_synthetic.png`.

So a bad dependency, a broken engine, or a `ggplot2` that can't draw the figures all fail
`docker build`, not a downstream user.

### Why the R 4.6.0 base + 2026 snapshot (matching the validated env)

The engine's parallel z-sweep ([`R/sensitivity.R`](../R/sensitivity.R)) calls `future_lapply(...,
future.globals = FALSE)` and relies on the worker closure carrying `cfg`. **`future.apply 1.11.x` has a
regression that drops that closure global** (`object 'cfg' not found`, failing even at `--threads 1`);
`1.20.x` — the version the engine was validated against — fixes it. An earlier attempt pinned the base
to `rocker/r-ver:4.4.2` (ADD §7.3) + a `2025-03-01` snapshot, which shipped `future.apply 1.11.3` and
`ggplot2 3.5.x` — i.e. CRAN packages *older* than the validated set, which is exactly what exposed the
regression. Rather than carry a two-date split pin onto an old R, the image now pins the base to
**R 4.6.0 (the host's validated R)** and uses **one current snapshot**, so `future.apply 1.20.x` and
`ggplot2 4.0.x` both come from a single coherent source matching what Phases 1–4 were validated on. The
engine self-test runs the sweep at `--threads 2`, so a `future.apply` regression fails the build.

## Build and publish

Tagged releases are published automatically to GitHub Container Registry by
[`publish-container.yml`](../.github/workflows/publish-container.yml). The version tag must match the
package version in `DESCRIPTION`:

```bash
git tag -a v0.1.2 -m "AnchorMap 0.1.2"
git push origin v0.1.2
```

The workflow builds `linux/amd64`, runs both Dockerfile self-tests, and publishes
`ghcr.io/micpreuss/anchormap:0.1.2` and `ghcr.io/micpreuss/anchormap:latest`. GitHub Actions uses its
built-in `GITHUB_TOKEN`; no PAT or repository secret is needed. After the package is published for
the first time, set its visibility to **Public** in the package settings so it can be pulled without
authentication. GitHub warns that making a package public cannot be reversed.

For a manual publication (for example, before enabling Actions), build and push with a classic PAT
carrying `write:packages`:

```bash
# Authenticate without putting the token in shell history.
read -rsp "GHCR token: " GHCR_TOKEN; echo
printf '%s' "$GHCR_TOKEN" | docker login ghcr.io -u micpreuss --password-stdin
unset GHCR_TOKEN

# Build for the release platform from the repository root.
docker build --platform linux/amd64 \
    -t ghcr.io/micpreuss/anchormap:0.1.2 \
    -f docker/Dockerfile .

docker push ghcr.io/micpreuss/anchormap:0.1.2
```

Bump the tag (never reuse `0.1.2`) when the engine or a pinned dependency changes; update the tag in
[`nextflow/params/gcp.yaml`](../nextflow/params/gcp.yaml) and any consuming configs.
