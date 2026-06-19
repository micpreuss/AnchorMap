# AnchorMap — pinned, reproducible image (THE tool). Phase 5.
#
# The image IS AnchorMap: the validated Phase 1-4 R engine + figures, baked in, self-validating at
# build time. Primary run interface:
#   docker run --rm -v "$PWD:/work" -w /work anchormap:0.1.0 \
#       Rscript /opt/anchormap/anchor_map.R --config <your.yaml> [--threads N]
# Figures:
#   ... Rscript /opt/anchormap/R/plot_anchors.R --config <plots.yaml>
#
# Pinning: rocker/r-ver:4.6.0 (= the host R 4.6.0 the engine was validated on) + a single dated Posit
# P3M CRAN snapshot, so R, future.apply, and ggplot2 all match the validated versions. GenomicSEM is
# NOT installed — the engine reads ldsc() .rds artifacts with base readRDS() and never runs ldsc()
# (see docker/README.md). Carries the parent container fixes: procps / USER root / ENTRYPOINT [].
#
# Build (release, for Google Batch):  docker build --platform linux/amd64 -t anchormap:0.1.0 -f docker/Dockerfile .
# Build (local validation, native):   docker build -t anchormap:0.1.0 -f docker/Dockerfile .

FROM rocker/r-ver:4.6.0

# ---- system deps: procps (Nextflow ps/trace) + build libs + figure-render libs ----
RUN apt-get update && apt-get install -y --no-install-recommends \
      procps \
      libcurl4-openssl-dev libssl-dev libxml2-dev \
      libcairo2-dev libfreetype6-dev libfontconfig1-dev \
      libharfbuzz-dev libfribidi-dev \
      libpng-dev libtiff5-dev libjpeg-dev \
      cmake git ca-certificates locales \
    && rm -rf /var/lib/apt/lists/*
# NB: ragg -> textshaping needs harfbuzz + fribidi (hb-ft.h); systemfonts needs fontconfig.
# Omitting these is the one non-obvious figure-stack build failure on a clean rocker base.

# ---- pin CRAN to a single dated P3M snapshot (binary install; codename detected from the base) ----
# Override the date with --build-arg P3M_SNAPSHOT=YYYY-MM-DD. 2026-06-01 reproduces the host's
# validated package set on R 4.6.0 — future.apply 1.20.x (the z-sweep needs >= 1.20: 1.11.x has a
# future.globals=FALSE regression that drops the worker closure's `cfg`) and ggplot2 4.0.x (so the
# figures match the validated ones). The verify + self-test steps below confirm the resolved versions.
ARG P3M_SNAPSHOT=2026-06-01
RUN CODENAME="$(. /etc/os-release && echo "$VERSION_CODENAME")" && \
    printf '\noptions(repos = c(CRAN = "https://packagemanager.posit.co/cran/__linux__/%s/%s"))\noptions(HTTPUserAgent = sprintf("R/%%s R (%%s)", getRversion(), paste(getRversion(), R.version$platform, R.version$arch, R.version$os)))\n' \
      "$CODENAME" "$P3M_SNAPSHOT" >> "$(R RHOME)/etc/Rprofile.site"

# ---- R deps (engine + figures); NO GenomicSEM, NO argparse/optparse ----
RUN R -e "install.packages(c('data.table','yaml','poolr','Matrix','future','future.apply', \
                             'ggplot2','patchwork','scales','ggrepel','ragg'))"

# ---- verify the full runtime set loads ----
RUN R -e "pk <- c('data.table','yaml','poolr','Matrix','future','future.apply','ggplot2', \
                  'patchwork','scales','ggrepel','ragg'); \
          ok <- vapply(pk, requireNamespace, logical(1), quietly=TRUE); \
          print(data.frame(pkg=pk, ver=vapply(pk, function(p) as.character(packageVersion(p)), ''), ok=ok)); \
          if (!all(ok)) quit(status=1)"

# ---- bake the tool in (pinned with the image tag) ----
WORKDIR /opt/anchormap
COPY R/ ./R/
COPY anchor_map.R ./anchor_map.R
COPY configs/ ./configs/
COPY ontology/ ./ontology/
COPY tests/fixtures/synthetic_ldsc_panel.rds tests/fixtures/synthetic_panel_ontology.tsv ./tests/fixtures/

# ---- build-time self-test #1 (engine): self-contained synthetic .rds run must recover ----
#      C5_sub0 -> anthro [sharp] + a FINISHED log. Fails the build on any engine/dep regression.
RUN Rscript anchor_map.R --config configs/synthetic_rds.yaml --threads 2 \
 && grep -Eq 'C5_sub0.*anthro.*sharp' results/synthetic_rds/cluster_anchor_labels.tsv \
 && grep -q '^FINISHED ok' results/synthetic_rds/anchormap.log

# ---- build-time self-test #2 (figures): the ggplot/ragg/cairo stack must render a PNG ----
#      from those TSVs. Fails the build if the pinned ggplot2 can't draw the figures.
RUN Rscript R/plot_anchors.R --config configs/synthetic_rds_plots.yaml \
 && ls results/synthetic_rds/figures/anchor_lollipops_synthetic.png \
 && rm -rf results/synthetic_rds

# ---- project container fixes (see ../UKBB_CLUSTER_GWAS/docker/postgwas/README.md) ----
# Docker has NO inline-comment syntax: anything after the instruction is parsed as its argument, so
# these comments MUST be on their own lines (an inline `# ...` silently mangles USER and ENTRYPOINT).
# USER root  -> write to GCS-FUSE-mounted work dirs on Google Batch.
USER root
# ENTRYPOINT []  -> keep it empty (exec form) so an upstream entrypoint can't intercept Nextflow's
# .command.run wrapper (which would exit 126). Must be `[]`, never `/bin/sh -c "[]"`.
ENTRYPOINT []
