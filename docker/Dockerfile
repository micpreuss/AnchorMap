# AnchorMap — pinned, version-constrained image (THE tool).
#
# The image IS AnchorMap: the `anchormap` R package (engine + figures) installed from source, plus
# `anchor_map` / `plot_anchors` PATH shims, self-validating at build time. Primary run interface:
#   docker run --rm -v "$PWD:/work" -w /work anchormap:0.1.2 \
#       anchor_map --config <your.yaml> --out-dir results --threads N
# Figures:
#   ... plot_anchors --config <plots.yaml> --out-dir results/figures
#
# Pinning: rocker/r-ver:4.6.0 (= the host R 4.6.0 the engine was validated on) + a single dated Posit
# P3M CRAN snapshot, so R and every CRAN dependency (incl. the bootstrap `remotes`, installed AFTER the
# snapshot repo is wired in below) resolve to fixed, validated versions. GenomicSEM is NOT installed —
# the engine reads ldsc() .rds artifacts with base readRDS() and never runs ldsc() (see docker/README.md).
# Carries the container fixes: procps / USER root / ENTRYPOINT [].
#
# Reproducibility scope: R + every CRAN package are version-pinned (base tag + snapshot date). The base
# tag and its apt packages are NOT digest-pinned by default, so they can drift as rocker rebuilds the
# tag. For a fully byte-reproducible release, pass the resolved base digest:
#   docker build --build-arg BASE_IMAGE=rocker/r-ver:4.6.0@sha256:<digest> -t anchormap:0.1.2 -f docker/Dockerfile .
# (resolve it once with `docker buildx imagetools inspect rocker/r-ver:4.6.0`).
#
# Build (release, for Google Batch):  docker build --platform linux/amd64 -t anchormap:0.1.2 -f docker/Dockerfile .
# Build (local validation, native):   docker build -t anchormap:0.1.2 -f docker/Dockerfile .

ARG BASE_IMAGE=rocker/r-ver:4.6.0
FROM ${BASE_IMAGE}

LABEL org.opencontainers.image.title="AnchorMap" \
      org.opencontainers.image.description="Correlation-aware ontology anchoring for latent genetic cluster factors" \
      org.opencontainers.image.source="https://github.com/micpreuss/AnchorMap" \
      org.opencontainers.image.licenses="MIT"

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
# validated package set on R 4.6.0 (future.apply 1.20.x — the z-sweep needs >= 1.20 — and ggplot2
# 4.0.x so the figures match the validated ones). The self-tests below confirm the resolved set works.
ARG P3M_SNAPSHOT=2026-06-01
RUN CODENAME="$(. /etc/os-release && echo "$VERSION_CODENAME")" && \
    printf '\noptions(repos = c(CRAN = "https://packagemanager.posit.co/cran/__linux__/%s/%s"))\noptions(HTTPUserAgent = sprintf("R/%%s R (%%s)", getRversion(), paste(getRversion(), R.version$platform, R.version$arch, R.version$os)))\n' \
      "$CODENAME" "$P3M_SNAPSHOT" >> "$(R RHOME)/etc/Rprofile.site"

# ---- install the anchormap package from source (deps resolved from DESCRIPTION) ----
# NO GenomicSEM (the engine reads ldsc() .rds with base readRDS, never runs ldsc()).
# `remotes::install_deps(dependencies = TRUE)` pulls Imports (incl. optparse for the CLI) +
# Suggests (poolr cross-check, ragg figure backend) against the pinned snapshot; then R CMD INSTALL.
# `install.packages('remotes')` itself resolves from the same snapshot (the snapshot repo was wired
# into Rprofile.site above), so the bootstrap is version-pinned too.
WORKDIR /src
COPY DESCRIPTION NAMESPACE LICENSE ./
COPY R/ ./R/
COPY inst/ ./inst/
COPY man/ ./man/
RUN R -e "install.packages('remotes'); remotes::install_deps('/src', dependencies = TRUE, upgrade = 'never')" \
 && R CMD INSTALL --no-multiarch --with-keep.source /src

# ---- PATH shims so the interface is `anchor_map` / `plot_anchors` ----
COPY docker/bin/anchor_map docker/bin/plot_anchors /usr/local/bin/
RUN chmod +x /usr/local/bin/anchor_map /usr/local/bin/plot_anchors

# ---- verify the runtime set loads ----
RUN R -e "pk <- c('data.table','yaml','optparse','poolr','future','future.apply','ggplot2', \
                  'patchwork','scales','ggrepel','ragg','anchormap'); \
          ok <- vapply(pk, requireNamespace, logical(1), quietly=TRUE); \
          print(data.frame(pkg=pk, ver=vapply(pk, function(p) as.character(packageVersion(p)), ''), ok=ok)); \
          if (!all(ok)) quit(status=1)"

# ---- build-time self-test #1 (engine): the shipped synthetic .rds run must recover ----
#      C5_sub0 -> anthro [sharp] + a FINISHED log. Fails the build on any engine/dep regression.
RUN mkdir -p /tmp/st \
 && anchor_map --config synthetic_rds --out-dir /tmp/st/out --threads 2 \
 && grep -Eq 'C5_sub0.*anthro.*sharp' /tmp/st/out/cluster_anchor_labels.tsv \
 && grep -q '^FINISHED ok' /tmp/st/out/anchormap.log

# ---- build-time self-test #2 (figures): the ggplot/ragg/cairo stack must render a PNG ----
#      from those TSVs (a plots config with absolute paths into self-test #1's output).
RUN Rscript -e 'cfg <- yaml::read_yaml(system.file("configs/synthetic_rds_plots.yaml", package="anchormap")); \
                cfg$out_dir <- "/tmp/st/out/figures"; \
                cfg$tracks[[1]]$scores <- "/tmp/st/out/category_anchor_scores.tsv"; \
                cfg$tracks[[1]]$labels <- "/tmp/st/out/cluster_anchor_labels.tsv"; \
                yaml::write_yaml(cfg, "/tmp/st/plots.yaml")' \
 && plot_anchors --config /tmp/st/plots.yaml \
 && ls /tmp/st/out/figures/anchor_lollipops_synthetic.png \
 && rm -rf /tmp/st

# ---- container fixes (see docker/README.md) ----
# Docker has NO inline-comment syntax: anything after the instruction is parsed as its argument, so
# these comments MUST be on their own lines (an inline `# ...` silently mangles USER and ENTRYPOINT).
# USER root  -> write to GCS-FUSE-mounted work dirs on Google Batch.
USER root
# ENTRYPOINT []  -> keep it empty (exec form) so an upstream entrypoint can't intercept Nextflow's
# .command.run wrapper (which would exit 126). Must be `[]`, never `/bin/sh -c "[]"`.
ENTRYPOINT []
