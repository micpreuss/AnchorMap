#!/usr/bin/env nextflow

// AnchorMap — Nextflow CONTAINER-VALIDATION HARNESS (Phase 5). NOT an orchestration layer.
//
// AnchorMap proper is the R engine + the pinned Docker image (run via `docker run … Rscript
// anchor_map.R --config …`). There is no DAG to orchestrate — the engine is one fast single-process
// job whose z-sweep/perm_p parallelism is internal (future/setDTthreads). This harness's ONLY job is
// to prove the image runs flawlessly under Nextflow's execution model, i.e. the three Nextflow-
// specific container failure modes a plain `docker run` can't catch:
//   ENTRYPOINT []  -> no exit-126 when Nextflow's .command.run wrapper invokes the image
//   procps/ps      -> Nextflow trace/metrics can poll the process
//   USER root      -> the engine can write the (GCS-FUSE) work dir   [only truly exercised on -profile gcp]
//
// Profiles (see nextflow.config):
//   -profile test  local docker, the CI gate (validates entrypoint + procps + output capture)
//   -profile gcp   Google Batch (spot), the one-time check that exercises USER-root / GCS-FUSE writes
//
// Payload = the image's self-contained synthetic .rds fixture (no external data). We only assert the
// engine produced its outputs + a FINISHED log where Nextflow captures them — values are validated by
// the engine's own Phase 1-3 tests, not here.

nextflow.enable.dsl = 2

params.anchormap_container = null   // required: image under test (set per-profile params file)
params.smoke_cpus          = 2

process ANCHORMAP_SMOKE {
    container params.anchormap_container
    cpus   params.smoke_cpus
    memory '4 GB'

    output:
    path "*.tsv",        emit: tsv
    path "anchormap.log", emit: log

    script:
    """
    export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
           VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1

    # Point the baked synthetic config at the absolute in-image fixtures and the ABSOLUTE work-dir as
    # out_dir, so the ENGINE itself writes into Nextflow's work dir (on Batch that dir is the GCS-FUSE
    # mount -> this is the real USER-root/FUSE write test). resolve_path() keeps absolute paths verbatim.
    Rscript -e '
      cfg <- yaml::read_yaml("/opt/anchormap/configs/synthetic_rds.yaml")
      cfg\$rds      <- "/opt/anchormap/tests/fixtures/synthetic_ldsc_panel.rds"
      cfg\$ontology <- "/opt/anchormap/tests/fixtures/synthetic_panel_ontology.tsv"
      cfg\$out_dir  <- getwd()
      yaml::write_yaml(cfg, "run_config.yaml")'

    Rscript /opt/anchormap/anchor_map.R --config run_config.yaml --threads ${task.cpus}

    # Fail the task (and the harness) if the engine didn't finish cleanly.
    grep -q '^FINISHED ok' anchormap.log
    """
}

workflow {
    main:
    if (!params.anchormap_container)
        error "Missing required param: --anchormap_container (set via -params-file nextflow/params/<profile>.yaml)"
    ANCHORMAP_SMOKE()

    publish:
    smoke_tsv = ANCHORMAP_SMOKE.out.tsv
    smoke_log = ANCHORMAP_SMOKE.out.log
}

output {
    smoke_tsv {
        path 'smoke'
        mode 'copy'
    }
    smoke_log {
        path 'smoke'
        mode 'copy'
    }
}
