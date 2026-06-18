---
name: create-rules
description: Create global rules (CLAUDE.md) for a scientific data/pipeline project by analyzing the repo
---

# Create Global Rules (scientific data pipelines)

Generate a `CLAUDE.md` that orients Claude on a **scientific computing / data-analysis** repo:
what the science is, how data flows, where compute runs, which schemas are load-bearing, and
how to run and validate the analysis.

This is **not** a web-app/library generator. The unit of work here is a *pipeline stage* or an
*analysis*, not a "feature"; the artifacts are *datasets, sumstats, matrices, manifests, figures*,
not "API endpoints". Classify and describe accordingly.

---

## Objective

Produce project-specific rules that tell Claude:
- **What** the project studies and what question each stage answers.
- **What data** it operates on — formats, schemas, and what level of claim each representation permits.
- **Where** compute runs (laptop / HPC scheduler / cloud batch) and where data lives (in-repo / object store / scratch).
- **Which conventions** are load-bearing (column schemas, units, sign preservation, config-over-CLI).
- **How** to run a stage and **how** to validate it (there is rarely a unit-test suite — validation is positive controls, smoke runs, and schema checks).

---

## Phase 1: DISCOVER — classify the project along five axes

Do **not** collapse a science repo into one "type". Classify it along five orthogonal axes; the
combination is the project type. Derive each axis by **reading signature files**, not by guessing
from a folder name.

### Axis A — Orchestration (how stages are wired together)

| Value | Signature indicators |
|---|---|
| Nextflow | `*.nf`, `main.nf`, `nextflow.config`, `-profile`, `params/*.json` |
| Snakemake | `Snakefile`, `workflow/rules/*.smk`, `config/config.yaml` |
| WDL / CWL | `*.wdl` + `inputs.json`; `*.cwl` |
| Make / shell DAG | `Makefile` targets that chain data steps; numbered `0x_*.sh` scripts |
| Notebook-driven | `*.ipynb` executed in order; `papermill`/`jupyter nbconvert` |
| Ad-hoc scripts | standalone `scripts/*.py` / `*.R` run by hand, no orchestrator |

### Axis B — Compute backend (where heavy steps execute)

| Value | Signature indicators |
|---|---|
| Local | runs on the laptop/workstation; no scheduler config |
| HPC scheduler | SLURM (`sbatch`, `#SBATCH`), SGE/UGE (`qsub`), PBS/LSF |
| Cloud batch | Google Batch (`google-batch` profile), AWS Batch, Azure Batch |
| Kubernetes | `k8s`/`eks`/`gke` executor, Helm/Argo |
| Mixed | heavy compute on a backend **+** local post-processing/plotting (very common) |

### Axis C — Analysis domain (what science / data)

| Value | Signature indicators |
|---|---|
| Statistical genetics / GWAS | GWAS sumstats (`*.tsv.gz`), PLINK/bgen, LDSC, GenomicSEM, MAGMA, FUMA, COJO, coloc/SuSiE, HapMap/1000G refs |
| Sequencing / bioinformatics | FASTQ/BAM/CRAM/VCF, nf-core, reference genome + GTF, samplesheets, aligners (STAR/bwa) |
| Single-cell / spatial | `*.h5ad`/`*.h5`/`*.rds`, scanpy/Seurat/scvi-tools, AnnData, cell × gene matrices |
| Imaging | TIFF/CZI/OME/DICOM, ImageJ/CellProfiler/napari, segmentation masks |
| Instrument / wet-lab data | instrument exports (CSV/XLSX/PDF), LIMS/ELN, plate maps, Allotrope/ASM |
| Tabular / epi / clinical | cohort tables, REDCap, survival/regression, `data dictionary` |
| Modeling / ML | training scripts, checkpoints, `configs/*.yaml`, metrics logs, sweeps |
| Simulation / numerics | parameter sweeps, solver configs, HDF5/NetCDF outputs |
| Method / tool package | a reusable library: `pyproject.toml` with `project.scripts`, or R `DESCRIPTION`/`NAMESPACE` |

> A repo can carry two domains (e.g. a GWAS pipeline that also vendors a clustering *tool*). Name both.

### Axis D — Data scale & locality (where the real data is)

| Value | Signature indicators |
|---|---|
| In-repo small | committed `data/`, test fixtures, `< ~100 MB` |
| External object store | `gs://`, `s3://`, `az://` paths in configs/scripts |
| HPC filesystem | absolute `/scratch`, `/projects`, `/work` paths; symlinked data dirs |
| Public/remote download | accessions (GEO `GSE…`, SRA `SRR…`, dbGaP, EGA), download scripts |

### Axis E — Reproducibility layer (envs & containers)

| Value | Signature indicators |
|---|---|
| Containers | `Dockerfile`, `*.def`/`Singularity`, image refs in configs (Artifact Registry / Docker Hub / quay) |
| Conda / mamba | `environment.yml`, `meta.yaml`, `micromamba` |
| Python venv | `requirements.txt`, `pyproject.toml`, `uv.lock`/`poetry.lock` |
| R | `renv.lock`, `DESCRIPTION`, `.Rprofile` |
| None pinned | bare interpreters; flag this as a reproducibility gap in the output |

### Detection sweep (run these, then read what they surface)

```bash
ls -la
# orchestration + repro signatures
find . -maxdepth 3 \( -name '*.nf' -o -name 'nextflow.config' -o -name 'Snakefile' \
  -o -name '*.smk' -o -name '*.wdl' -o -name '*.cwl' -o -name 'Makefile' \
  -o -name 'Dockerfile' -o -name '*.def' -o -name 'environment.yml' \
  -o -name 'requirements.txt' -o -name 'pyproject.toml' -o -name 'renv.lock' \
  -o -name 'DESCRIPTION' \) -not -path '*/.git/*' 2>/dev/null
# scheduler / backend hints
grep -rIl -E '#SBATCH|sbatch|qsub|google-batch|aws-batch|s3://|gs://|az://' . \
  --include='*.sh' --include='*.config' --include='*.nf' --include='*.smk' \
  --include='*.json' --include='*.yaml' 2>/dev/null | head
# data domain hints (extensions actually present)
find . -type f -not -path '*/.git/*' | sed -E 's/.*(\.[a-z0-9]+(\.gz)?)$/\1/' \
  | sort | uniq -c | sort -rn | head -30
# the existing entry points and any workflow index
find . -iname 'README*' -not -path '*/.git/*' | head
git log -1 --format='%ai  %s' 2>/dev/null
```

### Worked classification (example, end-to-end)

A repo that clusters UK Biobank phenotypes by GWAS signature then runs gene/disease follow-ups
classifies as:

- **A:** Nextflow (`scripts/*/main.nf`, `nextflow.config`, `-profile gcp`) **+** ad-hoc local Python/R for post-processing → *orchestrated pipeline with a local post-processing tail*.
- **B:** Cloud batch (Google Batch, head VM submits) **+** local — **Mixed**.
- **C:** Statistical genetics / GWAS (LDSC, GenomicSEM gPCA, MAGMA, coloc) **and** a vendored clustering *method package* (DIMPLE-GWAS). Two domains.
- **D:** External object store (`gs://…-usc1/`); `data/` holds only small refs.
- **E:** Containers (Artifact Registry images, pinned `postgwas:1.3`) + per-stage configs.

That five-tuple — not a single label — is the "project type". Lead the generated `CLAUDE.md`
with it.

---

## Phase 2: ANALYZE — extract the stack, the schemas, and the key files

### Scientific stack (fill what the repo actually uses)

- **Orchestration & profiles** — engine + profile names (e.g. `test` local vs `gcp`/`slurm`).
- **Compute** — backend, region/partition, how a run is submitted (head node? login node? laptop?).
- **Containers / envs** — image registry + tags (note any **pinned** image and why), conda/renv locks.
- **Storage** — bucket/scratch root and its subtree layout; how `workDir`/`outdir` are derived.
- **Languages & libraries** — Python (pandas/numpy/scipy/…), R (which Bioconductor/CRAN), Bash; versions if pinned.
- **Methods** — the statistical/algorithmic methods by name (clustering, gPCA, LDSC `rg`, MAGMA, SuSiE-coloc, MAD filtering, …). This is what a returning reader needs most.
- **External reference data** — reference panels, ontologies, public cohorts (HapMap3, 1000G, GTEx, FinnGen, …).

### Conventions — the load-bearing part for science

Read code and headers; do **not** infer schemas from filenames.

- **Data schemas** — for each interchange file, the exact column contract. Example precision to aim for:
  `cluster manifest = (phenotype_id, cluster_label, vcf_path)`;
  `gPCA sumstats = (SNP, EA, OA, EAF, BETA, SE, P, N_eff, N_obs)` where `N_obs` is a *fractional weight, not a sample size*.
  Capturing "what each column **is**, and what claim it licenses" is the whole point.
- **Units & sign** — log vs linear, `-log10(p)` vs `p`, OR vs beta, build GRCh37 vs 38, and any
  **sign-preservation** requirement (`b_SMR`, `rg`, `beta` signs that must survive aggregation).
- **Parameterization** — config/params files vs CLI flags; where the *canonical* run params live.
- **Naming** — dataset tags, run tags, `_with_features`/`_test` suffixes and what they mean.
- **Gotchas** — known data bugs, misdetections, crash modes (e.g. "`n_eff_mode=sum` → indefinite
  `cov_z` → MAGMA crash; use `max`"). One line each; these save hours.

### Key files to name in the output

- Pipeline **entry points** (`main.nf`, `Snakefile`, top driver scripts).
- **Canonical params/configs** (the one run that is "the answer").
- **Data dictionary** (IDP key, sample sheet, ontology maps).
- **Workflow-index README** (the top-level README that lists stages — see the `create-readme` skill).
- **Env/container** definitions and any handover/known-issue notes.

---

## Phase 3: GENERATE — write CLAUDE.md

**Output path:** `CLAUDE.md` (repo root).

Use the template below verbatim as the skeleton, then fill from Phases 1–2. Drop any section the
repo genuinely lacks (say so rather than leaving a stub). Keep it terse and mechanistic — match the
voice of the repo's own READMEs.

````markdown
# <project-name>

<1–3 sentences: the scientific goal and the shape of the work. State whether it's greenfield or
brownfield (most analytical work done vs. actively building), and who the reader is.>

---

## Project type

<The five-axis classification as prose or a short list:>
- **Orchestration:** <Nextflow | Snakemake | … (+ local tail?)>
- **Compute backend:** <local | SLURM | Google Batch | mixed …>
- **Domain:** <stat-genetics | single-cell | imaging | … (+ any vendored method package)>
- **Data locality:** <in-repo | gs:// | /scratch | accessions>
- **Reproducibility:** <containers + conda/renv; note pinned images>

<One line on build/test posture, e.g. "No build system or unit-test suite at root; each stage is a
self-contained subdir with its own config, params, and scripts.">

---

## Tech stack

| Layer | Tools |
|---|---|
| Orchestration | <engine + profiles> |
| Compute | <backend, region/partition, submission host> |
| Containers / envs | <registry:tags, conda/renv locks; note any pinned image + why> |
| Storage | <bucket/scratch root + subtree layout> |
| Languages | <Python libs, R libs, Bash> |
| Methods | <named statistical/algorithmic methods> |
| External reference data | <panels, ontologies, public cohorts> |

---

## Repository structure

```
<tree of the meaningful dirs with a one-line role each; mark the entry point "start here">
```

<Point at the workflow-index README: "Read the relevant stage README before editing that stage.">

---

## Data flow (high-level)

```
<ascii diagram: stageA ──► stageB ──► stageC ; show forks/joins and terminal stages>
```

<2–4 bullets on what flows between stages and what must NOT be crossed (e.g. reporting-only
variants that must never feed a pipeline).>

---

## Canonical datasets / runs

| Dataset | Canonical run / solution | Notes |
|---|---|---|
| <name> | <the params/config that is "the answer"> | <why this one; gotchas> |

---

## Conventions

### Compute / infrastructure
- <where to submit from; profile names; spot/preemption; output convention>

### Data schemas
- **<file>:** `(col, col, …)` — <what each non-obvious column IS and what claim it licenses>
- <units / sign / build conventions; sign-preservation requirements>

### Code style
- <config-over-CLI threshold; read real headers before parsing; long vs wide format rules; language idioms>

---

## Running stages

```bash
# Canonical, copy-pasteable invocation per major stage
<command>
```

---

## Validation / testing

<There is usually no unit-test suite. Document the real validation:>
- Smoke runs with small `--test` params.
- Pre-registered **positive controls** (e.g. "coloc gates on HMGCR↔LDL PP.H4 > 0.95").
- Schema / sanity checks (`N_eff` positive and physically plausible; column contracts hold).

---

## Key files to know

- <entry points, canonical params, data dictionary, workflow-index README, env/container defs, handover notes>

---

## Working with this repo via Claude

- **READMEs:** use the `create-readme` skill — one top-level index that lists each subproject as a
  linked list item, plus a canonical per-subproject README (Orientation / Inputs / Outputs /
  Workflow / How to run / Gotchas / Results / Related).
- **Finishing a subproject:** use the `report-findings` skill to write a `REPORT.md` (results,
  findings, and a clearly non-binding exploratory interpretation) and link it from the subproject
  README's Results section.
- **Orientation:** `prime`. **Committing:** `commit`.

---

## Out-of-scope / external links

- <external web tools with no API (FUMA, SMR portal, COJO), sibling repos, accession sources>
````

> **Encode the README + report workflow in the generated file.** The "Working with this repo via
> Claude" section above is mandatory — it is how a future session discovers the `create-readme` and
> `report-findings` skills.

---

## Phase 4: OUTPUT

```markdown
## Global Rules Created

**File:** `CLAUDE.md`

### Project type (five-axis)
- Orchestration / Compute / Domain / Data locality / Reproducibility: <one line each>

### Stack summary
<key tools + methods detected>

### Load-bearing schemas captured
<the interchange-file column contracts written into CLAUDE.md>

### Gaps / open questions
<anything you could not derive from the repo — unpinned envs, unclear canonical run, missing data dictionary>

### Next steps
1. Review `CLAUDE.md`; correct any schema column you had to infer.
2. Run `create-readme` for the top-level index + each subproject.
3. Confirm the canonical run/params table.
```

---

## Tips

- **Object before method.** Lead each section with the data object, then the method over it. A
  reader needs "what is `N_obs`" before "how gPCA computes it".
- **Schemas are the highest-value content.** A wrong column contract causes silent, costly errors;
  capture them exactly, from code, not filenames.
- **Don't claim status without evidence.** "complete / active / scaffolded" must come from prose or
  `git log`, not vibes.
- **Keep it scannable.** Link to stage READMEs and memory notes for depth; don't duplicate them.
- **Flag reproducibility gaps** (unpinned envs, uncontainerized steps) explicitly — that is useful
  signal, not a failure of the skill.
