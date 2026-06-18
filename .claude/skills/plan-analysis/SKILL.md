---
name: plan-analysis
description: "Create a comprehensive, context-rich plan for a new analysis or pipeline stage through repo + data + methods research"
---

# Plan an analysis / pipeline stage

## Analysis: $ARGUMENTS

## Mission

Turn an analysis request into a **complete implementation plan** through systematic inspection of
the repo, the data and its schemas, and the relevant methods/literature.

**Core principle**: we do NOT run the analysis or write production code in this phase. The goal is a
context-rich plan that lets an execution agent build the stage in one pass.

**Key philosophy**: *Context is king.* The plan must carry everything the executor needs — upstream
schemas, the exact method, reference-data versions, the canonical params, validation controls, and
runnable validation commands — so it succeeds on the first attempt without re-deriving context.

The unit here is an **analysis** or **pipeline stage**, not a software "feature". Artifacts are
datasets, sumstats, matrices, manifests, figures — not API endpoints.

## Planning process

### Phase 1: Frame the analysis

- State the **scientific question** and the inference the output licenses.
- Name the **object**: what data it operates on, and what representation (matrix? long table? VCF?).
- Classify: New stage / Extension of a stage / Re-analysis / Method swap / Bug fix.
- Assess complexity: Low / Medium / High, and the main risk (data, method, or compute).

Write it as:

```
As an analyst
I want to <analysis/transformation>
So that <inference / decision it supports>
```

### Phase 2: Repo + data intelligence

**1. Pipeline position**
- Where does this stage sit in the DAG? Upstream producer(s), downstream consumer(s).
- Read the upstream stage README + the workflow-index README. Mirror the orchestration pattern
  (Nextflow `main.nf` / Snakefile / driver script) used by sibling stages.

**2. Schema contracts (the load-bearing part)**
- For every input, read the **actual file header / code** — never infer columns from a filename.
- Record the exact column contract and what each non-obvious column *is* (e.g. `N_obs` = fractional
  gPCA weight, not a sample size; `p` may be `-log10`). Note units, sign conventions, and genome build.
- Define the **output schema** the downstream consumer expects.

**3. Reference data & methods**
- Catalog reference panels / ontologies / external cohorts the stage needs, **with versions**.
- Name the statistical/algorithmic method precisely; find its canonical implementation in the repo
  or the library docs. Note assumptions (e.g. LDSC `rg` observed vs liability scale) and failure modes.

**4. Compute & reproducibility**
- Backend (local / SLURM / cloud batch), profile, container/env, where it submits from.
- Canonical params/config location and naming.

**Clarify ambiguities now.** If the method, the input run, or the success criterion is unclear, ask
before proceeding.

### Phase 3: Methods / literature research

- Confirm the method's correct usage and current best practice; capture doc links with anchors.
- Note known pitfalls (batch effects, double-counting samples, multiple-testing, indefinite
  covariance, power) relevant to this analysis.
- Identify **positive/negative controls** the field uses to sanity-check this kind of result.

```markdown
## Relevant references
- [Method docs](https://…#section) — Why: exact call + assumptions for this stage
- [Paper / vignette](https://…) — Why: control expectation (e.g. HMGCR↔LDL should colocalize)
```

### Phase 4: Strategic thinking

- How does this stage fit the existing DAG and conventions?
- Critical dependencies and order of operations.
- What could go wrong scientifically (confounds, leakage, scale/locality) and computationally
  (preemption, memory, indefinite matrices)?
- How is correctness *demonstrated* — which controls, which schema checks?

### Phase 5: Write the plan

Fill this template for the execution agent:

````markdown
# Analysis: <name>

Validate the schema contracts and method usage against the actual code/data before implementing.
Pay special attention to column names, units, sign, and genome build.

## Question & object
<The scientific question; the data object and representation; the inference the output licenses.>

## Analyst story
As an analyst / I want to <analysis> / So that <inference>.

## Pipeline position
- **Upstream:** <stage(s) + the input run/manifest>
- **Downstream:** <consumer stage(s)>  ·  **Orchestration pattern to mirror:** <file:lines>

## CONTEXT REFERENCES — READ BEFORE IMPLEMENTING

### Input schemas (verified from code/headers)
- `<path>` — columns `(…)`; <what each non-obvious column is; units/sign/build>
- ...

### Output schema (contract for downstream)
- `<path>` — columns `(…)`; <semantics>

### Reference data & methods
- <panel/ontology/cohort> — version <…>, location <path/gs://…>
- Method: <name> — canonical impl `<file>` / [docs](…#anchor); assumptions: <…>; failure modes: <…>

### Files to read / create
- READ: `<file>` (lines) — Why: <pattern/schema to mirror>
- CREATE: `<file>` — <role>

## METHOD / IMPLEMENTATION PLAN

### Phase 1: Inputs & harmonization
<prepare/validate inputs; assert schema; align build/units>

### Phase 2: Core analysis
<the method step(s); params; per-unit processing>

### Phase 3: Outputs & integration
<write outputs in the contract schema; wire into the DAG / index README>

### Phase 4: Validation
<controls, smoke run, schema/sanity checks>

## STEP-BY-STEP TASKS (execute top to bottom; each atomic + checkable)

Keywords: CREATE / UPDATE / ADD / MIRROR / RUN.

### {ACTION} {target}
- **IMPLEMENT**: <detail>
- **PATTERN**: <file:line to mirror>
- **DATA/SCHEMA**: <input contract asserted; output contract produced>
- **GOTCHA**: <known pitfall to avoid>
- **VALIDATE**: `<runnable command>`

## VALIDATION STRATEGY
There is usually no unit-test suite. Validate by:
- **Smoke run** with a tiny `--test` params variant.
- **Positive/negative controls** (state the expected result, e.g. control rg > 0.6, FTO null expected).
- **Schema/sanity checks** (column contract holds; `N_eff` positive & physically plausible; no NaN explosion).
- **Sensitivity** (if a transform/threshold is involved, re-run under an alternative and compare).

## VALIDATION COMMANDS (run all; zero schema/control failures)
```bash
<orchestrator dry-run: nextflow … -preview / snakemake -n>
<smoke run command>
<schema/control check command>
```

## ACCEPTANCE CRITERIA
- [ ] Inputs assert their schema; outputs match the downstream contract
- [ ] Positive controls recovered; negative controls null
- [ ] Smoke run passes; full run command documented
- [ ] Stage wired into the DAG and anchored in the index README
- [ ] Provenance captured (params, reference-data versions, container/env, commit)

## NOTES
<design decisions, trade-offs, scale/locality considerations>
````

## Output format

**Filename**: `.agents/plans/{kebab-case-name}.md` (create the dir if absent).
Examples: `add-coloc-bridge-stage.md`, `reharmonise-finngen-r12.md`, `sensitivity-rint-vs-log.md`.

## Quality criteria

- **Schema completeness** — every input/output column contract is explicit and verified from code.
- **Control-anchored** — at least one positive (and ideally one negative) control with an expected result.
- **Runnable validation** — every task has an executable check; the plan passes the "no prior
  knowledge" test (someone new to the repo could execute it).
- **Convention-consistent** — mirrors the repo's orchestration, storage, and naming patterns.

## Report

After writing the plan, return: summary of the analysis + approach; full path to the plan; complexity
and main risk; and a confidence score (#/10) for one-pass execution.
