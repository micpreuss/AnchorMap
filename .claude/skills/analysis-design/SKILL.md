---
name: analysis-design
description: Create an Analysis Design Document (ADD) for a scientific analysis or project from the conversation
argument-hint: [output-filename]
---

# analysis-design: Generate an Analysis Design Document

## Overview

Generate a comprehensive **Analysis Design Document (ADD)** from the current conversation — the
scientific analogue of a product spec. It defines the question, the data and its schemas, the
methods and assumptions, the compute/reproducibility plan, the validation controls, and the
deliverables, before any analysis is run.

Use it to scope a new analysis or project; it pairs with `plan-analysis` (which turns one stage of
the ADD into a runnable implementation plan).

## Output file

Write the ADD to: `$ARGUMENTS` (default: `ANALYSIS_DESIGN.md`).

## ADD structure

Adapt depth to available information. Use ✅ / ❌ checkboxes for in/out-of-scope items.

### 1. Overview
- The scientific **question** in 2–3 sentences; the **hypothesis** (or "exploratory, no prior").
- What inference/decision the result will support; why it matters now.

### 2. Aims
- 2–4 concrete aims (specific, assessable). One sentence each.

### 3. Background & rationale
- Prior work / what's known; the gap this analysis fills; key references.

### 4. Data
- **Sources & access:** cohorts, accessions, panels — with versions and access path (`gs://`, dbGaP, etc.).
- **Schemas:** for each dataset, the column contract and what non-obvious columns *are* (units, sign,
  genome build). Read real headers; don't infer.
- **QC & inclusion:** filters, exclusions, sample/feature counts expected.

### 5. Scope
- **In scope** (✅) and **Out of scope / deferred** (❌), grouped (Data, Methods, Compute, Deliverables).

### 6. Methods & analysis plan
- The statistical/algorithmic methods, named precisely, in execution order.
- **Assumptions** each method makes and how they're checked (e.g. normality, independence, LD,
  observed vs liability scale, batch structure).
- Estimands / test statistics; multiple-testing strategy; effect-size + uncertainty reporting.

### 7. Pipeline / compute architecture
- Stage DAG (text), orchestration engine + profile, compute backend, container/env.
- Storage layout for intermediates and outputs.

### 8. Validation & controls
- **Positive/negative controls** with expected results (e.g. HMGCR↔LDL should colocalize; a null where
  none is expected).
- Smoke tests, sensitivity analyses, and schema/sanity checks.

### 9. Outputs & deliverables
- Tables, figures, latent GWAS/sumstats, exports (FUMA/SMR/COJO), and the final `REPORT.md`
  (see the `report-findings` skill). Specify formats and the consuming tool.

### 10. Reproducibility
- Pinned envs/containers, canonical params/config, seeds, provenance capture (commit, input versions).

### 11. Success criteria
- What "done and trustworthy" means (✅ checkboxes): controls recovered, criteria met, schema contracts hold.

### 12. Phases / milestones
- 3–4 phases; each with Goal, Deliverables (✅), and Validation gate.

### 13. Risks & mitigations
- 3–5 risks specific to *this* analysis: confounding, batch effects, power, leakage/double-counting,
  data gotchas, compute/preemption — each with a mitigation.

### 14. Future / follow-ups
- Extensions and downstream analyses enabled if the aims succeed.

### 15. Appendix
- Reference links, data-dictionary pointers, related docs/repos.

## Instructions

### 1. Extract requirements
- Review the whole conversation; capture the question, data, methods, constraints, and success criteria.

### 2. Synthesize
- Organize into the sections above; fill reasonable scientific defaults where details are missing and
  flag them as assumptions; keep methods internally consistent and feasible on the stated compute.

### 3. Write the ADD
- Clear, mechanistic language; markdown (headings, tables, code blocks, checkboxes).
- Use equations/notation where they add precision over prose.
- Put exact schema contracts and control expectations in — those are the highest-value content.

### 4. Quality checks
- ✅ All sections present (or explicitly N/A)
- ✅ Aims are assessable; success criteria measurable
- ✅ Each method names its assumptions + how they're checked
- ✅ At least one positive control with an expected result
- ✅ Data schemas explicit; reference data versioned
- ✅ Consistent terminology

## Style guidelines

- **Tone:** precise, mechanistic, object-first (state the data object, then the method over it).
- **Format:** markdown throughout; tables for schemas and datasets.
- **Checkboxes:** ✅ in-scope / met, ❌ out-of-scope.
- **Specificity:** concrete cohorts, versions, columns, and expected control values over abstractions.
- **Length:** comprehensive but scannable.

## Output confirmation

After writing the ADD:
1. Confirm the file path.
2. Summarize the question, aims, and headline method.
3. List assumptions made where information was missing.
4. Suggest next steps — typically: run `plan-analysis` on the first stage; set up the canonical params.

## Notes

- If a critical decision is missing (which cohort, which estimand, which control), ask before generating.
- For method-heavy analyses, emphasize Methods, Assumptions, and Validation.
- For data-heavy projects, emphasize Data schemas, QC, and Reproducibility.
- This skill contains the full ADD template — no external references needed.
