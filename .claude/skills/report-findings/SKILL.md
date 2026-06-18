---
name: report-findings
description: Write a short completion REPORT.md for a finished subproject — results, findings, and a clearly non-binding exploratory interpretation — and link it from the subproject README
argument-hint: [<subproject-dir>] [run-tag]
---

# report-findings

When a subproject (pipeline stage / analysis) reaches a reportable state, write a **short**
`REPORT.md` that a returning reader — or a collaborator who never saw the run — can read in two
minutes to know *what came out, what it shows, and what it might mean*.

The report is built around one hard distinction:

- **Findings** = statements the results **support**. Tied to a number or an output. Durable.
- **Exploratory interpretation** = **hypotheses** the findings *suggest*. Explicitly **non-binding** —
  "not consequences chiseled in stone", a quick way to build intuition, to be confirmed or dropped
  by later work. Hedged language only.

Keeping these two in separate, labeled sections is the point of the skill. Never let interpretation
masquerade as a finding.

---

## Where it goes (default: separate file + README pointer)

- **`REPORT.md` at the subproject root** — the canonical/latest run's report. This is the default.
- For a run-specific report, write `output/<run-tag>/REPORT.md` and have the root `REPORT.md`
  link to the per-run ones.
- **Update the subproject README's `## Results` section** with a 2–4 line headline + a link to
  `REPORT.md`. The README stays the scannable index; the report carries the detail. (Use the
  `create-readme` skill's Template B if the README has no `## Results` section yet.)

> Rationale for separate-file-plus-pointer: it matches the repo's existing pattern (a short Results
> section in the README, with the full narrative in `summary.md`/`validation_report.md` per run). The
> README stays short and durable; the longer, hedged interpretation lives where it won't bloat the
> index.

---

## Before writing — gather the evidence

Do not summarize from memory. Read the actual outputs.

- Read the subproject README (Purpose, Inputs, Outputs, How to run) for context.
- List and read the **output artifacts** (result tables, `summary.*`, validation reports, key figures).
  Pull the headline numbers directly from them.
- Recover **provenance**: which params/config produced this run, input dataset + version, container
  image/tag or env, the commit (`git rev-parse --short HEAD`), and the date (`git log -1` or run logs).
- Note any **positive/negative controls** the stage defined and whether they passed.

If a number isn't in an artifact you can read, don't put it in the report.

---

## REPORT.md template

````markdown
# <subproject> — Results report

**Run / dataset:** <run-tag> on <dataset + version>
**Status:** <final | interim>  ·  **Date:** <YYYY-MM-DD>  ·  **Commit:** <short-sha>
**Headline:** <one sentence: the single most important outcome.>

## Results

<The concrete outputs and numbers. A table is usually right. State the metric, the value, and where
it came from — not what it means yet.>

| Metric / output | Value | Source artifact |
|---|---|---|
| <e.g. clusters recovered> | <e.g. 15> | `output/<run>/…` |
| <e.g. positive control PP.H4 (HMGCR↔LDL)> | <e.g. 0.97 — PASS> | `output/<run>/…` |

<Controls: state explicitly which passed/failed.>

## Findings

<3–6 bullets. Each is a statement the results above SUPPORT, tied to a number/output. Factual,
mechanistic, present-tense. No speculation here.>

- <finding tied to a result>
- ...

## Exploratory interpretation  *(non-binding — hypotheses, not conclusions)*

> These readings are exploratory: a fast way to build intuition about what the findings *might*
> mean. They are **not** established results and should not be cited as conclusions. Treat each as a
> hypothesis to confirm, revise, or drop in later work.

- <hedged reading: "consistent with …", "may indicate …", "one explanation is …"> — *to confirm: <what would test it>.*
- ...

## Caveats / limitations

<What could undercut the findings: power, sample size, multiple testing, confounds, untested
assumptions, known data gotchas that apply to this run, scope the result does NOT cover.>

- ...

## Next steps  *(optional, non-committal)*

- <suggested follow-ups; clearly not promises>

## Provenance / reproducibility

- **Command:** `<exact invocation, e.g. nextflow run … -params-file … -profile …>`
- **Params/config:** `<path>`  ·  **Input:** `<dataset + version/path>`
- **Environment:** `<container image:tag or env lockfile>`  ·  **Commit:** `<short-sha>`
- **Outputs:** `<output dir / bucket path>`
````

### Subproject README `## Results` pointer (insert/refresh)

````markdown
## Results

<2–4 line headline of the canonical run's outcome, including whether controls passed.>
See [`REPORT.md`](REPORT.md) for full results, findings, and exploratory interpretation.
````

---

## Workflow when invoked

1. **Locate** the subproject (argument or cwd) and the run (argument or the canonical/latest run).
2. **Read** the README + output artifacts; extract headline numbers and control outcomes; recover provenance.
3. **Draft** `REPORT.md` from the template. Enforce the split:
   - Anything causal/speculative goes under **Exploratory interpretation** with hedged wording and a
     "to confirm" note — never under Findings.
   - Every Results-table number must trace to a named artifact.
4. **Write** `REPORT.md` (subproject root by default; `output/<run>/REPORT.md` for run-specific).
5. **Update the README:** insert/refresh the `## Results` headline + link. If the README lacks a
   `## Results` section, add it (per `create-readme` Template B) directly before `## Related`.
6. **Report back:** path(s) written; the one-sentence headline; and an explicit list of which
   statements you placed under *interpretation* vs *findings*, so the author can sanity-check the split.

## Anti-patterns to avoid

- **Interpretation creep** — writing "X causes Y" in Findings when the data only shows association.
  If you can't point at the number that proves it, it's interpretation.
- **Numbers without a source** — every value cites the artifact it came from.
- **Over-long reports** — this is a 1–2 minute read, not a manuscript. Push depth into the artifacts.
- **Silent overwrite** — if `REPORT.md` exists for a different run, don't clobber it; write the new
  run's report under `output/<run>/` and link both from the root.
- **Confident hedging** — interpretation bullets must actually read as tentative; pair each with what
  would test it.
