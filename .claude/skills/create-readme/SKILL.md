---
name: create-readme
description: Create or update the top-level project README (subproject index) or a subproject/stage README, using canonical scientific-pipeline templates
argument-hint: [top-level | <subproject-dir>]
---

# create-readme

Write and maintain two kinds of README for a scientific pipeline repo:

1. **Top-level project README** — the *index*. One paragraph of purpose, the inputs/outputs the
   project consumes/produces, the canonical runs, and a numbered **workflow order** that lists each
   subproject as a linked list item. This is the entry point a returning author reads first.
2. **Subproject (stage) README** — one per analysis stage, following a fixed canonical template so
   every stage is documented the same way.

The two are linked: every subproject **anchors itself** into the top-level README as one list item
under "Workflow order", and every subproject README links back up (`Related → Upstream/Downstream`).

---

## Operating rules (read before writing)

- **Inspect first, write second.** Read the folder you are documenting: `ls`/glob the files; read
  the orchestrator entry (`main.nf` / `Snakefile` / driver script), the params/config files, and
  the `bin/` scripts. Grep for cross-references in and out. **Derive inputs, outputs, and workflow
  steps from code — never invent them.**
- **Don't claim non-obvious facts without a source.** `Status:` and `Last touched:` come from prose
  already present or from `git log -1 --format=%ai -- <dir>`. Otherwise write `unknown`.
- **Read parents and siblings** so the new README threads into the workflow chain. The top-level
  README is the authoritative index; resolve every `Upstream:`/`Downstream:` link against a file
  that exists (verify with a read/glob) before writing it.
- **Match the house voice:** terse, mechanistic, present-tense. No emojis, no sales language.
- **One README per invocation.** Don't fan out to siblings unless asked.

---

## Which template?

- Argument `top-level` (or you're at repo root and asked for "the project README") → **Template A**.
- Argument is a subproject directory (or you're inside a stage folder) → **Template B**.
- No argument → infer from cwd: a folder that *contains* stage subfolders → A; a leaf stage → B.

---

## Template A — top-level project README (the index)

Render this structure. It mirrors the project's own `CLAUDE.md` but is reader-facing and
navigational; do not duplicate CLAUDE.md's deep convention text — link to it.

````markdown
# <project-name>

Top-level index for the project's analytical pipeline. Each subfolder is a self-contained stage
with its own README. Read them in the order below to catch up on the work.

## Purpose

<2–4 sentences: the scientific question and the shape of the work (greenfield vs brownfield).>

## Inputs and outputs

- **Inputs:** <raw data types + formats; reference panels; external cohorts>
- **Intermediate contracts:** <the interchange files between stages and their schema, e.g.
  "cluster manifest = (phenotype_id, cluster_label, vcf_path)">
- **Terminal outputs:** <final artifacts: latent GWAS, gene results, figures, exports>

## Canonical datasets and runs

| Dataset | Canonical run / solution | Notes |
|---|---|---|
| <name> | <the params/config that is "the answer"> | <why; key gotcha> |

## Workflow order

<Numbered list — ONE list item per subproject. Each item: linked folder, a one-line role, and
optionally its status. THIS is where a new subproject anchors itself.>

1. [setup/](setup/README.md) — <one-line role>. *one-time.*
2. [preprocessing/](preprocessing/README.md) — <one-line role>.
3. [pipeline/](pipeline/README.md) — <one-line role>.
4. ...

## Data flow at a glance

```
<ascii diagram of the stage DAG: forks, joins, terminal stages>
```

- <2–3 bullets on what flows between stages and what must never be crossed>

## Conventions

<Short pointers only — link to CLAUDE.md for the authoritative convention text.>
- **Compute:** <where stages submit from / run>
- **Storage:** <bucket/scratch root + subtree layout>
- **Schemas:** <the one or two most load-bearing column contracts, or "see CLAUDE.md">
````

### Anchoring a new subproject into the index

When a new stage is created, add (don't replace) one list item under **Workflow order** at the
correct position in the DAG, and, if the stage changes the data flow, update the diagram. Keep the
list item to a single line: `N. [<dir>/](<dir>/README.md) — <role>.` Verify the link resolves.

---

## Template B — subproject / stage README

Render every header below (permissive on *additions* — extra domain subsections are fine when the
stage has substantive content — but never omit a header). If a section is empty, say so explicitly
(`Downstream: none — terminal`) rather than dropping it.

````markdown
# <folder-name> — <one-line what-and-why>

**Stage:** <position in the pipeline, e.g. "step 6 of 8 — full-resolution gPCA + MAGMA gene tests">
**Status:** <complete | active | scaffolded | deferred | one-time>; <which datasets/runs are done>
**Last touched:** <YYYY-MM from `git log -1 -- <folder>`>

## Orientation

**Start reading:** <entry-point file> → <next> → <next>
**Common tasks here:** <2–4 things a returning author/Claude will do — "add a params variant",
"debug a failed batch task", "rerun for a new cohort">
**Don't touch without thinking:** <files/configs with non-obvious dependents>
**Non-obvious state:** <where data actually lives if not under output/; external state that matters>

## Purpose

<2–4 sentences: what this stage does and why it exists. If it diverges from a sibling stage, name
the sibling and the reason.>

## Inputs

| Item | Source | Format | Notes |
|---|---|---|---|
| <name> | <upstream stage or storage path> | <TSV/VCF/h5ad/…> | <schema, constraint, gotcha> |

## Outputs

| Item | Location | Format | Consumed by |
|---|---|---|---|
| <name> | `<outdir>/<path>` or `gs://…` | <format> | <downstream stage> |

## Workflow

1. `<PROCESS_OR_SCRIPT>` (`<file>`) — <one line>.
2. ...

<If a Nextflow/Snakemake pipeline, mirror the process/rule names exactly. If standalone scripts,
number them in execution order.>

## How to run

```bash
# Canonical, copy-pasteable invocation
<command>
```

<State the submission context: head/login node vs local; profile name.>

## Key decisions / gotchas

- **<decision>** — <one-line rationale; surface anything that would surprise a re-reader>.
- ...

## Results

<Until the stage produces results: "Pending — no canonical run yet.">
<Once it does: 2–4 lines of the HEADLINE result, then link to the full report:>
See [`REPORT.md`](REPORT.md) for the full results, findings, and exploratory interpretation.
<Generate/refresh REPORT.md with the `report-findings` skill.>

## Related

- Upstream: <link(s) to upstream stage README(s)>
- Downstream: <link(s) to downstream stage README(s)> — or `none — terminal`
- Sibling: <if applicable>
````

---

## Workflow when invoked

1. **Resolve mode** (A vs B) from the argument or cwd.
2. **Inspect** the target folder (or, for A, the set of stage subfolders): read orchestrator entry,
   params/configs, `bin/` scripts; grep for `gs://`/`s3://`/`/scratch` paths and cross-references.
3. **Read parent + siblings** so links thread correctly. For A, enumerate the stage subfolders in
   DAG order. For B, find this stage's place in the top-level index.
4. **Derive `Last touched:`** from `git log -1 --format=%ai -- <folder>`.
5. **Write/refresh:**
   - **New README:** render every section with real file names, real paths, real process names. No placeholders.
   - **Existing README:** if `## Orientation` (B) or `## Results` (B) is missing, backfill it in place;
     reconcile structural drift against the template but **preserve substantive existing prose** — don't
     rewrite accurate content just to match phrasing.
   - For B, also **anchor the stage** into Template A's "Workflow order" (add the one-line list item if absent).
6. **Cross-reference integrity:** verify every link resolves (read/glob) before writing it.
7. **Report:** path(s) written; one-line summary of change (new vs backfill vs reconcile); any open
   question you couldn't resolve from code (unclear status, no entry point, missing canonical run).

## Anti-patterns to avoid

- Inventing `Inputs`/`Outputs` rows from extension guessing — read the script.
- Linking to `Upstream:`/`Downstream:` paths you didn't verify.
- Writing `Status: complete` with no evidence from prose or git log.
- Adding emojis or marketing tone. Match the terse, present-tense house voice.
- Bypassing the template "because this stage is different" — keep the headers; mark empties explicitly.
