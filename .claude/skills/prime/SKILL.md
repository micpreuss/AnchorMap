---
name: prime
description: Prime agent with project context (scientific data/pipeline project)
---

# prime

Load project context quickly via `CLAUDE.md`, the workflow-index README, and minimal git state.

## Objective

Fast orientation on what the project studies, how data flows, and where things are. Leans on
`CLAUDE.md` as the authoritative source; supplements with the top-level README and a structure
snapshot.

## Process

### 1. Read CLAUDE.md (the source of truth)

Project purpose, project type, data flow, conventions, canonical datasets/runs, how to run, and
validation patterns all live here.

```bash
cat CLAUDE.md
```

### 2. Read the workflow-index README

The top-level README that lists the stages/subprojects and the order to read them. It's the entry
point to the analytical pipeline.

```bash
# whichever exists — the index is usually at repo root or one level down
cat README.md 2>/dev/null || cat scripts/README.md 2>/dev/null
```

### 3. Minimal structure check

```bash
ls -la | grep -E '^d'          # top-level dirs
ls scripts/ 2>/dev/null        # stages, if the project uses a scripts/ tree
```

### 4. Current state (only if needed)

```bash
git status
git branch -v
git log -1 --format='%ai  %s'
```

## Output report

Keep it scannable — bullets, not prose. Link to memory files / stage READMEs for depth.

### Project purpose
- The scientific goal; greenfield vs brownfield; key cohorts/datasets.

### Project type (from CLAUDE.md)
- Orchestration / compute backend / domain / data locality / reproducibility (the five-axis summary).

### Pipeline stages
- The high-level stage order (e.g. setup → preprocessing → … → terminal exports) and which run is canonical.

### Tech stack & methods
- Languages, orchestration, compute backend, containers/envs; the named statistical/analytical methods.

### Key conventions
- The load-bearing data schemas (interchange-file column contracts), units/sign rules, config-over-CLI,
  and any standout gotchas.

### Current focus (if git checked)
- Active branch; uncommitted changes, if present.
