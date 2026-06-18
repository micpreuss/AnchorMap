---
name: execute
description: Execute an implementation/analysis plan
argument-hint: [path-to-plan]
---

# Execute: Implement from Plan

## Plan to Execute

Read plan file: `$ARGUMENTS`

## Execution Instructions

### 1. Read and Understand

- Read the ENTIRE plan carefully
- Understand all tasks and their dependencies
- Note the validation commands to run
- Review the testing / validation strategy

### 2. Execute Tasks in Order

For EACH task in "Step by Step Tasks":

#### a. Navigate to the task
- Identify the file and action required
- Read existing related files if modifying (and the real data headers before writing any parser)

#### b. Implement the task
- Follow the detailed specifications exactly
- Maintain consistency with existing code patterns and the project's conventions (schemas, units, sign preservation)
- Document non-obvious choices; log where the project logs
- Prefer config/params over hard-coded flags where the project does

#### c. Verify as you go
- After each change, check syntax (and that the orchestrator graph still parses, e.g. `nextflow ... -preview` / `snakemake -n`)
- Ensure imports/dependencies resolve
- Confirm any new data schema matches the documented column contract

### 3. Implement the Validation Strategy

After completing implementation tasks:

- Create the tests / smoke configs specified in the plan
- Implement the positive/negative controls the plan names
- Cover the edge cases the plan lists

### 4. Run Validation Commands

Execute ALL validation commands from the plan in order:

```bash
# Run each command exactly as specified in the plan
```

If any command fails:
- Fix the issue
- Re-run the command
- Continue only when it passes

### 5. Final Verification

Before completing:

- ✅ All tasks from the plan completed
- ✅ Tests / smoke runs created and passing
- ✅ All validation commands pass
- ✅ Code follows project conventions and schema contracts
- ✅ Documentation (stage README / CLAUDE.md) updated as needed

## Output Report

### Completed Tasks
- Tasks completed; files created (paths); files modified (paths)

### Validation Added
- Tests / smoke configs / controls implemented and their results

### Validation Results
```bash
# Output from each validation command
```

### Ready for Commit
- Confirm changes complete and validations pass; ready for the `commit` skill

## Notes

- If you hit issues not addressed in the plan, document them
- If you must deviate from the plan, explain why
- If validation fails, fix until it passes — don't skip validation steps
