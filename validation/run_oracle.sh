#!/usr/bin/env bash
# run_oracle.sh — end-to-end cross-language parity check for the AnchorMap R engine (Phase 1).
# Runs the R engine on the anthro + disease configs and compares each against the committed Python
# reference output in the sibling cluster_anchoring repo. Exits non-zero on any deterministic mismatch.
#
# (Does NOT re-run the Python reference, to avoid mutating the sibling repo's committed output. To
#  regenerate it: cd <parent>/scripts/cluster_anchoring && python3 anchor_categories.py --config configs/<name>.yaml)
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PAR="/Users/mpreuss/Library/Mobile Documents/com~apple~CloudDocs/Desktop/projects/Northwell/project/UKBB_CLUSTER_GWAS"
ORA="$PAR/scripts/cluster_anchoring/output"
cd "$HERE"

# The Carey/FinnGen configs are machine-specific (absolute input paths) and live in the gitignored
# local/configs/. The engine is the installed `anchormap` package, driven via its CLI wrapper.
status=0
for run in carey_rint15_anthro carey_rint15; do
  echo "================ $run ================"
  Rscript inst/scripts/anchor_map.R --config "local/configs/$run.yaml" --out-dir "results/$run" >/dev/null
  Rscript validation/compare_oracle.R \
    --r-out  "results/$run/category_anchor_scores.tsv" \
    --oracle "$ORA/$run/category_anchor_scores.tsv" || status=1
  echo
done

echo "================ unit tests ================"
Rscript -e 'testthat::test_local()' | tail -3

exit $status
