---
name: cross-tool-check
description: Verify anti-salami scope boundary and consistency with sister tools
---

Check that bact-mdr-profiler stays within its defined scope.

## Steps

1. **Read scope boundary**: Parse `docs/anti_salami_checklist.json`

2. **Grep for scope violations** in `src/bactmdrprofiler/`:
   - Association terms: `mutual_information`, `phi_coefficient`, `transfer_entropy`, `simplicial`
   - Phylo terms: `pagel`, `fritz_purvis`, `ancestral_state`, `newick`, `phylo_logistic`
   - Clustering terms: `k_modes`, `kmodes`, `consensus_matrix`, `NVI`, `CKA`, `kernel_fusion`

3. **Check shared interfaces** (if sister repos at `../bact-*`)

Report: PASS / WARN / FAIL with details.
