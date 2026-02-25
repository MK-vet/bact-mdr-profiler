---
name: validate-on-real
description: Validate pipeline on real data from examples/
---

```bash
python -m bactmdrprofiler.cli --config examples/config.yaml
```

## Check invariants:
- No crash (exit code 0)
- Output files exist
- `run_manifest.json` and `config_used.yaml` generated
- No NA→0 conversion
- MDR/XDR/PDR classifications valid (mutually exclusive, ontology-consistent)
- Causal edges in CPDAG format (some may be undirected)
- Shapley values sum correctly per sample
- All p-values in [0, 1]

Report: structured checklist PASS/FAIL.
