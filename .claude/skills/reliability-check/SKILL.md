---
name: reliability-check
description: Run reliability preflight QC and quality gate check
---

```bash
python -m bactmdrprofiler.cli --config examples/config.yaml --reliability-only
```

## What to check:
1. Drug column-to-class ontology mapping (all columns mapped?)
2. NA proportion per drug (flag >50% missing)
3. Degenerate columns (zero variance, all-NA)
4. Minimum 3 antimicrobial classes represented (for MDR definition)
5. Quality gate status (PASS/WARN/FAIL)

Report quality gate result and warnings with remediation suggestions.
