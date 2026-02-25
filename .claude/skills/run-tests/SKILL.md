---
name: run-tests
description: Run the test suite for bact-mdr-profiler
---

Run the full test suite with verbose output and short tracebacks.

```bash
python -m pytest tests/ -v --tb=short
```

Report: total passed/failed/skipped, any failures with file:line, confirm "All tests PASS" if clean.
