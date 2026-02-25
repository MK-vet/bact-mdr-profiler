---
name: lint
description: Run ruff linter and formatter check for bact-mdr-profiler
---

```bash
ruff check src/
ruff format --check src/
```

Report: lint errors/warnings, format check result, confirm "Lint PASS" if clean.
