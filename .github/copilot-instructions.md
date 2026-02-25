# bact-mdr-profiler — AI Coding Assistant Instructions

## Context
This is **bact-mdr-profiler**, a Python tool for causal co-resistance discovery with hypergraph MDR modelling. It implements the PC algorithm for CPDAG structure learning, MDR/XDR/PDR classification per Magiorakos et al. (2012), Shapley values for pattern contribution, and EVPI decision theory. Published as an individual SoftwareX paper.

## CRITICAL: Scientific Correctness Rules (NEVER violate)
1. NEVER fabricate p-values, test statistics, CI, or effect sizes
2. NEVER generate fake citations (valid: Spirtes et al. 2000, Magiorakos et al. 2012, Shapley 1953)
3. NEVER convert NA to 0 (NA = "not tested" — critical for MDR classification)
4. NEVER apply continuous-variable methods to binary resistance data
5. NEVER skip FDR correction on PC algorithm edge tests
6. NEVER assume sample independence without phylogenetic justification
7. NEVER add functionality from sister tools (see scope below)

## Scope Boundary (Anti-Salami)
- THIS tool: PC-algorithm causal, hypergraph MDR, Shapley, EVPI
- NOT this tool: pairwise association networks, phylogenetic methods, clustering/consensus

## Code Conventions
- Python 3.10+, ruff, pytest, `logging.getLogger(__name__)`
- Data: binary 0/1/NA CSV, nullable Int8, drug-to-class ontology YAML

## Dev Commands
```bash
pip install -e ".[dev]"
ruff check src/ && ruff format --check src/
python -m pytest tests/ -v --tb=short
python -m bactmdrprofiler.cli --self-check
```

## Read-Only Files (DO NOT MODIFY)
- CITATION.cff, LICENSE, docs/anti_salami_checklist.json, examples/*
