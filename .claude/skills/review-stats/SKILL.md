---
name: review-stats
description: Review statistical methods for correctness — causal discovery and MDR audit
---

Review causal discovery and MDR classification implementations.

## Checklist

1. **MDR/XDR/PDR classification**:
   - MDR = resistant to >= 3 antimicrobial CLASSES (not drugs)?
   - XDR = resistant to all but <= 2 classes?
   - PDR = resistant to ALL classes?
   - Ontology correctly maps drug columns to classes?
   - NA handled explicitly (not counted as susceptible)?

2. **PC algorithm**:
   - Output is CPDAG (some edges may be undirected)?
   - CMH conditioning used for binary data?
   - Skeleton discovery correct (adjacency test)?
   - V-structure orientation correct?
   - FDR correction on edge tests?

3. **Shapley values**:
   - Pattern contribution (not ML feature importance)?
   - Efficiency axiom: sum of Shapley values = prediction?
   - Symmetry: equal contributors get equal value?

4. **EVPI**:
   - Requires probabilistic model?
   - Correctly computes expected value of perfect information?

5. **Hypergraph**:
   - Hyperedges connect >= 2 drugs?
   - Centrality measures (betweenness, closeness) correct on hypergraph?

## Key files
- `src/bactmdrprofiler/mdr/classification.py`, `src/bactmdrprofiler/mdr/decision_v3.py`
- `src/bactmdrprofiler/causal/discovery.py`
- `src/bactmdrprofiler/network/hypergraph.py`

Report: method | check | status (PASS/WARN/FAIL) | evidence (file:line)
