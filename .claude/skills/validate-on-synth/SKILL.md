---
name: validate-on-synth
description: Validate signal recovery on synthetic MDR/causal data
---

Run pipeline on synthetic data and verify signal recovery.

## Validation per scenario

1. **MDR/XDR/PDR**: 100% classification accuracy on clean data
2. **PC algorithm**: CPDAG skeleton matches true DAG adjacency (SHD = 0 or low)
3. **V-structure**: Correct orientation A→C←B detected
4. **Shapley**: |sum(values) - prediction| < 1e-6 for each sample
5. **Hypergraph**: Centrality ranking matches ground truth (Spearman rho > 0.8)
6. **Ontology**: Correct MDR/XDR/PDR boundary classification
7. **Missing data**: Classification degrades gracefully; no false PDR from missingness

## Multi-seed
For PC algorithm: seeds {42,43,44,45,46}, compute SHD mean/variance.

Report: scenario | metric | expected | observed | PASS/FAIL
