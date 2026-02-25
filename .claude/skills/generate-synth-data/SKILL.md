---
name: generate-synth-data
description: Generate synthetic resistance data with known MDR/causal ground truth
---

Generate synthetic binary resistance matrices with known properties.

## Scenarios

1. **MDR/XDR/PDR classification** (seed=42, n=100, 15 drugs in 6 classes):
   - Known resistance patterns → MDR/XDR/PDR labels deterministic
   - Ground truth: exact classification for each sample

2. **PC algorithm DAG recovery** (seed=43, n=200, 5 drugs):
   - Known DAG (A→B→C, A→D, E independent)
   - Generate data from known conditional probabilities
   - Ground truth: adjacency matrix of true DAG

3. **V-structure** (seed=44, n=300):
   - A→C←B with known conditional probabilities
   - Ground truth: correct V-structure orientation

4. **Shapley additivity** (seed=45, n=200):
   - Known additive model for MDR status
   - Ground truth: Shapley values must sum to prediction

5. **Hypergraph centrality** (seed=46, n=150, 10 drugs):
   - Known co-resistance patterns forming specific hyperedges
   - Ground truth: expected centrality ranking

6. **Ontology edge cases** (seed=47):
   - Drug names mapping to 1, 2, 3+ classes
   - Boundary cases: exactly 3 classes (MDR), all classes (PDR)
   - Ground truth: exact MDR/XDR/PDR boundary

7. **Missing data** (seed=48, n=200):
   - Random vs systematic missingness in resistance profiles
   - Ground truth: MDR status under complete data

## Implementation
Python script, deterministic seeds, output to `synth_data/`.
