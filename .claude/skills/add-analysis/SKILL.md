---
name: add-analysis
description: Guided workflow to add a new analysis step to the MDR pipeline
---

## Pre-flight

1. **Scope check**: Does this analysis involve causal discovery, MDR classification, or hypergraph modelling?
   - Pairwise associations/networks → bact-assoc-net
   - Phylogenetic methods → bact-phylo-trait
   - Clustering/consensus → bact-trait-cluster

2. **Method validation**: Appropriate for binary resistance data?
3. **Ontology check**: Does it need the drug-to-class ontology?

## Implementation checklist

4. Add parameters to `src/bactmdrprofiler/config.py`
5. Add step to `src/bactmdrprofiler/pipeline.py`
6. Add tests with known resistance patterns
7. Run `/validate-config`, `/run-tests`, `/selfcheck`, `/reliability-check`
8. Run `/cross-tool-check`

## Post-implementation

9. Update README.md if user-facing
10. Verify provenance
