# bact-mdr-profiler

[![CI](https://github.com/MK-vet/bact-mdr-profiler/actions/workflows/ci.yaml/badge.svg)](https://github.com/MK-vet/bact-mdr-profiler/actions/workflows/ci.yaml)
[![PyPI](https://img.shields.io/pypi/v/bact-mdr-profiler.svg)](https://pypi.org/project/bact-mdr-profiler/)
[![Python](https://img.shields.io/pypi/pyversions/bact-mdr-profiler.svg)](https://pypi.org/project/bact-mdr-profiler/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)


Causal co-resistance discovery with hypergraph multi-drug resistance modelling.

## Novel Contributions

1. **PC algorithm for causal co-resistance** — constraint-based skeleton discovery with Cochran–Mantel–Haenszel stratified conditional independence tests. Distinguishes direct co-resistance from confounded associations.
2. **Hypergraph MDR patterns** — each co-resistance profile is a hyperedge, preserving multi-way structure lost in pairwise graphs (Battiston et al. 2020)
3. **Interaction information** — *k*-wise synergy/redundancy detection via inclusion-exclusion on joint entropy (McGill 1954)
4. **Bayesian prevalence** — Jeffrey's prior Beta(0.5, 0.5) for small-sample credible intervals

## Installation

```bash
pip install bact-mdr-profiler
pip install bact-mdr-profiler[gui]         # interactive dashboard (marimo)
```

## Interactive dashboard (marimo)

```bash
bact-mdr-profiler-dashboard
```

Edit mode:

```bash
bact-mdr-profiler-dashboard --edit
```

## Quick Start

```bash
bact-mdr-profiler config.yaml -v
```

### Python API

```python
from bactmdrprofiler.config import Config
from bactmdrprofiler.pipeline import Pipeline

cfg = Config.from_yaml("config.yaml")
results = Pipeline(cfg).run()
```

## Configuration

The YAML config contains a user-supplied **antibiotic ontology** mapping drug classes to column names — **zero hardcoding of species or antibiotics**.

```yaml
ontology:
  classes:
    Penicillins: [PEN, AMP, AMX]
    Tetracyclines: [TET, DOX]
    Macrolides: [ERY]
  mdr_threshold: 3
```

See `examples/config.yaml` for a complete template.

## Outputs

| File | Description |
|------|-------------|
| `mdr_spectrum.csv` | MDR/XDR/PDR per isolate |
| `prevalence.csv` | Bayesian prevalence per class |
| `network_edges.csv` | Typed network edges (pheno–pheno, gene–gene, pheno–gene) |
| `co_resistance.graphml` | Network for Cytoscape/Gephi |
| `causal_edges.csv` | PC algorithm edge decisions |
| `causal_graph.graphml` | Partially oriented causal graph |
| `hyperedges.csv` | MDR patterns with support |
| `hypergraph_centrality.csv` | Node importance in hypergraph |
| `interaction_information.csv` | Synergy/redundancy per feature subset |
| `motifs.csv` | Network motif census |

## Testing

```bash
pip install -e ".[dev]"
pytest tests/ -v
```

## License

MIT
