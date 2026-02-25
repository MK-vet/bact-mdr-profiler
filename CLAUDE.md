# bact-mdr-profiler — Claude Code Context

## Project Identity

- **Tool**: bact-mdr-profiler v1.0.0
- **Package**: `bactmdrprofiler` (in `src/bactmdrprofiler/`)
- **Author**: Maciej Kochanowski, National Veterinary Research Institute (PIWet-PIB), Pulawy, Poland
- **Domain**: Antimicrobial resistance (AMR) epidemiology, *Streptococcus suis* and related bacteria
- **Purpose**: Causal co-resistance discovery with hypergraph MDR modelling — PC algorithm, Shapley values, EVPI, hypergraph centrality for MDR/XDR/PDR classification
- **Publication**: Individual SoftwareX paper — do NOT mix scope with sister tools

## Sister Tools (Anti-Salami Boundaries)

| Tool | Scope | THIS tool does NOT do this |
|------|-------|---------------------------|
| bact-assoc-net | MI, phi, OR, CMH, PID, TE, simplicial topology | NO pairwise association networks |
| bact-phylo-trait | Pagel lambda, Fritz-Purvis D, Fitch ASR, phylo-logistic | NO phylogenetic comparative methods |
| **bact-mdr-profiler** (THIS) | PC-algorithm causal, hypergraph MDR, Shapley, EVPI | — |
| bact-trait-cluster | K-Modes consensus, CKA fusion, NVI, TDA, SHAP | NO clustering, NO consensus fusion |

Reference: `docs/anti_salami_checklist.json`

## Architecture

```
src/bactmdrprofiler/
  config.py              — YAML config v1.1, dataclass schema, strict mode
  pipeline.py            — Main orchestrator: load → MDR classify → causal → hypergraph → manifest
  cli.py                 — CLI entry: --config, --self-check, --benchmark, --reliability-*
  io/loader.py           — CSV loading, wide/long auto-detect, NA-preserving Int8
  mdr/
    classification.py    — MDR/XDR/PDR classification by drug class ontology
    decision_v3.py       — Decision framework with Shapley and EVPI
  causal/
    discovery.py         — PC algorithm for CPDAG structure learning
  network/
    hypergraph.py        — Hypergraph construction and centrality measures
  reliability/core.py    — Preflight QC, consistency checks, marimo warnings
  selfcheck.py           — Internal validation
  benchmark.py           — Synthetic benchmark
  gui/                   — Marimo interactive app
  dashboards/            — Plotly dashboard
```

**Key patterns**: config-driven pipeline, NA-preserving (never NA→0), provenance via `run_manifest.json` + `config_used.yaml`, deterministic seeds, ontology-driven MDR classification.

## Development Quickstart

```bash
pip install -e ".[dev]"
ruff check src/
ruff format --check src/
python -m pytest tests/ -v --tb=short
python -m bactmdrprofiler.cli --self-check
python -m bactmdrprofiler.cli --benchmark
python -m bactmdrprofiler.cli --config examples/config.yaml
```

## CRITICAL: Scientific Correctness Rules

**These rules are ABSOLUTE. Violating any of them produces scientifically invalid results.**

1. **NEVER fabricate** p-values, test statistics, confidence intervals, or effect sizes. All must come from actual computation on actual data.
2. **NEVER generate fake citations**. Valid references for this tool:
   - Spirtes, Glymour & Scheines (2000) — PC algorithm
   - Magiorakos et al. (2012) — MDR/XDR/PDR definitions
   - Shapley (1953) — Shapley value theory
   - Raiffa & Schlaifer (1961) — EVPI decision theory
3. **NEVER convert NA to 0**. NA means "not tested"; 0 means "tested and susceptible". MDR classification with NA→0 underestimates resistance.
4. **NEVER apply continuous-variable methods** to binary (0/1) resistance data.
5. **NEVER skip multiple testing correction** when testing multiple causal edges. Apply FDR (Benjamini-Hochberg).
6. **NEVER assume independence** of samples with shared phylogenetic origin without explicit justification.
7. **NEVER add functionality belonging to a sister tool** (see anti-salami table above).

## Statistical Constraints — THIS TOOL ONLY

- **MDR classification**: MDR requires resistance to >= 3 antimicrobial CLASSES (not individual drugs); ontology MUST map drug columns to classes correctly
- **PC algorithm**: Output is a CPDAG (not DAG); some edges may be undirected; CMH conditioning for binary data
- **Shapley values**: Measure contribution of resistance PATTERNS (not ML feature importance); must sum to prediction value (efficiency axiom)
- **EVPI**: Requires a probabilistic MDR classification model; measures value of perfect information for decision-making
- **Hypergraph**: Hyperedges connect >= 2 drugs in co-resistance patterns; centrality measures (betweenness, closeness) on hypergraph structure

## Data Format Constraints

- Input: CSV with binary resistance features (0 = susceptible, 1 = resistant, NA = not tested)
- Drug-to-class ontology: YAML mapping of drug column names to antimicrobial classes
- Internal dtype: nullable Int8 (`pd.Int8Dtype()`)
- Formats: wide (samples x drugs) or long (ID, drug, value) — auto-detected
- NA handling: preserve throughout; MDR classification handles NA explicitly
- Output: CSV with classifications + causal edges + hypergraph metrics, provenance JSON

## Testing Conventions

- **Real data**: Run on `examples/` CSV files with example configs
- **Synthetic data**: Generated with known resistance patterns (`ground_truth.json`)
- **Signal recovery**: Verify MDR/XDR/PDR classification matches known patterns; PC recovers known DAG skeleton
- **Multi-seed**: Run with multiple seeds for permutation-based tests
- **Edge cases**: All-NA drug columns, single-class resistance, ontology mismatches
- **Reproducibility**: Same seed + config → identical output (bit-for-bit)

## Publication Context (SoftwareX)

- CITATION.cff, LICENSE, README.md exist — **DO NOT MODIFY** these files
- Scope boundary defined in `docs/anti_salami_checklist.json` — **DO NOT MODIFY**
- Reproducibility artifacts: `run_manifest.json` (input SHA-256), `config_used.yaml` (config snapshot)

## Code Style

- Python 3.10+ (type hints with `X | Y` union syntax)
- Linter: ruff (check + format)
- Logging: `logging.getLogger(__name__)`, never print()
- Tests: pytest in `tests/`, parametrized where possible
