"""YAML-driven configuration — user-supplied antibiotic ontology, any organism.

v1.1 changes
----------
- Adds schema_version and config validation with explicit reporting of unknown keys.
- Unknown keys are reported (WARN by default) and ignored; enable strict mode to fail.
- Adds `min_classes_matched` hard gate to prevent silent MDR degeneracy when ontology
  does not match any input columns.

Strict mode triggers
--------------------
- YAML key: config_strict: true
- Env var: SSUIS_CONFIG_STRICT=1
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List
import dataclasses as _dc
import logging
import os
from typing import get_args, get_origin

import yaml

logger = logging.getLogger(__name__)

SCHEMA_VERSION = "1.1"
SUPPORTED_SCHEMA_VERSIONS = {"1.0", "1.1"}


def _is_dataclass_type(tp: Any) -> bool:
    try:
        return hasattr(tp, "__dataclass_fields__")
    except Exception:
        return False


def _dataclass_allowed_fields(dc_cls: Any) -> set[str]:
    return {f.name for f in _dc.fields(dc_cls)}


def _validate_and_filter(dc_cls: Any, raw: Any, prefix: str, unknown_paths: List[str]) -> Any:
    if raw is None:
        raw = {}
    if isinstance(raw, dict):
        allowed = _dataclass_allowed_fields(dc_cls)
        filtered: dict = {}
        for k, v in raw.items():
            if k not in allowed:
                unknown_paths.append(f"{prefix}.{k}" if prefix else str(k))
                continue
            f = next((ff for ff in _dc.fields(dc_cls) if ff.name == k), None)
            if f is None:
                filtered[k] = v
                continue
            tp = f.type
            if _is_dataclass_type(tp) and isinstance(v, dict):
                filtered[k] = _validate_and_filter(tp, v, f"{prefix}.{k}" if prefix else k, unknown_paths)
            else:
                origin = get_origin(tp)
                args = get_args(tp)
                if origin in (list, List) and args and _is_dataclass_type(args[0]) and isinstance(v, list):
                    out_list = []
                    for i, item in enumerate(v):
                        if isinstance(item, dict):
                            out_list.append(_validate_and_filter(args[0], item, f"{prefix}.{k}[{i}]" if prefix else f"{k}[{i}]", unknown_paths))
                        else:
                            out_list.append(item)
                    filtered[k] = out_list
                else:
                    filtered[k] = v
        return filtered
    return raw


@dataclass
class OntologySpec:
    """Maps antibiotic class names → column names in the input CSV."""
    classes: Dict[str, List[str]] = field(default_factory=dict)
    mdr_threshold: int = 3  # ≥ N classes resistant = MDR

    @property
    def all_drugs(self) -> List[str]:
        return [d for ds in self.classes.values() for d in ds]


@dataclass
class CausalSpec:
    enabled: bool = True
    algorithm: str = "pc"          # 'pc' or 'fci'
    alpha: float = 0.05
    max_cond_set: int = 3


@dataclass
class HypergraphSpec:
    enabled: bool = True
    min_pattern: int = 2
    max_pattern: int = 6
    min_support: float = 0.05
    max_interaction_order: int = 3


@dataclass
class ProbabilisticSpec:
    enabled: bool = True
    model: str = "bernoulli_by_class"    # for missing class calls
    n_mc: int = 5000


@dataclass
class GeneLayerSpec:
    format: str = "auto"                # auto | wide | long
    feature_column: str | None = None
    value_column: str | None = None


@dataclass
class CompressionSpec:
    enabled: bool = True
    redundancy_jaccard: float = 0.9     # drop patterns with >= this overlap and <= support
    keep_top: int = 100


@dataclass
class NetworkSpec:
    alpha: float = 0.05
    fdr_method: str = "fdr_bh"
    motif_sizes: List[int] = field(default_factory=lambda: [3, 4])


@dataclass
class Config:
    # schema + validation
    schema_version: str = SCHEMA_VERSION
    config_strict: bool = False

    input_csv: str = ""
    gene_csv: str = ""                # optional: AMR gene presence/absence
    output_dir: str = "mdr_results"
    id_column: str = "Strain_ID"

    ontology: OntologySpec = field(default_factory=OntologySpec)
    causal: CausalSpec = field(default_factory=CausalSpec)
    hypergraph: HypergraphSpec = field(default_factory=HypergraphSpec)
    network: NetworkSpec = field(default_factory=NetworkSpec)
    gene_layer: GeneLayerSpec = field(default_factory=GeneLayerSpec)
    probabilistic: ProbabilisticSpec = field(default_factory=ProbabilisticSpec)
    compression: CompressionSpec = field(default_factory=CompressionSpec)

    # hard gates
    min_classes_matched: int = 1

    n_bootstrap: int = 5000
    confidence: float = 0.95
    n_jobs: int = -1
    seed: int = 42
    test_costs: Dict[str, float] = field(default_factory=dict)

    @classmethod
    def from_yaml(cls, path: str | Path) -> "Config":
        raw = yaml.safe_load(Path(path).read_text()) or {}

        schema_in = str(raw.get("schema_version", "1.0"))
        strict = bool(raw.get("config_strict", False)) or (os.environ.get("SSUIS_CONFIG_STRICT", "0") == "1")
        if schema_in not in SUPPORTED_SCHEMA_VERSIONS:
            msg = f"Unsupported schema_version={schema_in!r}. Supported: {sorted(SUPPORTED_SCHEMA_VERSIONS)}"
            if strict:
                raise ValueError(msg)
            logger.warning(msg)

        unknown_paths: List[str] = []
        filtered = _validate_and_filter(cls, raw, prefix="", unknown_paths=unknown_paths)

        ont = OntologySpec(**(filtered.pop("ontology", {}) or {}))
        cau = CausalSpec(**(filtered.pop("causal", {}) or {}))
        hyp = HypergraphSpec(**(filtered.pop("hypergraph", {}) or {}))
        net = NetworkSpec(**(filtered.pop("network", {}) or {}))
        gl = GeneLayerSpec(**(filtered.pop("gene_layer", {}) or {}))
        pb = ProbabilisticSpec(**(filtered.pop("probabilistic", {}) or {}))
        cp = CompressionSpec(**(filtered.pop("compression", {}) or {}))

        cfg = cls(ontology=ont, causal=cau, hypergraph=hyp, network=net, gene_layer=gl, probabilistic=pb, compression=cp, **filtered)

        cfg._config_validation = {
            "schema_version_in": schema_in,
            "schema_version_effective": cfg.schema_version,
            "supported_schema_versions": sorted(SUPPORTED_SCHEMA_VERSIONS),
            "unknown_keys": sorted(set(unknown_paths)),
            "strict": strict,
            "status": "PASS" if (schema_in in SUPPORTED_SCHEMA_VERSIONS and not unknown_paths) else ("WARN" if not strict else "FAIL"),
        }

        if unknown_paths:
            msg = f"Unknown config keys ignored ({len(set(unknown_paths))}): {sorted(set(unknown_paths))[:20]}" + (" ..." if len(set(unknown_paths)) > 20 else "")
            if strict:
                raise ValueError(msg)
            logger.warning(msg)

        return cfg

    def to_yaml(self, path: str | Path) -> None:
        import dataclasses as dc

        def ser(o):
            if dc.is_dataclass(o):
                return {k: ser(v) for k, v in dc.asdict(o).items()}
            if isinstance(o, list):
                return [ser(v) for v in o]
            return o

        Path(path).write_text(yaml.dump(ser(self), default_flow_style=False, sort_keys=False))
