"""Pipeline: classify → network → causal → hypergraph → export.

Defensibility constraints:
- No implicit imputation: missing stays NA (not silently mapped to 0).
- Optional gene layer loaded as wide/long with explicit coverage reporting.
- A run_manifest.json captures input hashes for reproducibility.
"""

from __future__ import annotations
from . import __version__
import logging
from pathlib import Path
from typing import Dict
import json
import networkx as nx
import numpy as np
import pandas as pd
from .config import Config
from .mdr.classification import (
    build_class_matrix,
    identify_mdr,
    mdr_spectrum,
    prevalence_table,
    mdr_probability,
)
from .causal.discovery import pc_skeleton, motif_census, build_hybrid_network
from .network.hypergraph import (
    extract_hyperedges,
    hypergraph_centrality,
    interaction_information,
)
from .mdr.decision_v3 import (
    next_best_test_evpi,
    pattern_mdl_compression,
    shapley_pattern_contributions,
)
from .io import load_layer_csv, layer_coverage, qc_binary_features, sha256_file

logger = logging.getLogger(__name__)


def _save(df, p):
    num = df.select_dtypes(include=[np.number]).columns
    out = df.copy()
    out[num] = out[num].round(4)
    out.to_csv(p, index=False)
    logger.info("  → %s (%d rows)", p, len(out))


class Pipeline:
    def __init__(self, cfg: Config):
        self.cfg = cfg
        self.out = Path(cfg.output_dir)
        self.out.mkdir(parents=True, exist_ok=True)

    def _write_manifest(self, extra_inputs: list[str]) -> None:
        import platform
        import sys

        inp = [
            {
                "name": "phenotype",
                "path": self.cfg.input_csv,
                "sha256": sha256_file(self.cfg.input_csv),
            }
        ]
        if self.cfg.gene_csv:
            inp.append(
                {
                    "name": "gene_layer",
                    "path": self.cfg.gene_csv,
                    "sha256": sha256_file(self.cfg.gene_csv),
                }
            )
        for p in extra_inputs:
            inp.append({"name": "extra", "path": p, "sha256": sha256_file(p)})
        import dataclasses as dc

        try:
            config_snapshot = dc.asdict(self.cfg)
        except Exception:
            config_snapshot = None
        manifest = {
            "tool": "bact-mdr-profiler",
            "python": sys.version.split()[0],
            "platform": platform.platform(),
            "seed": self.cfg.seed,
            "tool_version": __version__,
            "inputs": inp,
        }
        manifest["config_snapshot"] = config_snapshot
        try:
            import numpy
            import pandas

            manifest["package_versions"] = {
                "numpy": numpy.__version__,
                "pandas": pandas.__version__,
            }
            try:
                import scipy

                manifest["package_versions"]["scipy"] = scipy.__version__
            except Exception:
                pass
            try:
                import networkx

                manifest["package_versions"]["networkx"] = networkx.__version__
            except Exception:
                pass
        except Exception:
            pass
        (self.out / "run_manifest.json").write_text(json.dumps(manifest, indent=2))

    def run(self) -> Dict[str, pd.DataFrame]:
        c = self.cfg
        res = {}

        # 1. Load phenotype
        logger.info("=== Loading data ===")
        raw = pd.read_csv(c.input_csv)
        raw[c.id_column] = raw[c.id_column].astype(str)
        raw = raw.set_index(c.id_column)

        # 2. Class-level resistance
        logger.info("=== MDR classification ===")
        class_df = build_class_matrix(raw, c.ontology.classes)

        # 2a. Ontology coverage report + hard gate against silent degeneracy
        declared = c.ontology.classes or {}
        raw_cols = set(raw.columns)
        matched = {
            cls: [d for d in (drugs or []) if d in raw_cols]
            for cls, drugs in declared.items()
        }
        unmatched = {
            cls: [d for d in (drugs or []) if d not in raw_cols]
            for cls, drugs in declared.items()
        }
        n_classes_declared = int(len(declared))
        n_classes_matched = int(sum(1 for cls, ds in matched.items() if len(ds) > 0))
        n_drugs_declared = int(sum(len(v) for v in declared.values()))
        n_drugs_matched = int(sum(len(v) for v in matched.values()))
        coverage_rep = {
            "n_classes_declared": n_classes_declared,
            "n_classes_matched": n_classes_matched,
            "n_drugs_declared": n_drugs_declared,
            "n_drugs_matched": n_drugs_matched,
            "matched_drug_columns_per_class": matched,
            "unmatched_drug_columns_per_class": unmatched,
        }
        (self.out / "ontology_mapping_report.json").write_text(
            json.dumps(coverage_rep, indent=2, ensure_ascii=False)
        )
        # also store a compact CSV for quick inspection
        map_rows = []
        for cls in sorted(declared):
            map_rows.append(
                {
                    "Class": cls,
                    "N_Matched": int(len(matched.get(cls, []))),
                    "Matched_Drugs": ";".join(matched.get(cls, [])),
                    "N_Unmatched": int(len(unmatched.get(cls, []))),
                    "Unmatched_Drugs": ";".join(unmatched.get(cls, [])),
                }
            )
        if map_rows:
            pd.DataFrame(map_rows).to_csv(
                self.out / "ontology_mapping_report.csv", index=False
            )
        if class_df.shape[1] == 0 or n_classes_matched < int(
            getattr(c, "min_classes_matched", 1)
        ):
            raise ValueError(
                f"Ontology mapping matched {n_classes_matched}/{n_classes_declared} classes and {n_drugs_matched}/{n_drugs_declared} drugs. "
                "This would make MDR classification degenerate. Fix ontology.classes to match MIC column names."
            )
        if n_classes_matched < int(c.ontology.mdr_threshold):
            logger.warning(
                "Only %d antibiotic classes matched (threshold=%d) — MDR calls may be impossible or rare.",
                n_classes_matched,
                int(c.ontology.mdr_threshold),
            )

        spectrum = mdr_spectrum(class_df, c.ontology.mdr_threshold)
        spectrum.insert(0, c.id_column, class_df.index.values)
        _save(spectrum, str(self.out / "mdr_spectrum.csv"))
        res["spectrum"] = spectrum

        prev = prevalence_table(class_df, c.confidence)
        _save(prev, str(self.out / "prevalence.csv"))
        res["prevalence"] = prev

        if c.probabilistic.enabled:
            pm = mdr_probability(class_df, c.ontology.mdr_threshold)
            _save(pm, str(self.out / "mdr_probability.csv"))
            res["p_mdr"] = pm

        # v3 decision-theoretic add-on: next best test (EVPI-like)
        try:
            evpi_costs = (
                getattr(c, "test_costs", None) if hasattr(c, "test_costs") else None
            )
            evpi = next_best_test_evpi(
                class_df, c.ontology.mdr_threshold, test_costs=evpi_costs
            )
            if not evpi.empty:
                _save(evpi, str(self.out / "next_best_test_evpi.csv"))
                res["evpi"] = evpi
        except Exception as e:
            logger.warning("EVPI block skipped: %s", e)

        mdr_mask = identify_mdr(class_df, c.ontology.mdr_threshold)
        n_def = int(mdr_mask.notna().sum())
        n_mdr = int(mdr_mask.eq(True).sum())
        n_unc = int(mdr_mask.isna().sum())
        rate = (n_mdr / n_def * 100) if n_def else float("nan")
        logger.info(
            "  MDR definite rate: %d/%d (%.1f%%), uncertain: %d",
            n_mdr,
            n_def,
            rate,
            n_unc,
        )

        # 3. Hybrid network
        logger.info("=== Hybrid co-resistance network ===")
        gene_df = None
        cov_rows = []
        qc_rows = []
        if c.gene_csv:
            lr = load_layer_csv(
                "AMR_genes",
                c.gene_csv,
                id_col=c.id_column,
                all_ids=class_df.index.tolist(),
                feature_col=c.gene_layer.feature_column,
                value_col=c.gene_layer.value_column,
            )
            gene_df = lr.data.select_dtypes(include=[np.number])
            # optional QC (drop invariant/ultra-rare/mostly missing features); no imputation
            gene_df, rep = qc_binary_features(
                gene_df,
                observed=lr.observed,
                min_prev=0.01,
                max_prev=0.99,
                max_missing_frac=0.5,
            )
            rep.insert(0, "Layer", "AMR_genes")
            qc_rows.append(rep)
            cov_rows.append(layer_coverage(lr))

        if cov_rows:
            cov_df = pd.DataFrame(cov_rows)
            _save(cov_df, str(self.out / "layer_coverage.csv"))
            res["coverage"] = cov_df
        if qc_rows:
            qc_df = pd.concat(qc_rows, ignore_index=True)
            _save(qc_df, str(self.out / "feature_qc.csv"))
            res["qc"] = qc_df

        G, edge_df = build_hybrid_network(
            class_df, gene_df, c.network.alpha, c.network.fdr_method
        )
        cols_to_save = [c for c in edge_df.columns if c != "_set"]
        _save(
            edge_df[cols_to_save] if "_set" in edge_df.columns else edge_df,
            str(self.out / "network_edges.csv"),
        )
        res["edges"] = edge_df
        nx.write_graphml(G, str(self.out / "co_resistance.graphml"))
        logger.info(
            "  Network: %d nodes, %d edges", G.number_of_nodes(), G.number_of_edges()
        )

        # Motifs
        if G.number_of_edges() > 2:
            mdf = motif_census(G, c.network.motif_sizes)
            _save(mdf, str(self.out / "motifs.csv"))
            res["motifs"] = mdf

        # 4. Causal discovery
        if c.causal.enabled:
            logger.info("=== Causal skeleton (PC) ===")
            feats = class_df.columns.tolist()
            if gene_df is not None:
                feats += gene_df.columns.tolist()
                combined = pd.concat(
                    [class_df.reset_index(drop=True), gene_df.reset_index(drop=True)],
                    axis=1,
                )
            else:
                combined = class_df.reset_index(drop=True)
            DG, det = pc_skeleton(
                combined, feats, c.causal.alpha, c.causal.max_cond_set
            )
            _save(det, str(self.out / "causal_edges.csv"))
            res["causal"] = det
            nx.write_graphml(DG, str(self.out / "causal_graph.graphml"))
            logger.info("  Causal edges: %d", DG.number_of_edges())

        # 5. Hypergraph MDR
        if c.hypergraph.enabled:
            logger.info("=== Hypergraph MDR ===")
            hedges = extract_hyperedges(
                class_df,
                c.hypergraph.min_pattern,
                c.hypergraph.max_pattern,
                c.hypergraph.min_support,
            )
            save_cols = [c2 for c2 in hedges.columns if c2 != "_set"]
            _save(hedges[save_cols], str(self.out / "hyperedges.csv"))
            res["hyperedges"] = hedges
            try:
                hmdl = pattern_mdl_compression(hedges)
                if not hmdl.empty:
                    _save(
                        hmdl[[c for c in hmdl.columns if c != "_set"]],
                        str(self.out / "hyperedges_mdl_compression.csv"),
                    )
                    res["hyperedges_mdl"] = hmdl
            except Exception as e:
                logger.warning("Hyperedge MDL skipped: %s", e)
            hc = hypergraph_centrality(hedges, class_df.columns.tolist())
            _save(hc, str(self.out / "hypergraph_centrality.csv"))
            res["hg_cent"] = hc

            # Interaction information
            ii = interaction_information(
                class_df,
                class_df.columns.tolist(),
                c.hypergraph.max_interaction_order,
                c.hypergraph.min_support,
            )
            _save(ii, str(self.out / "interaction_information.csv"))
            res["ii"] = ii

        # v3 explanatory output: class contributions to MDR probability
        try:
            shapley = shapley_pattern_contributions(class_df, c.ontology.mdr_threshold)
            if not shapley.empty:
                _save(shapley, str(self.out / "mdr_shapley_contributions.csv"))
                res["shapley"] = shapley
        except Exception as e:
            logger.warning("Shapley contribution block skipped: %s", e)

        self.cfg.to_yaml(str(self.out / "config_used.yaml"))
        self._write_manifest(extra_inputs=[str(self.out / "config_used.yaml")])
        logger.info("Done — %s", self.out)
        return res
