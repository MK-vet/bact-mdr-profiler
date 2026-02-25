"""Tests for bact-mdr-profiler causal, hypergraph, and classification modules."""

import numpy as np
import pandas as pd
import pytest


@pytest.fixture
def pheno_df():
    """60 isolates, 6 drugs, 3 antibiotic classes with co-resistance signal."""
    rng = np.random.RandomState(42)
    n = 60
    df = pd.DataFrame(
        {
            "Strain_ID": [f"S{i:03d}" for i in range(n)],
            "AMP": rng.binomial(1, 0.6, n),
            "AMX": rng.binomial(1, 0.55, n),
            "TET": rng.binomial(1, 0.5, n),
            "DOX": rng.binomial(1, 0.45, n),
            "ERY": rng.binomial(1, 0.3, n),
            "AZI": rng.binomial(1, 0.25, n),
        }
    )
    # Plant co-resistance: AMP→TET co-occur
    df.loc[df["AMP"] == 1, "TET"] = rng.binomial(1, 0.85, (df["AMP"] == 1).sum())
    return df


@pytest.fixture
def ontology():
    return {
        "Penicillins": ["AMP", "AMX"],
        "Tetracyclines": ["TET", "DOX"],
        "Macrolides": ["ERY", "AZI"],
    }


class TestClassification:
    def test_build_class_matrix(self, pheno_df, ontology):
        from bactmdrprofiler.mdr.classification import build_class_matrix

        cm = build_class_matrix(pheno_df.drop(columns=["Strain_ID"]), ontology)
        assert cm.shape == (60, 3)
        assert set(cm.columns) == {"Penicillins", "Tetracyclines", "Macrolides"}
        assert all(pd.api.types.is_integer_dtype(cm[c].dtype) for c in cm.columns)

    def test_mdr_spectrum(self, pheno_df, ontology):
        from bactmdrprofiler.mdr.classification import build_class_matrix, mdr_spectrum

        cm = build_class_matrix(pheno_df.drop(columns=["Strain_ID"]), ontology)
        sp = mdr_spectrum(cm, threshold=3)
        assert "Category" in sp.columns
        assert set(sp["Category"].unique()).issubset(
            {"susceptible", "MDR", "XDR", "PDR"}
        )

    def test_bayesian_prevalence(self):
        from bactmdrprofiler.mdr.classification import bayesian_prevalence

        r = bayesian_prevalence(10, 100)
        assert 0 < r["Posterior_Mean"] < 1
        assert r["CI_Lo"] < r["CI_Hi"]
        # Edge case: 0 positives
        r0 = bayesian_prevalence(0, 50)
        assert r0["Posterior_Mean"] > 0  # Jeffrey's prior never gives 0


class TestCausal:
    def test_pc_skeleton_returns_graph(self, pheno_df):
        from bactmdrprofiler.causal.discovery import pc_skeleton

        feats = ["AMP", "AMX", "TET", "DOX", "ERY", "AZI"]
        G, det = pc_skeleton(pheno_df, feats, alpha=0.05, max_cond=1)
        assert G.number_of_nodes() == 6
        assert "Node1" in det.columns

    def test_motif_census(self, pheno_df):
        from bactmdrprofiler.causal.discovery import pc_skeleton, motif_census

        feats = ["AMP", "TET", "ERY"]
        G, _ = pc_skeleton(pheno_df, feats, alpha=0.2, max_cond=0)
        ug = G.to_undirected()
        if ug.number_of_edges() >= 3:
            mdf = motif_census(ug, [3])
            assert "Motif" in mdf.columns

    def test_hybrid_network(self, pheno_df):
        from bactmdrprofiler.causal.discovery import build_hybrid_network

        pdf = pheno_df[["AMP", "TET", "ERY"]]
        G, edf = build_hybrid_network(pdf, None, alpha=0.2)
        assert "Feature1" in edf.columns
        assert "Type" in edf.columns


class TestHypergraph:
    def test_extract_hyperedges(self, pheno_df, ontology):
        from bactmdrprofiler.mdr.classification import build_class_matrix
        from bactmdrprofiler.network.hypergraph import extract_hyperedges

        cm = build_class_matrix(pheno_df.drop(columns=["Strain_ID"]), ontology)
        he = extract_hyperedges(cm, min_support=0.01)
        assert "Size" in he.columns
        assert (he["Size"] >= 2).all()

    def test_interaction_information_range(self, pheno_df, ontology):
        from bactmdrprofiler.mdr.classification import build_class_matrix
        from bactmdrprofiler.network.hypergraph import interaction_information

        cm = build_class_matrix(pheno_df.drop(columns=["Strain_ID"]), ontology)
        ii = interaction_information(
            cm, cm.columns.tolist(), max_order=3, min_support=0.01
        )
        if not ii.empty:
            assert "II" in ii.columns
            assert set(ii["Type"].unique()).issubset(
                {"synergistic", "redundant", "independent"}
            )


class TestConfig:
    def test_roundtrip(self, tmp_path):
        from bactmdrprofiler.config import Config, OntologySpec

        cfg = Config(
            input_csv="/tmp/x.csv", ontology=OntologySpec(classes={"Beta": ["AMP"]})
        )
        p = tmp_path / "c.yaml"
        cfg.to_yaml(p)
        cfg2 = Config.from_yaml(p)
        assert cfg2.ontology.classes == {"Beta": ["AMP"]}
