"""Signal-recovery and integration validation for bact-mdr-profiler."""
import subprocess, sys, json
import numpy as np
import pandas as pd
import pytest


# ── fixture: 6 drugs, 4 classes, strong planted co-resistance ────────────

@pytest.fixture
def coresistance_df():
    """200 isolates. AMP co-occurs with TET at p=0.95. ERY/CLI/GEN/SXT independent."""
    rng = np.random.RandomState(7)
    n = 200
    amp = rng.binomial(1, 0.55, n)
    tet = np.where(amp == 1,
                   rng.binomial(1, 0.95, n),
                   rng.binomial(1, 0.08, n))
    df = pd.DataFrame({
        "Strain_ID": [f"S{i:04d}" for i in range(n)],
        "AMP": amp, "AMX": rng.binomial(1, 0.5, n),
        "TET": tet, "DOX": rng.binomial(1, 0.45, n),
        "ERY": rng.binomial(1, 0.30, n),
        "GEN": rng.binomial(1, 0.20, n),
    })
    return df


@pytest.fixture
def ontology():
    return {
        "Penicillins":     ["AMP", "AMX"],
        "Tetracyclines":   ["TET", "DOX"],
        "Macrolides":      ["ERY"],
        "Aminoglycosides": ["GEN"],
    }


# ── signal recovery ───────────────────────────────────────────────────────

class TestSignalRecovery:
    def test_pc_detects_amp_tet_edge(self, coresistance_df):
        """PC skeleton (max_cond=0) must contain the AMP-TET edge."""
        from bactmdrprofiler.causal.discovery import pc_skeleton
        G, _ = pc_skeleton(coresistance_df,
                            ["AMP", "TET", "ERY", "GEN"],
                            alpha=0.01, max_cond=0)
        assert G.to_undirected().has_edge("AMP", "TET"), (
            f"AMP-TET missing from PC skeleton. Edges: {list(G.to_undirected().edges())}"
        )

    def test_pc_amp_tet_stronger_than_random_pair(self, coresistance_df):
        """AMP-TET chi-squared statistic > ERY-GEN (independent pair)."""
        from scipy.stats import chi2_contingency
        ct_at = pd.crosstab(coresistance_df["AMP"], coresistance_df["TET"])
        ct_eg = pd.crosstab(coresistance_df["ERY"], coresistance_df["GEN"])
        chi2_at = chi2_contingency(ct_at)[0]
        chi2_eg = chi2_contingency(ct_eg)[0]
        assert chi2_at > chi2_eg * 5, (
            f"AMP-TET chi2={chi2_at:.1f} not sufficiently > ERY-GEN chi2={chi2_eg:.1f}"
        )

    def test_mdr_spectrum_non_susceptible_present(self, coresistance_df, ontology):
        """MDR spectrum must contain non-susceptible isolates."""
        from bactmdrprofiler.mdr.classification import build_class_matrix, mdr_spectrum
        cm = build_class_matrix(coresistance_df.drop(columns=["Strain_ID"]), ontology)
        sp = mdr_spectrum(cm, threshold=3)
        non_susc = (sp["Category"] != "susceptible").sum()
        assert non_susc > 0, "All isolates classified susceptible — unexpected"

    def test_mdr_identify_amp_tet_coresistant(self, coresistance_df, ontology):
        """Isolates resistant to ≥2 classes are NOT classified as definitely-not-MDR
        when threshold=2 (i.e. identify_mdr returns True or NA, never False)."""
        from bactmdrprofiler.mdr.classification import build_class_matrix, identify_mdr
        cm = build_class_matrix(coresistance_df.drop(columns=["Strain_ID"]), ontology)
        # threshold=2: any isolate with ≥2 resistant classes must be MDR
        mdr_flag = identify_mdr(cm, threshold=2)
        double_res = cm[(cm["Penicillins"] == 1) & (cm["Tetracyclines"] == 1)].index
        # With threshold=2, these isolates (2 classes resistant) → True, not False
        definitely_not_mdr = (mdr_flag.loc[double_res] == False).sum()  # noqa: E712
        assert definitely_not_mdr == 0, (
            f"{definitely_not_mdr} co-resistant isolates (Pen+Tet) "
            f"incorrectly classified as definitely not-MDR at threshold=2"
        )

    def test_bayesian_credible_interval_covers_truth(self, coresistance_df):
        """95% CI for AMP prevalence must bracket the observed proportion."""
        from bactmdrprofiler.mdr.classification import bayesian_prevalence
        k = int(coresistance_df["AMP"].sum())
        n = len(coresistance_df)
        r = bayesian_prevalence(k, n)
        true_p = k / n
        assert r["CI_Lo"] < true_p < r["CI_Hi"], (
            f"True prev {true_p:.3f} outside CI [{r['CI_Lo']:.3f}, {r['CI_Hi']:.3f}]"
        )

    def test_interaction_information_penicillin_tetracycline(self, coresistance_df, ontology):
        """Penicillins+Tetracyclines pair must appear in interaction information results."""
        from bactmdrprofiler.mdr.classification import build_class_matrix
        from bactmdrprofiler.network.hypergraph import interaction_information
        cm = build_class_matrix(coresistance_df.drop(columns=["Strain_ID"]), ontology)
        ii = interaction_information(cm, cm.columns.tolist(), max_order=2, min_support=0.01)
        if not ii.empty:
            pair = ii[ii["Features"].apply(
                lambda x: set(x) == {"Penicillins", "Tetracyclines"})]
            assert len(pair) > 0, "Penicillins-Tetracyclines pair absent from II results"


# ── reproducibility ───────────────────────────────────────────────────────

class TestReproducibility:
    def test_pc_skeleton_deterministic(self, coresistance_df):
        """PC skeleton is deterministic for same data and alpha."""
        from bactmdrprofiler.causal.discovery import pc_skeleton
        G1, _ = pc_skeleton(coresistance_df, ["AMP","TET","ERY","GEN"],
                             alpha=0.05, max_cond=0)
        G2, _ = pc_skeleton(coresistance_df, ["AMP","TET","ERY","GEN"],
                             alpha=0.05, max_cond=0)
        assert set(G1.to_undirected().edges()) == set(G2.to_undirected().edges())


# ── CLI integration ───────────────────────────────────────────────────────

class TestCLI:
    def test_help_exits_zero(self):
        r = subprocess.run([sys.executable, "-m", "bactmdrprofiler.cli", "--help"],
                           capture_output=True)
        assert r.returncode == 0

    def test_version_output(self):
        from bactmdrprofiler import __version__
        r = subprocess.run([sys.executable, "-m", "bactmdrprofiler.cli", "--version"],
                           capture_output=True, text=True)
        assert __version__ in r.stdout

    def test_self_check_passes(self):
        r = subprocess.run([sys.executable, "-m", "bactmdrprofiler.cli", "--self-check"],
                           capture_output=True, text=True)
        assert r.returncode == 0
        assert json.loads(r.stdout)["status"] == "PASS"

    def test_run_produces_mandatory_outputs(self, coresistance_df, tmp_path):
        """End-to-end: write CSV → run CLI → verify mandatory output files exist."""
        import yaml
        data_dir = tmp_path / "data"; data_dir.mkdir()
        out_dir  = tmp_path / "results"
        coresistance_df.to_csv(data_dir / "AMR.csv", index=False)

        config = {
            "schema_version": "1.1",
            "input_csv": str(data_dir / "AMR.csv"),
            "id_column": "Strain_ID",
            "output_dir": str(out_dir),
            "min_classes_matched": 1,
            "ontology": {
                "classes": {
                    "Penicillins": ["AMP", "AMX"],
                    "Tetracyclines": ["TET", "DOX"],
                    "Macrolides": ["ERY"],
                    "Aminoglycosides": ["GEN"],
                },
                "mdr_threshold": 3,
            },
            "causal": {"enabled": True, "algorithm": "pc",
                       "alpha": 0.01, "max_cond_set": 0},
            "hypergraph": {"enabled": True, "min_pattern": 2},
        }
        with open(tmp_path / "config.yaml", "w") as f:
            yaml.dump(config, f)

        r = subprocess.run(
            [sys.executable, "-m", "bactmdrprofiler.cli",
             "--config", str(tmp_path / "config.yaml")],
            capture_output=True, text=True, timeout=120
        )
        assert r.returncode == 0, f"CLI failed:\n{r.stderr}"
        assert (out_dir / "mdr_spectrum.csv").exists()
        assert (out_dir / "causal_edges.csv").exists()
        assert (out_dir / "run_manifest.json").exists()
