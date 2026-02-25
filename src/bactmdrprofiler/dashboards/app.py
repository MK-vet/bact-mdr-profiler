import marimo
from pathlib import Path


app = marimo.App(width="full")


@app.cell
def _():
    import os
    import sys
    import json
    import time
    import subprocess
    from pathlib import Path

    import marimo as mo
    import pandas as pd
    import yaml

    return Path, json, mo, os, pd, subprocess, sys, time, yaml


@app.cell
def _(mo):
    mo.md(
        """
        # bact-mdr-profiler — interactive dashboard (marimo)

        This dashboard focuses on MDR/XDR/PDR spectra, Bayesian prevalence, and gene–phenotype networks.
        Runs are started in a separate process to keep the UI responsive.

        Key safety gates implemented by the tool:
        - Ontology mapping coverage (`min_classes_matched`): prevents degenerate MDR outputs when MIC columns don't match.
        """
    )
    return


@app.cell
def _(mo):
    data_root = mo.ui.text(value=str(Path.cwd()), label="Data root")
    out_root = mo.ui.text(value=str(Path.cwd() / "bactmdrprofiler_out"), label="Output directory")

    mic = mo.ui.text(value="MIC.csv", label="MIC CSV")
    amr = mo.ui.text(value="AMR_genes.csv", label="AMR genes CSV")

    min_classes = mo.ui.number(value=1, label="min_classes_matched", step=1)
    mdr_thr = mo.ui.number(value=3, label="MDR threshold (#classes)", step=1)
    alpha = mo.ui.number(value=0.05, label="FDR alpha", step=0.01)
    config_strict = mo.ui.switch(value=False, label="Config strict mode")
    schema_version = mo.ui.text(value="1.1", label="schema_version")

    mo.md("## Config builder")
    mo.vstack(
        [
            data_root,
            out_root,
            mo.hstack([mic, amr]),
            mo.hstack([min_classes, mdr_thr, alpha]),
            mo.hstack([config_strict, schema_version]),
        ]
    )
    return alpha, amr, config_strict, data_root, mdr_thr, mic, min_classes, out_root, schema_version


@app.cell
def _(Path, mo, yaml, alpha, amr, config_strict, data_root, mdr_thr, mic, min_classes, out_root, schema_version):
    def _abs(p: str) -> str:
        p = p.strip()
        pp = Path(p)
        if pp.is_absolute():
            return str(pp)
        return str(Path(data_root.value) / pp)

    cfg = {
        "schema_version": schema_version.value,
        "config_strict": bool(config_strict.value),
        "output_dir": str(Path(out_root.value)),
        "mic_path": _abs(mic.value),
        "amr_path": _abs(amr.value),
        "min_classes_matched": int(min_classes.value),
        "mdr_threshold": int(mdr_thr.value),
        "fdr_alpha": float(alpha.value),
        "reliability": {"enabled": True, "fail_fast": False},
    }
    editor = mo.ui.text_area(value=yaml.safe_dump(cfg, sort_keys=False), label="Config YAML (editable)", rows=18)
    mo.md("## Config editor")
    editor
    return editor


@app.cell
def _(Path, json, mo, subprocess, sys, time, editor, out_root):
    mo.md("## Run controls")
    start = mo.ui.button(label="Start analysis", kind="success")
    refresh = mo.ui.button(label="Refresh status")
    mo.hstack([start, refresh])

    out_dir = Path(out_root.value)
    out_dir.mkdir(parents=True, exist_ok=True)
    cfg_path = out_dir / "config_used.yaml"
    log_path = out_dir / "run.log"
    state_path = out_dir / "run_state.json"

    if start.value:
        cfg_path.write_text(editor.value, encoding="utf-8")
        cmd = [sys.executable, "-m", "bactmdrprofiler.cli", "--config", str(cfg_path)]
        with open(log_path, "ab") as f:
            proc = subprocess.Popen(cmd, stdout=f, stderr=subprocess.STDOUT, cwd=str(out_dir))
        state_path.write_text(
            json.dumps({"pid": proc.pid, "cmd": cmd, "started_at": time.time(), "cwd": str(out_dir)}, indent=2),
            encoding="utf-8",
        )
    return cfg_path, log_path, refresh, state_path


@app.cell
def _(mo, log_path, refresh):
    if not log_path.exists():
        mo.md("No log yet.")
    else:
        with open(log_path, "rb") as f:
            f.seek(0, 2)
            size = f.tell()
            f.seek(max(0, size - 20000), 0)
            tail = f.read().decode("utf-8", errors="replace")
        mo.md("### Log (tail)")
        mo.code(tail)
    return


@app.cell
def _(Path, mo, pd, out_root):
    mo.md("## Results explorer")
    out_dir = Path(out_root.value)
    pattern = mo.ui.text(value="*.csv", label="File glob")
    limit = mo.ui.number(value=200, label="Preview rows", step=50)
    label_map = mo.ui.text_area(value="{}", label="Label map (JSON) — applied to displayed tables", rows=4)
    mo.vstack([mo.hstack([pattern, limit]), label_map])
    files = sorted([p for p in out_dir.rglob(pattern.value) if p.is_file()]) if out_dir.exists() else []
    if not files:
        mo.md("No files yet.")
        return
    display = files[:200]
    selector = mo.ui.dropdown(
        options=[str(p.relative_to(out_dir)) for p in display],
        value=str(display[0].relative_to(out_dir)),
        label=f"Select a file ({len(display)}/{len(files)})",
    )
    file_path = out_dir / selector.value
    selector
    if file_path.suffix.lower() == ".csv":
        df = pd.read_csv(file_path, nrows=int(limit.value))
        try:
            mapping = __import__("json").loads(label_map.value or "{}")
            if isinstance(mapping, dict) and mapping:
                df = df.rename(columns={k: v for k, v in mapping.items() if k in df.columns})
                for c in df.columns:
                    if df[c].dtype == object:
                        df[c] = df[c].replace(mapping)
        except Exception:
            pass
        mo.ui.table(df, pagination=True, label=f"Preview: {selector.value}")
    else:
        mo.md(f"Selected: `{selector.value}`")
    return


@app.cell
def _(Path, mo, out_root, pd):
    mo.md("## Quick plots")
    out_dir = Path(out_root.value)
    prev = out_dir / "prevalence.csv"
    if not prev.exists():
        mo.md("`prevalence.csv` not found yet.")
        return
    try:
        df = pd.read_csv(prev)
    except Exception as e:
        mo.md(f"Failed to read prevalence.csv: `{e}`")
        return

    # Expected columns: Class, posterior_mean, ci_low, ci_high (names may vary)
    cls = "Class" if "Class" in df.columns else df.columns[0]
    mean = "posterior_mean" if "posterior_mean" in df.columns else ("mean" if "mean" in df.columns else None)
    lo = "ci_low" if "ci_low" in df.columns else None
    hi = "ci_high" if "ci_high" in df.columns else None
    if mean is None:
        mo.md("Could not infer posterior mean column.")
        return

    df2 = df.sort_values(mean, ascending=False).head(20)

    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(8, 4))
    ax.bar(df2[cls].astype(str), df2[mean].astype(float))
    ax.set_ylabel("Prevalence (posterior mean)")
    ax.set_title("Top classes by prevalence")
    ax.tick_params(axis="x", labelrotation=45)
    mo.mpl.interactive(fig)
    return


@app.cell
def _(mo):
    mo.md("## Output reference")
    artifacts = [
        {"file": "ontology_mapping_report.csv/json", "meaning": "How MDR classes mapped to MIC columns.", "interpretation": "If mapping is empty/low, fix ontology or MIC column names."},
        {"file": "mdr_spectrum.csv", "meaning": "Per-isolate #classes and MDR/XDR/PDR labels.", "interpretation": "Verify that class counts match expectations from MIC calls."},
        {"file": "mdr_probability.csv", "meaning": "Probabilistic MDR outputs (if enabled).", "interpretation": "Use as uncertainty-aware summary, not ground truth."},
        {"file": "prevalence.csv", "meaning": "Bayesian prevalence per class.", "interpretation": "Use posterior mean + CrI; compare classes across cohorts."},
        {"file": "network_edges.csv", "meaning": "Edges among genes/phenotypes.", "interpretation": "Association graph; tune thresholds and compare to phylogeny."},
        {"file": "motifs.csv", "meaning": "Repeated local patterns.", "interpretation": "Highlights recurring co-resistance/marker motifs."},
        {"file": "hyperedges.csv", "meaning": "Frequent multi-feature co-occurrence patterns.", "interpretation": "Use for MDR signatures and cluster-defining sets."},
        {"file": "hypergraph_centrality.csv", "meaning": "Central features in the hypergraph.", "interpretation": "Candidate hubs/diagnostic markers."},
        {"file": "interaction_information.csv", "meaning": "Synergy/redundancy of feature triplets.", "interpretation": "Synergy suggests combinatorial signals; control for prevalence."},
        {"file": "config_validation.json", "meaning": "Config schema + unknown-key report.", "interpretation": "Fix config if strict mode fails."},
        {"file": "run_manifest.json", "meaning": "Input hashes + environment.", "interpretation": "Reproducibility."},
    ]
    mo.ui.table(artifacts, pagination=True, label="Artifacts")

    mo.md(
        """
        ## Statistical assumptions (practical)
        - Prevalence uses Jeffreys prior Beta(0.5,0.5): stabilises intervals for small n.
        - MDR labels depend on ontology mapping; inspect `ontology_mapping_report` first.
        - Gene–phenotype edges are association-based and can be confounded by clonality.
        """
    )
    return


if __name__ == "__main__":
    app.run()
