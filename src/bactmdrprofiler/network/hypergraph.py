"""
Hypergraph representation of multi-drug resistance patterns.

Novel contributions
-------------------
1.  **MDR patterns as hyperedges** — each unique co-resistance pattern is a
    hyperedge connecting drug classes, preserving multi-way structure that
    pairwise graphs lose (Battiston et al. 2020).
2.  **Hypergraph centrality** — node participation weighted by hyperedge
    size and support.
3.  **Interaction information** — detects *k*-wise synergistic vs redundant
    resistance using inclusion-exclusion on joint entropy (McGill 1954).
    Positive II = synergistic (pairwise misses it); negative = redundant.
"""
from __future__ import annotations
from collections import Counter
from itertools import combinations
from typing import List
import numpy as np
import pandas as pd

# ── hyperedge extraction ─────────────────────────────────────────────────

def extract_hyperedges(
    class_df: pd.DataFrame,
    min_size: int = 2,
    max_size: int = 6,
    min_support: float = 0.05,
) -> pd.DataFrame:
    """Each MDR pattern (set of resistant classes) → one hyperedge."""
    n = len(class_df)
    pats = []
    for _, row in class_df.iterrows():
        fs = frozenset(c for c in class_df.columns if row[c] == 1)
        if min_size <= len(fs) <= max_size:
            pats.append(fs)
    counts = Counter(pats)
    rows = []
    for pat, cnt in counts.most_common():
        sup = cnt / n
        if sup >= min_support:
            rows.append({"Hyperedge": ", ".join(sorted(pat)),
                         "_set": pat, "Size": len(pat),
                         "Count": cnt, "Support": round(sup, 4)})
    return pd.DataFrame(rows)


# ── hypergraph centrality ────────────────────────────────────────────────

def hypergraph_centrality(
    hedges_df: pd.DataFrame,
    all_classes: List[str],
) -> pd.DataFrame:
    """Degree / weighted-degree / size-weighted centrality per drug class."""
    if not all_classes:
        return pd.DataFrame(columns=[
            "Class", "Degree", "Weighted_Degree", "Size_Weighted",
            "Degree_Norm", "Weighted_Degree_Norm", "Size_Weighted_Norm",
        ])
    deg = {c: 0 for c in all_classes}
    wd  = {c: 0.0 for c in all_classes}
    swd = {c: 0.0 for c in all_classes}
    for _, r in hedges_df.iterrows():
        for m in r["_set"]:
            if m in deg:
                deg[m] += 1
                wd[m] += r["Support"]
                swd[m] += r["Support"] * r["Size"]
    rows = [{"Class": c, "Degree": deg[c],
             "Weighted_Degree": round(wd[c], 4),
             "Size_Weighted": round(swd[c], 4)} for c in all_classes]
    df = pd.DataFrame(rows).sort_values("Size_Weighted", ascending=False)
    for col in ["Degree", "Weighted_Degree", "Size_Weighted"]:
        mx = df[col].max()
        if mx > 0: df[f"{col}_Norm"] = (df[col] / mx).round(4)
    return df


# ── interaction information ──────────────────────────────────────────────

def _joint_H(cols: list[np.ndarray]) -> float:
    """Joint Shannon entropy (bits)."""
    combo = np.column_stack(cols)
    _, cnts = np.unique(combo, axis=0, return_counts=True)
    p = cnts / cnts.sum()
    return float(-np.sum(p * np.log2(p + 1e-15)))


def interaction_information(
    df: pd.DataFrame,
    features: List[str],
    max_order: int = 3,
    min_support: float = 0.05,
) -> pd.DataFrame:
    """
    Interaction information (co-information) for feature subsets.

    II > 0  → synergistic (jointly informative beyond pairwise)
    II < 0  → redundant (pairwise captures relationship)
    """
    rows = []
    for order in range(3, max_order + 1):
        for sub in combinations(features, order):
            joint_prev = (df[list(sub)].sum(axis=1) == order).mean()
            if joint_prev < min_support: continue
            cols = [df[f].values for f in sub]
            n = len(sub)
            ii = 0.0
            for sz in range(1, n + 1):
                sign = (-1) ** (n - sz)
                for idx in combinations(range(n), sz):
                    ii += sign * _joint_H([cols[i] for i in idx])
            typ = "synergistic" if ii > 0.01 else ("redundant" if ii < -0.01 else "independent")
            rows.append({"Features": ", ".join(sub), "Order": order,
                         "II": round(ii, 4), "Type": typ,
                         "Joint_Prevalence": round(joint_prev, 4)})
    result = pd.DataFrame(rows)
    if not result.empty:
        result = result.sort_values("II", key=abs, ascending=False)
    return result
