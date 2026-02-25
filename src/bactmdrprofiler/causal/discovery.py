"""
Conditional independence structure learning for antimicrobial co-resistance.

Novel contributions
-------------------
1.  **PC algorithm on binary AMR data** — constraint-based conditional
    independence skeleton discovery using Cochran–Mantel–Haenszel stratified
    conditional independence tests.  Identifies which pairs of resistance
    classes remain associated after conditioning on others, consistent with
    shared mobile genetic elements or co-selective pressures.
    First application to AMR co-resistance data.

    IMPORTANT — causal interpretation caveat
    ----------------------------------------
    The PC algorithm recovers the Markov equivalence class of a causal DAG
    under the assumption of *causal sufficiency* (no latent common causes).
    In AMR data, mobile genetic elements, clonal population structure, and
    shared antimicrobial exposure are potential hidden confounders that violate
    causal sufficiency.  We therefore report only the **undirected conditional
    independence skeleton** and interpret edges as conditional dependence
    relationships, not directed causal effects.

2.  **Network motif census** — subgraph isomorphism class counts (triangles,
    stars, paths) revealing co-resistance architecture (Milo et al. 2002).

References
----------
Spirtes et al. (2000) Causation, Prediction, and Search. MIT Press.
Milo et al. (2002) Science 298:824.
"""
from __future__ import annotations
import logging
from itertools import combinations
from typing import Dict, List, Tuple
import networkx as nx
import numpy as np
import pandas as pd
from scipy.stats import chi2 as chi2_dist, chi2_contingency, fisher_exact

logger = logging.getLogger(__name__)


# ── conditional independence via CMH ─────────────────────────────────────

def _ci_test(
    df: pd.DataFrame, i: str, j: str, S: List[str], alpha: float,
) -> Tuple[bool, float]:
    """Test X_i ⊥ X_j | X_S using CMH for binary data."""
    if not S:
        tab = pd.crosstab(df[i], df[j])
        if tab.shape == (2, 2):
            exp = np.outer(tab.sum(axis=1), tab.sum(axis=0)) / tab.values.sum()
            if (exp < 5).any():
                _, p = fisher_exact(tab)
            else:
                _, p, _, _ = chi2_contingency(tab)
        else:
            try:
                _, p, _, _ = chi2_contingency(tab)
            except Exception:
                p = 1.0
        return p > alpha, float(p)

    # Stratified CMH
    key = df[S].apply(tuple, axis=1)
    num = den = 0.0
    for _, st in df.groupby(key):
        tab = pd.crosstab(st[i], st[j])
        if tab.shape != (2, 2): continue
        a, b, c, d = tab.values[0, 0], tab.values[0, 1], tab.values[1, 0], tab.values[1, 1]
        n = a + b + c + d
        if n < 2: continue
        num += a - (a + b) * (a + c) / n
        den += (a + b) * (c + d) * (a + c) * (b + d) / (n**2 * (n - 1))
    if den <= 0:
        return True, 1.0
    stat = num**2 / den
    p = 1 - chi2_dist.cdf(stat, 1)
    return p > alpha, float(p)


# ── PC algorithm ─────────────────────────────────────────────────────────

def pc_skeleton(
    df: pd.DataFrame,
    features: List[str],
    alpha: float = 0.05,
    max_cond: int = 3,
) -> Tuple[nx.DiGraph, pd.DataFrame]:
    """
    PC algorithm: conditional independence skeleton for binary AMR data.

    Recovers the undirected skeleton of conditional independence relationships
    using CMH-based tests. V-structures are oriented for completeness but the
    primary output for biological interpretation is the **undirected skeleton**
    (see module-level docstring for causal interpretation caveat).

    The algorithm assumes causal sufficiency (no latent common causes). In AMR
    data, latent confounders such as mobile genetic elements and clonal
    population structure may be present. Edges should therefore be interpreted
    as conditional dependence relationships rather than direct causal effects.

    Parameters
    ----------
    df : pd.DataFrame
        Binary AMR data (rows = isolates, columns = features).
    features : list of str
        Feature column names to include.
    alpha : float
        Significance level for conditional independence tests.
    max_cond : int
        Maximum conditioning set size.

    Returns
    -------
    G : nx.DiGraph
        Partially oriented conditional independence graph (skeleton).
    details : pd.DataFrame
        Edge-level test decisions with columns Node1, Node2, Removed,
        Cond_Set, P, Cond_Size.
    """
    G = nx.Graph()
    G.add_nodes_from(features)
    G.add_edges_from(combinations(features, 2))

    sep: Dict[Tuple[str, str], List[str]] = {}
    details = []

    # Phase 1: edge removal
    for sz in range(max_cond + 1):
        rm = []
        for u, v in list(G.edges()):
            nbrs = list(set(G.neighbors(u)) | set(G.neighbors(v)) - {u, v})
            if len(nbrs) < sz: continue
            for S in combinations(nbrs, sz):
                indep, p = _ci_test(df, u, v, list(S), alpha)
                if indep:
                    rm.append((u, v))
                    sep[(u, v)] = list(S); sep[(v, u)] = list(S)
                    details.append({"Node1": u, "Node2": v, "Removed": True,
                                    "Cond_Set": ",".join(S) if S else "∅",
                                    "P": round(p, 4), "Cond_Size": sz})
                    break
            else:
                if sz == 0:
                    _, p = _ci_test(df, u, v, [], alpha)
                    details.append({"Node1": u, "Node2": v, "Removed": False,
                                    "Cond_Set": "—", "P": round(p, 4), "Cond_Size": -1})
        for u, v in rm:
            if G.has_edge(u, v): G.remove_edge(u, v)

    # Phase 2: orient v-structures
    DG = nx.DiGraph()
    DG.add_nodes_from(G.nodes)
    for u, v in G.edges(): DG.add_edge(u, v); DG.add_edge(v, u)

    for k in G.nodes:
        for i, j in combinations(G.neighbors(k), 2):
            if not G.has_edge(i, j):
                if k not in sep.get((i, j), []):
                    if DG.has_edge(k, i): DG.remove_edge(k, i)
                    if DG.has_edge(k, j): DG.remove_edge(k, j)

    return DG, pd.DataFrame(details)


# ── motif census ─────────────────────────────────────────────────────────

def motif_census(G: nx.Graph, sizes: List[int] = None) -> pd.DataFrame:
    """Count connected subgraph isomorphism classes of given sizes."""
    from collections import Counter
    if sizes is None: sizes = [3, 4]
    counts: Dict[str, int] = Counter()
    nodes = list(G.nodes)
    for sz in sizes:
        if sz > len(nodes): continue
        for sub in combinations(nodes, sz):
            sg = G.subgraph(sub)
            if not nx.is_connected(sg): continue
            ne = sg.number_of_edges()
            degs = tuple(sorted((d for _, d in sg.degree()), reverse=True))
            counts[f"n{sz}_e{ne}_{degs}"] += 1
    return pd.DataFrame([{"Motif": m, "Count": c} for m, c in counts.most_common()])


# ── hybrid phenotype-genotype network ────────────────────────────────────

def build_hybrid_network(
    pheno_df: pd.DataFrame,
    gene_df: pd.DataFrame | None,
    alpha: float = 0.05,
    fdr_method: str = "fdr_bh",
) -> Tuple[nx.Graph, pd.DataFrame]:
    """
    Construct typed co-resistance network with edges:
    pheno↔pheno, gene↔gene, pheno↔gene (if gene_df provided).
    """
    from statsmodels.stats.multitest import multipletests

    edges = []
    def _test_pair(x, y, cat):
        tab = pd.crosstab(x, y)
        if tab.shape != (2, 2): return
        a, b, c, d = tab.values[0, 0], tab.values[0, 1], tab.values[1, 0], tab.values[1, 1]
        n = a + b + c + d
        phi = (a*d - b*c) / np.sqrt(max((a+b)*(c+d)*(a+c)*(b+d), 1))
        exp = np.outer(tab.sum(axis=1), tab.sum(axis=0)) / n
        if (exp < 5).any():
            _, p = fisher_exact(tab)
        else:
            _, p, _, _ = chi2_contingency(tab)
        edges.append({"Feature1": x.name, "Feature2": y.name,
                       "Phi": round(phi, 4), "P": p, "Type": cat})

    # Pheno-pheno
    pcols = pheno_df.columns.tolist()
    for i, j in combinations(pcols, 2):
        _test_pair(pheno_df[i], pheno_df[j], "pheno-pheno")

    if gene_df is not None:
        gcols = gene_df.columns.tolist()
        # Gene-gene
        for i, j in combinations(gcols, 2):
            _test_pair(gene_df[i], gene_df[j], "gene-gene")
        # Pheno-gene
        for p in pcols:
            for g in gcols:
                _test_pair(pheno_df[p], gene_df[g], "pheno-gene")

    edf = pd.DataFrame(edges)
    if len(edf) == 0:
        return nx.Graph(), edf

    _, padj, _, _ = multipletests(edf["P"], alpha=alpha, method=fdr_method)
    edf["P_adj"] = padj
    edf["Significant"] = padj < alpha

    G = nx.Graph()
    for _, r in edf[edf["Significant"]].iterrows():
        G.add_edge(r["Feature1"], r["Feature2"],
                   phi=r["Phi"], p=r["P_adj"], edge_type=r["Type"])
    return G, edf
