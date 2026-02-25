"""MDR classification and prevalence estimation (NA-aware).

Key design choice (defensibility): missingness is *not* treated as susceptibility.
If an isolate was not tested against any drug in a class, the class call is NA.

Intended input
--------------
The phenotype input is assumed to be **binary** per drug: 1=resistant, 0=susceptible,
with missing values allowed (NA/empty cells). If the input contains numeric MICs,
convert to binary first (e.g., EUCAST ECOFF / clinical breakpoint decision) before
using this tool.
"""

from __future__ import annotations

from typing import Dict, List, Literal

import numpy as np
import pandas as pd
from scipy.stats import beta as beta_dist


def _to_binary(df: pd.DataFrame) -> pd.DataFrame:
    """Coerce numeric-ish values to {0,1} Int8 with NA preserved."""
    out = df.copy()
    for c in out.columns:
        if pd.api.types.is_bool_dtype(out[c]):
            out[c] = out[c].astype("Int8")
        elif pd.api.types.is_numeric_dtype(out[c]):
            out[c] = out[c].where(out[c].isna(), (out[c] > 0).astype("Int8"))
        else:
            num = pd.to_numeric(out[c], errors="coerce")
            if num.notna().any():
                out[c] = num.where(num.isna(), (num > 0).astype("Int8"))
    for c in out.columns:
        if out[c].dtype == "int64":
            out[c] = out[c].astype("Int8")
    return out


def build_class_matrix(
    pheno_df: pd.DataFrame,
    class_to_drugs: Dict[str, List[str]],
    class_rule: Literal["any", "all", "majority"] = "any",
) -> pd.DataFrame:
    """Map drug-level phenotypes to antibiotic-class phenotypes with NA-aware logic.

    Parameters
    ----------
    pheno_df : pd.DataFrame
        DataFrame indexed by isolate IDs. Columns are drugs (binary 0/1 with NA allowed).
    class_to_drugs : dict
        Mapping: class_name -> list of drug column names.
    class_rule : {"any", "all", "majority"}, default "any"
        Rule for classifying a class as resistant from individual drug calls:

        ``"any"`` *(default, conservative)*
            Class = 1 if at least one tested drug in the class is resistant.
            Class = 0 only if all tested drugs are susceptible and all were tested.
            This is the WHO/EUCAST-aligned definition and is the most sensitive
            (minimises missed MDR).

        ``"all"``
            Class = 1 only if *all* tested drugs in the class are resistant.
            Class = 0 if at least one drug is susceptible and all were tested.
            More specific; appropriate when pan-class resistance is the question
            (e.g., extended-spectrum beta-lactamases covering all penicillins).

        ``"majority"``
            Class = 1 if strictly more than half of tested drugs are resistant.
            Class = 0 if at most half of tested drugs are resistant and all tested.
            Intermediate between ``"any"`` and ``"all"``.

        In all rules: NA is returned if the tested information is insufficient
        to make a definitive call (i.e., the classification would change if
        untested drugs were known).

    Returns
    -------
    pd.DataFrame indexed as pheno_df with one column per class, values in {0, 1, NA}.

    Examples
    --------
    >>> build_class_matrix(pheno_df, class_map, class_rule="any")   # default WHO definition
    >>> build_class_matrix(pheno_df, class_map, class_rule="majority")  # stricter rule
    """
    if pheno_df.index is None:
        raise ValueError("pheno_df must be indexed by isolate IDs")

    _valid_rules = ("any", "all", "majority")
    if class_rule not in _valid_rules:
        raise ValueError(
            f"class_rule must be one of {_valid_rules}, got {class_rule!r}"
        )

    ph = _to_binary(pheno_df)
    out = pd.DataFrame(index=ph.index)

    for cls, drugs in class_to_drugs.items():
        present = [d for d in drugs if d in ph.columns]
        if not present:
            continue

        vals = ph[present]
        n_drugs = len(present)

        n_resistant = (vals == 1).sum(axis=1, skipna=True).astype(int)
        n_susceptible = (vals == 0).sum(axis=1, skipna=True).astype(int)
        n_tested = n_resistant + n_susceptible
        all_tested = n_tested == n_drugs

        s = pd.Series(pd.NA, index=ph.index, dtype="Int8")

        if class_rule == "any":
            # Resistant if any drug tested is resistant
            any_pos = n_resistant >= 1
            all_neg = all_tested & (n_resistant == 0)
            s.loc[any_pos] = 1
            s.loc[all_neg & (~any_pos)] = 0

        elif class_rule == "all":
            # Resistant only if all tested drugs are resistant
            # Susceptible only if at least one tested drug is susceptible and all tested
            all_pos = all_tested & (n_resistant == n_drugs)
            any_neg = all_tested & (n_susceptible >= 1)
            s.loc[all_pos] = 1
            s.loc[any_neg & (~all_pos)] = 0
            # Uncertain: not all tested, or some resistant and some susceptible

        elif class_rule == "majority":
            # Resistant if >50% of tested drugs are resistant
            majority_pos = (n_tested > 0) & (n_resistant > n_tested / 2)
            majority_neg = all_tested & (n_resistant <= n_tested / 2)
            s.loc[majority_pos] = 1
            s.loc[majority_neg & (~majority_pos)] = 0

        out[cls] = s

    return out


def identify_mdr(class_df: pd.DataFrame, threshold: int = 3) -> pd.Series:
    """MDR indicator with explicit uncertainty under missingness.

    Returns:
    - True  : definitely MDR (min_resistant >= threshold)
    - False : definitely not MDR (max_resistant < threshold)
    - NA    : uncertain due to missing class calls
    """
    min_res = (class_df == 1).sum(axis=1, skipna=True).astype(int)
    n_missing = class_df.isna().sum(axis=1).astype(int)
    max_res = min_res + n_missing

    out = pd.Series(False, index=class_df.index, dtype="object")
    out[min_res >= threshold] = True
    out[max_res < threshold] = False
    out[(min_res < threshold) & (max_res >= threshold)] = pd.NA
    return out


def mdr_spectrum(class_df: pd.DataFrame, threshold: int = 3) -> pd.DataFrame:
    """MDR/XDR/PDR spectrum with NA-aware bounds."""
    n_classes = int(class_df.shape[1])
    min_res = (class_df == 1).sum(axis=1, skipna=True).astype(int)
    n_missing = class_df.isna().sum(axis=1).astype(int)
    max_res = min_res + n_missing

    def cat(s: int) -> str:
        if s >= n_classes:
            return "PDR"
        if s >= n_classes - 2:
            return "XDR"
        if s >= threshold:
            return "MDR"
        return "susceptible"

    cat_def = min_res.apply(cat)
    cat_pos = max_res.apply(cat)

    return pd.DataFrame(
        {
            "N_Classes_Resistant_Min": min_res,
            "N_Classes_Resistant_Max": max_res,
            # Backward-compatible alias (legacy): Category == definite classification
            "Category": cat_def,
            "Category_Definite": cat_def,
            "Category_Possible": cat_pos,
            "N_Classes_Total": n_classes,
            "N_Classes_Missing": n_missing,
        }
    )


def bayesian_prevalence(
    positives: int,
    total: int,
    prior_a: float = 0.5,
    prior_b: float = 0.5,
    credible: float = 0.95,
) -> Dict[str, float]:
    """Bayesian prevalence with Jeffrey's prior Beta(0.5,0.5)."""
    if total <= 0:
        return {
            "Posterior_Mean": np.nan,
            "Posterior_Mode": np.nan,
            "CI_Lo": np.nan,
            "CI_Hi": np.nan,
        }

    pa = prior_a + positives
    pb = prior_b + (total - positives)
    mean = pa / (pa + pb)
    mode = (pa - 1) / (pa + pb - 2) if pa > 1 and pb > 1 else mean
    tail = (1 - credible) / 2

    return {
        "Posterior_Mean": float(round(mean, 6)),
        "Posterior_Mode": float(round(mode, 6)),
        "CI_Lo": float(round(float(beta_dist.ppf(tail, pa, pb)), 6)),
        "CI_Hi": float(round(float(beta_dist.ppf(1 - tail, pa, pb)), 6)),
    }


def prevalence_table(class_df: pd.DataFrame, credible: float = 0.95) -> pd.DataFrame:
    """Bayesian prevalence for each drug class (conditioning on being tested)."""
    rows = []
    n_all = int(len(class_df))
    for col in class_df.columns:
        obs = int(class_df[col].notna().sum())
        pos = int((class_df[col] == 1).sum(skipna=True))
        miss = n_all - obs
        bp = bayesian_prevalence(pos, obs, credible=credible)
        rows.append(
            {
                "Class": col,
                "Positive": pos,
                "Total_Observed": obs,
                "Total_All": n_all,
                "Missing": miss,
                **bp,
            }
        )
    return pd.DataFrame(rows)


def mdr_probability(class_df: pd.DataFrame, threshold: int = 3) -> pd.DataFrame:
    """Probabilistic MDR under missing class calls."""
    n_classes = int(class_df.shape[1])
    p_class = (
        class_df.mean(axis=0, skipna=True).astype(float).fillna(0.0).clip(0.0, 1.0)
    )

    min_res = (class_df == 1).sum(axis=1, skipna=True).astype(int)
    missing_mask = class_df.isna()

    out_rows = []
    for sid in class_df.index:
        miss_cols = missing_mask.loc[sid]
        probs = p_class[miss_cols].to_numpy(dtype=float)

        dist = np.array([1.0], dtype=float)
        for p in probs:
            dist = np.convolve(dist, np.array([1.0 - p, p], dtype=float))

        need_mdr = max(0, threshold - int(min_res.loc[sid]))
        need_xdr = max(0, (n_classes - 2) - int(min_res.loc[sid]))
        need_pdr = max(0, n_classes - int(min_res.loc[sid]))

        p_mdr = float(dist[need_mdr:].sum()) if need_mdr < len(dist) else 0.0
        p_xdr = float(dist[need_xdr:].sum()) if need_xdr < len(dist) else 0.0
        p_pdr = float(dist[need_pdr:].sum()) if need_pdr < len(dist) else 0.0

        out_rows.append(
            {
                "Strain_ID": sid,
                "Observed_Resistant_Classes": int(min_res.loc[sid]),
                "Missing_Classes": int(miss_cols.sum()),
                "P_MDR": p_mdr,
                "P_XDR": p_xdr,
                "P_PDR": p_pdr,
            }
        )
    return pd.DataFrame(out_rows)
