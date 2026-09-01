"""
================================================================================
Script 07 - Exact permutation inference (n = 6) and COG functional-category
            aggregation
================================================================================
Replaces the asymptotic Spearman p-values used in earlier scripts, which are
invalid at n = 6, and re-derives all COG annotations from the authoritative
NCBI COG-2024 release (Galperin et al. 2025, Nucleic Acids Res 53:D356-D363,
https://doi.org/10.1093/nar/gkae983).

Reference tables required (see data/reference_COG/):
    wget https://ftp.ncbi.nih.gov/pub/COG/COG2024/data/cog-24.def.tab
    wget https://ftp.ncbi.nih.gov/pub/COG/COG2024/data/cog-24.fun.tab
Access date recorded in data/reference_COG/PROVENANCE.txt: 25 August 2026.

Rationale
---------
1) With n = 6 the permutation space is finite (6! = 720). The exact two-sided
   p-value is obtainable by complete enumeration and the smallest attainable
   value is 2/720 = 0.00278. scipy.stats.spearmanr returns p = 0.0 for
   |rho| = 1 (the t-approximation diverges), which propagates into the
   Benjamini-Hochberg step as q = 0 and produces spurious "FDR-confirmed" hits.
2) Permuting the productivity vector against the observed abundance vector
   handles ties (zero-inflated CPM) correctly.
3) COG identifiers are orthologous groups, not functional categories. The
   25 COG functional categories are obtained from cog-24.def.tab; where a COG
   carries multiple category letters, the first letter is retained.

Inputs : data/Join_KEGG.csv, data/Join_COG_eggNOG.csv,
         data/reference_COG/cog-24.def.tab, data/reference_COG/cog-24.fun.tab
Outputs: KEGG_exact_correlations.csv, COG_exact_correlations.csv,
         COG_functional_categories_CPM.csv, COG_functional_categories_stats.csv
================================================================================
"""

# Author: Luan Daniel Silva Ferreira (ORCID 0000-0001-9187-6988)
# Federal University of Para (UFPA), Belem, PA, Brazil
# Repository: https://github.com/luandanbio/floodplain-cacao-mocajuba
# Archived (functional analysis): https://doi.org/10.5281/zenodo.21345125
# Archived (decontamination): https://doi.org/10.5281/zenodo.17498295
# Last updated: 31 August 2026


import itertools
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

SAMPLES = ["P1", "P2", "P3", "P4", "P5", "P6"]
PROD = np.array([2000., 600., 450., 1035., 1000., 1500.])


def _resolve(primary: str, fallback: str) -> Path:
    """Locate a COG reference table in the expected repository layout."""
    roots = [Path("."), Path("data") / "reference_COG",
             Path("..") / "data" / "reference_COG", Path("referencia_COG")]
    for root in roots:
        for name in (primary, fallback):
            candidate = root / name
            if candidate.exists():
                return candidate
    raise FileNotFoundError(
        f"Neither {primary} (NCBI) nor {fallback} (COGclassifier wheel) found. "
        f"Download from https://ftp.ncbi.nih.gov/pub/COG/COG2024/data/")


COG_DEF = _resolve("cog-24.def.tab", "cog_definition.tsv")
COG_FUN = _resolve("cog-24.fun.tab", "cog_func_category.tsv")


# -- Exact permutation machinery ---------------------------------------------
_ry = stats.rankdata(PROD)
_ry = (_ry - _ry.mean()) / np.sqrt(((_ry - _ry.mean()) ** 2).sum())
PERM = np.array([_ry[list(p)] for p in itertools.permutations(range(6))])  # 720 x 6


def exact_spearman(matrix: np.ndarray):
    """Spearman rho and exact two-sided permutation p-value (720 permutations)."""
    R = np.apply_along_axis(stats.rankdata, 1, matrix)
    R = R - R.mean(axis=1, keepdims=True)
    sd = np.sqrt((R ** 2).sum(axis=1))
    sd[sd == 0] = np.nan
    R = R / sd[:, None]
    allr = R @ PERM.T
    rho = allr[:, 0]
    p = (np.abs(allr) >= np.abs(rho)[:, None] - 1e-9).mean(axis=1)
    return rho, p


def benjamini_hochberg(p: np.ndarray) -> np.ndarray:
    p = np.asarray(p, float)
    m = len(p)
    order = np.argsort(p)
    q = p[order] * m / (np.arange(m) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty(m)
    out[order] = np.minimum(q, 1.0)
    return out


# -- Load ---------------------------------------------------------------------
def _find(name):
    for root in [Path("."), Path("data"), Path("..") / "data"]:
        if (root / name).exists():
            return root / name
    raise FileNotFoundError(name)


kegg = pd.read_csv(_find("Join_KEGG.csv"), index_col=0)
eggnog = pd.read_csv(_find("Join_COG_eggNOG.csv"), index_col=0)
cog = eggnog[[str(i).startswith("COG") for i in eggnog.index]]

defs = pd.read_csv(COG_DEF, sep="\t", header=None, dtype=str, engine="python",
                   on_bad_lines="skip",
                   names=["COG", "LETTERS", "NAME", "GENE", "PATHWAY",
                          "PMID", "PDB", "x8", "x9", "x10"]).set_index("COG")
funcat = pd.read_csv(COG_FUN, sep="\t", header=None, dtype=str,
                     names=["LETTER", "GROUP", "COLOR", "DESCRIPTION"]).set_index("LETTER")

print(f"eggNOG regrouped table: {eggnog.shape[0]} identifiers "
      f"(COG = {len(cog)}; the remainder are ENOG/arCOG/KOG)")
print(f"KEGG table: {kegg.shape[0]} KOs\n")


# -- Ortholog-level correlations ---------------------------------------------
def run(df, label, annotate=False):
    rho, p = exact_spearman(df[SAMPLES].to_numpy(float))
    t = pd.DataFrame({"ID": df.index, "rho": rho.round(4),
                      "p_exact_two_sided": p.round(5),
                      "mean_CPM": df[SAMPLES].mean(axis=1).round(3).values})
    t = t.dropna(subset=["rho"])
    t["q_BH"] = benjamini_hochberg(t["p_exact_two_sided"].values).round(4)
    if annotate:
        t["COG_name"] = [defs.loc[i, "NAME"] if i in defs.index else "" for i in t.ID]
        t["gene"] = [defs.loc[i, "GENE"] if i in defs.index else "" for i in t.ID]
        t["category"] = [str(defs.loc[i, "LETTERS"])[0] if i in defs.index else "" for i in t.ID]
        t["category_description"] = [funcat.loc[c, "DESCRIPTION"]
                                     if c in funcat.index else "" for c in t.category]
    n_sig = int((t["p_exact_two_sided"] <= 0.05).sum())
    print(f"{label}: {len(t)} tested | p <= 0.05: {n_sig} "
          f"(expected by chance {round(len(t)*0.03333)}) | "
          f"q_BH < 0.05: {int((t.q_BH < 0.05).sum())} | min q = {t.q_BH.min():.3f}")
    return t.sort_values("p_exact_two_sided")


run(kegg, "KEGG KO").to_csv("KEGG_exact_correlations.csv", index=False)
run(cog, "COG", annotate=True).to_csv("COG_exact_correlations.csv", index=False)


# -- Functional-category aggregation -----------------------------------------
letter = pd.Series({i: str(defs.loc[i, "LETTERS"])[0]
                    for i in cog.index if i in defs.index})
print(f"\nCOG identifiers mapped to a functional category: "
      f"{len(letter)}/{len(cog)} ({100*len(letter)/len(cog):.1f}%)")

agg = cog.loc[letter.index, SAMPLES].groupby(letter).sum()
agg = agg.reindex([l for l in funcat.index if l in agg.index])
rel = (100 * agg / agg.sum()).round(3)

rho_c, p_c = exact_spearman(agg.to_numpy(float))
cat = pd.DataFrame({
    "letter": agg.index,
    "description": [funcat.loc[l, "DESCRIPTION"] for l in agg.index],
    "group": [funcat.loc[l, "GROUP"] for l in agg.index],
    "mean_CPM": agg.mean(axis=1).round(1).values,
    "mean_percent": rel.mean(axis=1).round(2).values,
    "rho": rho_c.round(3),
    "p_exact_two_sided": p_c.round(4),
})
cat["q_BH"] = benjamini_hochberg(cat["p_exact_two_sided"].values).round(4)
cat = cat.sort_values("mean_percent", ascending=False)

print(f"\nFunctional categories: {len(cat)} tested; "
      f"minimum attainable q with m = {len(cat)} is {0.002778*len(cat):.3f}\n")
print(cat.to_string(index=False))

agg.round(2).to_csv("COG_functional_categories_CPM.csv")
cat.to_csv("COG_functional_categories_stats.csv", index=False)
print("\nDone.")
