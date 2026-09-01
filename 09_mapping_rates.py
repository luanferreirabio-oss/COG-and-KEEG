"""
================================================================================
Script 09 - Read mapping rates and sequencing-depth confounding
================================================================================
Extracts UNMAPPED and mapped RPK totals from the raw HUMAnN gene-family output
of each island, computes per-sample mapping rates, and tests whether detected
functional richness is confounded by sequencing depth.

Rationale
---------
HUMAnN reports UNMAPPED as the first row of the gene-family table, in RPK units.
Downstream scripts discard this row before CPM normalisation, so the mapping
rate is not recoverable from the joined abundance matrices and must be taken
from the raw per-sample tables.

Rows containing "|" are species-stratified duplicates of the corresponding
unstratified row and are excluded to avoid double counting.

Inputs : raw HUMAnN gene-family .tabular files (one per island)
Output : Supplementary_Table_S8_mapping_rates.csv
================================================================================
"""

# Author: Luan Daniel Silva Ferreira (ORCID 0000-0001-9187-6988)
# Federal University of Para (UFPA), Belem, PA, Brazil
# Repository: https://github.com/luandanbio/floodplain-cacao-mocajuba
# Archived (functional analysis): https://doi.org/10.5281/zenodo.21345125
# Archived (decontamination): https://doi.org/10.5281/zenodo.17498295
# Last updated: 1 September 2026

import itertools
import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

SAMPLES = ["P1", "P2", "P3", "P4", "P5", "P6"]
NAMES = {"P1": "Santana", "P2": "Santaninha", "P3": "Angapijo",
         "P4": "Conceicao", "P5": "Sao Joaquim", "P6": "Tauare"}
PROD = np.array([2000., 600., 450., 1035., 1000., 1500.])

# Galaxy dataset numbers 28-33 correspond to P1-P6 in order.
FILES = {
    "P1": "Galaxy181-_HUMAnN_on_dataset_28_and_113__Gene_families_and_their_abundance_.tabular",
    "P2": "Galaxy190-_HUMAnN_on_dataset_29_and_157__Gene_families_and_their_abundance_.tabular",
    "P3": "Galaxy316-_HUMAnN_on_dataset_30_and_161__Gene_families_and_their_abundance_.tabular",
    "P4": "Galaxy325-_HUMAnN_on_dataset_31_and_165__Gene_families_and_their_abundance_.tabular",
    "P5": "Galaxy222-_HUMAnN_on_dataset_32_and_169__Gene_families_and_their_abundance_.tabular",
    "P6": "Galaxy231-_HUMAnN_on_dataset_33_and_173__Gene_families_and_their_abundance_.tabular",
}


def exact_spearman(x, y):
    """Spearman rho with exact two-sided p-value by complete enumeration (n = 6)."""
    rx, ry = stats.rankdata(x), stats.rankdata(y)

    def corr(a, b):
        am, bm = a - a.mean(), b - b.mean()
        d = np.sqrt((am ** 2).sum() * (bm ** 2).sum())
        return np.dot(am, bm) / d if d > 0 else 0.0

    rho = corr(rx, ry)
    count = sum(1 for perm in itertools.permutations(ry)
                if abs(corr(rx, np.array(perm))) >= abs(rho) - 1e-9)
    return rho, count / math.factorial(len(x))


def _find(name):
    for root in [Path("."), Path("data"), Path("..") / "data",
                 Path("data") / "raw_humann"]:
        if (root / name).exists():
            return root / name
    raise FileNotFoundError(name)


rows = []
for sample, filename in FILES.items():
    d = pd.read_csv(_find(filename), sep="\t", header=0, names=["ID", "RPK"])
    d["RPK"] = pd.to_numeric(d["RPK"], errors="coerce")
    d = d[~d["ID"].astype(str).str.contains("|", regex=False)]   # drop stratified rows

    unmapped = float(d.loc[d.ID == "UNMAPPED", "RPK"].sum())
    total = float(d["RPK"].sum())
    mapped = total - unmapped
    families = int(((d.ID != "UNMAPPED") & (d.RPK > 0)).sum())

    rows.append({"Island": sample, "Name": NAMES[sample],
                 "Total_RPK": int(round(total)),
                 "UNMAPPED_RPK": int(round(unmapped)),
                 "Mapped_RPK": int(round(mapped)),
                 "Mapped_percent": round(100 * mapped / total, 2),
                 "UniRef90_families": families})

s8 = pd.DataFrame(rows)

# Attach richness and CPM totals from the joined matrices
kegg = pd.read_csv(_find("Join_KEGG.csv"), index_col=0)
egg = pd.read_csv(_find("Join_COG_eggNOG.csv"), index_col=0)
s8["KEGG_KOs_detected"] = [int((kegg[p] > 0).sum()) for p in SAMPLES]
s8["eggNOG_IDs_detected"] = [int((egg[p] > 0).sum()) for p in SAMPLES]
s8["KEGG_total_CPM"] = [round(float(kegg[p].sum()), 1) for p in SAMPLES]
s8["eggNOG_total_CPM"] = [round(float(egg[p].sum()), 1) for p in SAMPLES]

s8.to_csv("Supplementary_Table_S8_mapping_rates.csv", index=False)
print(s8.to_string(index=False))

print(f"\nMapping rate: {s8.Mapped_percent.min():.2f}-{s8.Mapped_percent.max():.2f}% "
      f"(mean {s8.Mapped_percent.mean():.2f}%, "
      f"CV {s8.Mapped_percent.std(ddof=1) / s8.Mapped_percent.mean() * 100:.2f}%)")

print("\nSequencing-depth confounding (exact permutation, n = 6):")
for label, y in [("eggNOG identifiers", s8.eggNOG_IDs_detected.values),
                 ("KEGG KOs", s8.KEGG_KOs_detected.values),
                 ("UniRef90 families", s8.UniRef90_families.values)]:
    rho, p = exact_spearman(s8.Total_RPK.values, y)
    print(f"  Total RPK x {label:<20} rho = {rho:+.3f}  p = {p:.4f}")

rho, p = exact_spearman(s8.Mapped_percent.values, PROD)
print(f"\nMapping rate x production:  rho = {rho:+.3f}  p = {p:.4f}")
