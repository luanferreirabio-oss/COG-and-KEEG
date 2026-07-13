"""
================================================================================
Script 04 — Metagenomic Data Normalisation and Quality Control
================================================================================
Study: Integration of metagenomic profiles, soil chemical attributes and
       socioenvironmental variables in the determination of floodplain cacao
       productivity in Mocajuba, Pará, Brazil.

Author: [Your name]
Institution: [Your institution]
Date: 2024

Description:
    Reads individual HUMAnN3 tabular output files (one per island, already
    regrouped into KEGG KO and COG categories by humann_regroup_table),
    performs:
        1. Removal of non-informative rows (UNGROUPED, UNMAPPED, UNINTEGRATED)
        2. Removal of species-stratified rows (rows containing "|")
        3. Construction of joint abundance matrices (categories × islands)
        4. Duplicate detection: checks whether any pair of samples has
           ≥ 95% identical values — raises an error if so (catches the
           P3/P4 duplication bug fixed during this study)
        5. Generation of Top-10 tables per island
        6. Export of all outputs

Input:
    Individual HUMAnN3 tabular files (already in CPM units from the Galaxy
    pipeline, or raw RPK — set ALREADY_CPM accordingly):
        Galaxy342-_KEGG_P1_.tabular  ...  Galaxy347-_KEGG_P6_.tabular
        Galaxy336-_COG_P1_.tabular   ...  Galaxy341-_COG_P6_.tabular

Output:
    - Join_KEGG_corrigido.csv
    - Join_COG_corrigido.csv
    - Top10_KEGG_corrigido.csv
    - Top10_COG_corrigido.csv

References:
    Beghini et al. (2021) eLife 10:e65088.
================================================================================
"""

import pandas as pd
import numpy as np
from pathlib import Path

# ── Configuration ─────────────────────────────────────────────────────────────
ALREADY_CPM = True    # Set False if input files are in RPK units

SAMPLES = ["P1", "P2", "P3", "P4", "P5", "P6"]

FILES = {
    "KEGG": {
        "P1": "Galaxy342-_KEGG_P1_.tabular",
        "P2": "Galaxy343-_KEGG_P2_.tabular",
        "P3": "Galaxy344-_KEGG_P3_.tabular",
        "P4": "Galaxy345-_KEGG_P4_.tabular",
        "P5": "Galaxy346-_KEGG_P5_.tabular",
        "P6": "Galaxy347-_KEGG_P6_.tabular",
    },
    "COG": {
        "P1": "Galaxy336-_COG_P1_.tabular",
        "P2": "Galaxy337-_COG_P2_.tabular",
        "P3": "Galaxy338-_COG_P3_.tabular",
        "P4": "Galaxy339-_COG_P4_.tabular",
        "P5": "Galaxy340-_COG_P5_.tabular",
        "P6": "Galaxy341-_COG_P6_.tabular",
    }
}

EXCLUDE_IDS = {"UNMAPPED", "UNGROUPED", "UNINTEGRATED"}
DUP_THRESHOLD = 0.95   # Fraction of identical values to flag as duplicate


def load_and_normalise(filepath: str, already_cpm: bool = True) -> pd.Series:
    """
    Load a single HUMAnN3 tabular file.
    Filters out non-informative rows and species-stratified lines.
    Returns a Series indexed by category ID.
    """
    df = pd.read_csv(filepath, sep="\t", skiprows=1,
                     header=0, names=["ID", "value"])
    # Remove non-informative categories
    df = df[~df["ID"].isin(EXCLUDE_IDS)]
    # Remove species-stratified rows (format: "CATEGORY|species")
    df = df[~df["ID"].str.contains("|", regex=False)]
    df = df.set_index("ID")["value"].astype(float)

    if not already_cpm:
        # Normalise RPK → CPM
        total = df.sum()
        if total > 0:
            df = df / total * 1_000_000

    return df


def check_duplicates(matrix: pd.DataFrame, threshold: float = 0.95) -> None:
    """
    Check all pairwise sample comparisons for suspiciously high identity.
    Raises ValueError if any pair exceeds the threshold.
    """
    cols = matrix.columns.tolist()
    duplicates_found = []
    for i in range(len(cols)):
        for j in range(i + 1, len(cols)):
            pct = (matrix[cols[i]] == matrix[cols[j]]).mean()
            if pct >= threshold:
                duplicates_found.append(
                    f"  {cols[i]} vs {cols[j]}: {pct*100:.1f}% identical"
                )
    if duplicates_found:
        raise ValueError(
            "DUPLICATE SAMPLES DETECTED (≥ 95% identical values):\n"
            + "\n".join(duplicates_found)
            + "\n\nAction required: re-check the HUMAnN3 pipeline for this sample pair."
        )
    print("  No duplicate samples detected. OK.")


# ── Process each database ────────────────────────────────────────────────────
for db_name, file_map in FILES.items():
    print(f"\n{'='*55}")
    print(f"  Processing {db_name}")
    print(f"{'='*55}")

    # Load all samples
    frames = {}
    for sample, filepath in file_map.items():
        if not Path(filepath).exists():
            raise FileNotFoundError(f"File not found: {filepath}")
        frames[sample] = load_and_normalise(filepath, ALREADY_CPM)
        print(f"  {sample}: {len(frames[sample])} categories loaded")

    # Build joint matrix
    joint = (
        pd.DataFrame(frames, columns=SAMPLES)
        .fillna(0)
        .sort_index()
    )
    joint.index.name = "# Gene Family"
    print(f"\n  Joint matrix: {joint.shape[0]} categories × {joint.shape[1]} samples")

    # Duplicate check
    print(f"\n  Checking for duplicate samples...")
    check_duplicates(joint, DUP_THRESHOLD)

    # Summary statistics per sample
    print(f"\n  Summary per sample:")
    for s in SAMPLES:
        total     = joint[s].sum()
        n_nonzero = (joint[s] > 0).sum()
        print(f"    {s}: total CPM = {total:>10,.1f} | "
              f"categories detected = {n_nonzero:>6,}")

    # Top-10 per sample
    top10_list = []
    for s in SAMPLES:
        top = (
            joint[s]
            .sort_values(ascending=False)
            .head(10)
            .reset_index()
        )
        top.columns = ["ID", "CPM"]
        top["Sample"] = s
        top["Rank"]   = range(1, 11)
        top10_list.append(top)
    top10 = pd.concat(top10_list, ignore_index=True)

    # Export
    joint.to_csv(f"Join_{db_name}_corrigido.csv")
    top10.to_csv(f"Top10_{db_name}_corrigido.csv", index=False)
    print(f"\n  Saved: Join_{db_name}_corrigido.csv")
    print(f"  Saved: Top10_{db_name}_corrigido.csv")

print("\nDone.")
