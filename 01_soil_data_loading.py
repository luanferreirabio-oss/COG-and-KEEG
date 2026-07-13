"""
================================================================================
Script 01 — Soil Chemistry Data Loading and Cleaning
================================================================================
Study: Integration of metagenomic profiles, soil chemical attributes and
       socioenvironmental variables in the determination of floodplain cacao
       productivity in Mocajuba, Pará, Brazil.

Author: [Your name]
Institution: [Your institution]
Date: 2024

Description:
    Loads raw soil chemistry data (36 samples: 6 islands × 2 depths × 3 replicates),
    standardises island naming conventions (1S → P1, etc.), and exports a clean
    CSV for downstream analyses.

Input:
    - solo_quimica.csv  (36 rows × 22 columns; raw lab output)

Output:
    - solo_quimica_clean.csv
================================================================================
"""

import pandas as pd
import numpy as np

# ── Variables analysed ────────────────────────────────────────────────────────
VARS_QUIM = [
    "pH", "Carbono", "MO", "N", "CN", "P", "Al", "Acidez",
    "Na", "K", "Ca", "Mg", "S", "CTC", "V", "Cu", "Zn", "Mn", "Fe"
]

UNITS = {
    "pH": "dimensionless",
    "Carbono": "g/kg", "MO": "g/kg", "N": "g/kg",
    "CN": "ratio",
    "P": "mg/kg", "Cu": "mg/kg", "Zn": "mg/kg", "Mn": "mg/kg", "Fe": "mg/kg",
    "Al": "cmolc/kg", "Acidez": "cmolc/kg", "Na": "cmolc/kg", "K": "cmolc/kg",
    "Ca": "cmolc/kg", "Mg": "cmolc/kg", "S": "cmolc/kg", "CTC": "cmolc/kg",
    "V": "%"
}

ISLAND_NAMES = {
    "P1": "Santana",
    "P2": "Santaninha",
    "P3": "Angapijó",
    "P4": "Conceição",
    "P5": "São Joaquim",
    "P6": "Tauaré"
}

PRODUCTIVITY = {
    "P1": 2000, "P2": 600, "P3": 450,
    "P4": 1035, "P5": 1000, "P6": 1500
}

# ── Load ──────────────────────────────────────────────────────────────────────
df = pd.read_csv("solo_quimica.csv")

# Standardise island naming: "1S" → "P1"
df["Ponto"] = df["Ponto"].str.extract(r"(\d+)")[0].apply(lambda x: f"P{x}")

# Add island names and productivity
df["Island_name"] = df["Ponto"].map(ISLAND_NAMES)
df["Productivity_kg_yr"] = df["Ponto"].map(PRODUCTIVITY)

# Basic validation
assert df.shape[0] == 36, f"Expected 36 rows, got {df.shape[0]}"
assert set(df["Ponto"].unique()) == {"P1","P2","P3","P4","P5","P6"}, "Missing islands"
assert set(df["Profundidade"].unique()) == {"0-10","10-20"}, "Unexpected depth values"

print(f"Dataset loaded: {df.shape[0]} samples, {df.shape[1]} columns")
print(f"Islands: {sorted(df['Ponto'].unique())}")
print(f"Depths: {sorted(df['Profundidade'].unique())}")
print(f"Replicates per island×depth: {df.groupby(['Ponto','Profundidade']).size().unique()}")

# ── Descriptive statistics by island and depth ────────────────────────────────
desc = df.groupby(["Ponto","Profundidade"])[VARS_QUIM].agg(["mean","std"]).round(3)
print("\nDescriptive statistics (mean ± SD per island × depth):")
print(desc)

# ── Aggregate to island-level mean (for integration with metagenomics) ─────────
island_mean = df.groupby("Ponto")[VARS_QUIM].mean().reset_index()
island_mean["Island_name"] = island_mean["Ponto"].map(ISLAND_NAMES)
island_mean["Productivity_kg_yr"] = island_mean["Ponto"].map(PRODUCTIVITY)
island_mean = island_mean.sort_values("Ponto").reset_index(drop=True)

print("\nIsland-level means (aggregated across depths and replicates):")
print(island_mean[["Ponto","Island_name","Productivity_kg_yr","Ca","Mg","pH","MO"]])

# ── Export ─────────────────────────────────────────────────────────────────────
df.to_csv("solo_quimica_clean.csv", index=False)
island_mean.to_csv("tabela_mestre.csv", index=False)
print("\nFiles saved: solo_quimica_clean.csv, tabela_mestre.csv")
