"""
================================================================================
Script 03 — Soil Chemistry × Productivity: Spearman Correlations and PCA
================================================================================
Study: Integration of metagenomic profiles, soil chemical attributes and
       socioenvironmental variables in the determination of floodplain cacao
       productivity in Mocajuba, Pará, Brazil.

Author: [Your name]
Institution: [Your institution]
Date: 2024

Description:
    1. Computes Spearman rank correlations between each of the 19 soil chemical
       variables and cacao productivity (n = 6 islands).
    2. Performs PCA on z-score standardised soil chemistry; productivity is
       projected as a supplementary variable (Pearson r with PCA axes).
    3. Generates and saves the PCA biplot (Figure 3 in the manuscript).

Input:
    - tabela_mestre.csv   (output of Script 01; 1 row per island, 19 chemical vars)

Output:
    - spearman_quimica_produtividade.csv
    - fig_PCA.png   (Figure 3)
================================================================================
"""

import pandas as pd
import numpy as np
from scipy import stats
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ── Configuration ─────────────────────────────────────────────────────────────
VARS_QUIM = [
    "pH", "Carbono", "MO", "N", "CN", "P", "Al", "Acidez",
    "Na", "K", "Ca", "Mg", "S", "CTC", "V", "Cu", "Zn", "Mn", "Fe"
]
PONTOS       = ["P1", "P2", "P3", "P4", "P5", "P6"]
ISLAND_NAMES = ["Santana","Santaninha","Angapijó","Conceição","S. Joaquim","Tauaré"]
PROD_VALS    = np.array([2000, 600, 450, 1035, 1000, 1500])

# ── Load ──────────────────────────────────────────────────────────────────────
mestre = pd.read_csv("tabela_mestre.csv")

# ── 1. Spearman correlations ──────────────────────────────────────────────────
rows = []
for v in VARS_QUIM:
    rho, p = stats.spearmanr(mestre[v], PROD_VALS)
    sig = "**" if p < 0.01 else ("*" if p < 0.05 else "ns")
    rows.append({
        "variable":  v,
        "rho":       round(rho, 3),
        "p":         round(p, 4),
        "sig":       sig,
        "direction": "Positive" if rho > 0 else "Negative"
    })

spearman_df = (
    pd.DataFrame(rows)
    .sort_values("rho", key=abs, ascending=False)
    .reset_index(drop=True)
)

print("Spearman correlations — soil chemistry vs. productivity (n = 6):")
print(spearman_df.to_string(index=False))
spearman_df.to_csv("spearman_quimica_produtividade.csv", index=False)

# ── 2. PCA ────────────────────────────────────────────────────────────────────
X   = mestre[VARS_QUIM].values
Z   = StandardScaler().fit_transform(X)
pca = PCA()
sc  = pca.fit_transform(Z)
ve  = pca.explained_variance_ratio_
pc1, pc2 = sc[:, 0], sc[:, 1]

# Biplot arrows: Pearson r between each standardised variable and each PC axis
arrows = {
    v: (
        stats.pearsonr(Z[:, i], pc1)[0],
        stats.pearsonr(Z[:, i], pc2)[0]
    )
    for i, v in enumerate(VARS_QUIM)
}

# Productivity as supplementary vector
rp1, _ = stats.pearsonr(PROD_VALS, pc1)
rp2, _ = stats.pearsonr(PROD_VALS, pc2)

print(f"\nPCA variance explained:")
for i, v in enumerate(ve[:4]):
    print(f"  PC{i+1}: {v*100:.1f}%  (cumulative: {ve[:i+1].sum()*100:.1f}%)")
print(f"\nProductivity correlation with PCA axes:")
print(f"  PC1: r = {rp1:.3f}  |  PC2: r = {rp2:.3f}")

# ── 3. PCA Biplot ─────────────────────────────────────────────────────────────
cmap_p   = plt.cm.RdYlGn
norm_p   = (PROD_VALS - PROD_VALS.min()) / (PROD_VALS.max() - PROD_VALS.min())
colors_p = [cmap_p(v) for v in norm_p]
SCALE    = 2.0
y_mid    = np.median(pc2)
y_rng    = abs(pc2.max() - pc2.min())
x_rng    = abs(pc1.max() - pc1.min())

fig, ax = plt.subplots(figsize=(11, 9))

# Sample points
for i in range(6):
    ax.scatter(pc1[i], pc2[i], c=[colors_p[i]], s=260, zorder=5,
               edgecolors="k", linewidths=1.1)

# Point labels
for i, p in enumerate(PONTOS):
    ty = pc2[i] + y_rng * (0.19 if pc2[i] >= y_mid else -0.23)
    tx = pc1[i] + x_rng * 0.09
    ax.annotate(
        f"{p} ({ISLAND_NAMES[i]})\n{int(PROD_VALS[i])} kg yr⁻¹",
        xy=(pc1[i], pc2[i]), xytext=(tx, ty),
        fontsize=8.5, fontweight="bold", color="#1a1a2e",
        arrowprops=dict(arrowstyle="->", color="#777", lw=0.7, mutation_scale=9),
        bbox=dict(fc="white", ec="none", alpha=0.8, pad=1.5), zorder=8
    )

# Variable vectors (biplot)
VARS_PLOT = ["Ca","Mg","P","N","MO","pH","Al","Fe","K","V","Mn","Cu"]
for v in VARS_PLOT:
    dx, dy = arrows[v]
    ax.annotate("", xy=(dx*SCALE, dy*SCALE), xytext=(0, 0),
                arrowprops=dict(arrowstyle="->", color="#457b9d",
                                lw=1.3, mutation_scale=10))
    ax.text(dx*SCALE*1.09, dy*SCALE*1.09, v, fontsize=8.5, color="#457b9d",
            ha="left" if dx >= 0 else "right", fontweight="semibold")

# Productivity supplementary vector
ax.annotate("", xy=(rp1*SCALE, rp2*SCALE), xytext=(0, 0),
            arrowprops=dict(arrowstyle="->", color="#e63946",
                            lw=2.5, mutation_scale=14))
ax.text(rp1*SCALE*1.13, rp2*SCALE*1.13, "Productivity",
        fontsize=10, color="#e63946", fontweight="bold",
        bbox=dict(fc="white", ec="#e63946", alpha=0.9, lw=1.2, pad=2.5))

ax.axhline(0, color="grey", lw=0.5, ls="--", alpha=0.5)
ax.axvline(0, color="grey", lw=0.5, ls="--", alpha=0.5)
ax.set_xlabel(f"PC1 ({ve[0]*100:.1f}% of variance)", fontsize=11)
ax.set_ylabel(f"PC2 ({ve[1]*100:.1f}% of variance)", fontsize=11)
ax.set_title(
    "Principal Component Analysis of soil chemical attributes\n"
    "(n = 6 floodplain islands; red arrow = productivity vector)",
    fontsize=11, fontweight="bold", pad=12
)
sm = plt.cm.ScalarMappable(
    cmap=cmap_p, norm=plt.Normalize(PROD_VALS.min(), PROD_VALS.max())
)
sm.set_array([])
plt.colorbar(sm, ax=ax, shrink=0.50, pad=0.03, label="Productivity (kg yr⁻¹)")
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.set_xlim(-3.9, 3.9)
ax.set_ylim(-3.9, 3.9)
plt.tight_layout()
plt.savefig("fig_PCA.png", dpi=220, bbox_inches="tight")
plt.close()

print("\nFiles saved: spearman_quimica_produtividade.csv, fig_PCA.png")
