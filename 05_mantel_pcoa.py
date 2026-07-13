"""
================================================================================
Script 05 — Mantel Tests and PCoA of Functional Profiles
================================================================================
Study: Integration of metagenomic profiles, soil chemical attributes and
       socioenvironmental variables in the determination of floodplain cacao
       productivity in Mocajuba, Pará, Brazil.

Author: [Your name]
Institution: [Your institution]
Date: 2024

Description:
    1. Computes Bray-Curtis dissimilarity matrices for KEGG and COG profiles.
    2. Computes Euclidean distance matrix for soil chemistry (z-score).
    3. Runs Mantel tests (Pearson r; 9,999 permutations) for three comparisons:
           Soil chemistry  ×  KEGG
           Soil chemistry  ×  COG
           KEGG  ×  COG
    4. Performs PCoA (classical double-centering) on KEGG and COG Bray-Curtis
       matrices; projects productivity as a supplementary vector.
    5. Generates and saves the PCoA biplot (Figure 5 in the manuscript).

Input:
    - tabela_mestre.csv            (output of Script 01)
    - Join_KEGG_corrigido.csv      (output of Script 04)
    - Join_COG_corrigido.csv       (output of Script 04)

Output:
    - mantel_results.csv
    - bray_curtis_pairwise.csv
    - fig_PCoA.png   (Figure 5)

References:
    Mantel (1967) Cancer Res. 27:209–220.
    Legendre & Legendre (2012) Numerical Ecology, 3rd ed. Elsevier.
================================================================================
"""

import pandas as pd
import numpy as np
from scipy import stats
from sklearn.preprocessing import StandardScaler
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ── Configuration ─────────────────────────────────────────────────────────────
RANDOM_SEED  = 42
N_PERM       = 9999
SAMPLES      = ["P1","P2","P3","P4","P5","P6"]
ISLAND_NAMES = ["Santana","Santaninha","Angapijó","Conceição","S. Joaquim","Tauaré"]
PROD_VALS    = np.array([2000, 600, 450, 1035, 1000, 1500])

VARS_QUIM = [
    "pH","Carbono","MO","N","CN","P","Al","Acidez",
    "Na","K","Ca","Mg","S","CTC","V","Cu","Zn","Mn","Fe"
]

# ── Load data ─────────────────────────────────────────────────────────────────
mestre = pd.read_csv("tabela_mestre.csv")
kegg   = pd.read_csv("Join_KEGG_corrigido.csv", index_col=0)
cog    = pd.read_csv("Join_COG_corrigido.csv",  index_col=0)

# ── Distance / dissimilarity functions ───────────────────────────────────────
def euclidean_dist(Z: np.ndarray) -> np.ndarray:
    """Euclidean distance matrix from z-score standardised matrix."""
    n = Z.shape[0]
    D = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            d = np.sqrt(((Z[i] - Z[j]) ** 2).sum())
            D[i, j] = D[j, i] = d
    return D


def bray_curtis(df: pd.DataFrame, samples: list) -> np.ndarray:
    """Bray-Curtis dissimilarity matrix between samples (columns)."""
    arr = df[samples].values.T    # shape: n_samples × n_features
    n   = arr.shape[0]
    D   = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            num = np.abs(arr[i] - arr[j]).sum()
            den = arr[i].sum() + arr[j].sum()
            D[i, j] = D[j, i] = num / den if den > 0 else 0
    return D


def vec_upper(D: np.ndarray) -> np.ndarray:
    """Upper triangle (excluding diagonal) as a vector."""
    return D[np.triu_indices(D.shape[0], k=1)]


def mantel_test(D1: np.ndarray, D2: np.ndarray,
                n_perm: int = 9999, seed: int = 42) -> tuple:
    """
    Mantel test: Pearson r between upper triangles of two distance matrices.
    p-value by random permutation of rows/columns of D1.
    Returns (r_observed, p_value).
    """
    rng   = np.random.default_rng(seed)
    v1, v2 = vec_upper(D1), vec_upper(D2)
    r_obs, _ = stats.pearsonr(v1, v2)
    n     = D1.shape[0]
    count = 0
    for _ in range(n_perm):
        perm  = rng.permutation(n)
        D1p   = D1[np.ix_(perm, perm)]
        r_p, _ = stats.pearsonr(vec_upper(D1p), v2)
        if r_p >= r_obs:
            count += 1
    p_val = (count + 1) / (n_perm + 1)
    return r_obs, p_val


def pcoa(D: np.ndarray) -> tuple:
    """
    Classical Principal Coordinates Analysis (PCoA) via double-centering.
    Returns (coordinates, proportion of variance explained) for positive eigenvalues.
    """
    n    = D.shape[0]
    H    = np.eye(n) - np.ones((n, n)) / n
    B    = -0.5 * H @ (D ** 2) @ H
    eigvals, eigvecs = np.linalg.eigh(B)
    idx     = np.argsort(eigvals)[::-1]
    eigvals = eigvals[idx]; eigvecs = eigvecs[:, idx]
    pos     = eigvals > 0
    coords  = eigvecs[:, pos] * np.sqrt(eigvals[pos])
    var_exp = eigvals[pos] / eigvals[pos].sum()
    return coords, var_exp


# ── Compute distance matrices ─────────────────────────────────────────────────
Z       = StandardScaler().fit_transform(mestre[VARS_QUIM].values)
D_quim  = euclidean_dist(Z)
D_kegg  = bray_curtis(kegg, SAMPLES)
D_cog   = bray_curtis(cog,  SAMPLES)

# ── Pairwise Bray-Curtis table ────────────────────────────────────────────────
bc_rows = []
for i in range(len(SAMPLES)):
    for j in range(i + 1, len(SAMPLES)):
        bc_rows.append({
            "Sample_A":  SAMPLES[i],
            "Sample_B":  SAMPLES[j],
            "BC_KEGG":   round(D_kegg[i, j], 4),
            "BC_COG":    round(D_cog[i, j],  4),
        })
bc_df = pd.DataFrame(bc_rows)
print("Pairwise Bray-Curtis dissimilarities:")
print(bc_df.to_string(index=False))
bc_df.to_csv("bray_curtis_pairwise.csv", index=False)

# ── Mantel tests ──────────────────────────────────────────────────────────────
print(f"\nMantel tests (Pearson r; {N_PERM} permutations):")

mantel_results = []
for label, D1, D2 in [
    ("Soil chemistry × KEGG", D_quim, D_kegg),
    ("Soil chemistry × COG",  D_quim, D_cog),
    ("KEGG × COG",            D_kegg, D_cog),
]:
    r, p = mantel_test(D1, D2, n_perm=N_PERM, seed=RANDOM_SEED)
    sig  = "**" if p < 0.01 else ("*" if p < 0.05 else "ns")
    print(f"  {label:<30}  r = {r:.3f}   p = {p:.4f}  {sig}")
    mantel_results.append({"Comparison": label, "r": round(r, 3),
                            "p": round(p, 4), "sig": sig})

pd.DataFrame(mantel_results).to_csv("mantel_results.csv", index=False)

# ── PCoA ──────────────────────────────────────────────────────────────────────
cmap_p   = plt.cm.RdYlGn
norm_p   = (PROD_VALS - PROD_VALS.min()) / (PROD_VALS.max() - PROD_VALS.min())
colors_p = [cmap_p(v) for v in norm_p]

fig, axes = plt.subplots(1, 2, figsize=(17, 8))
fig.subplots_adjust(wspace=0.40)

for ax, (nome, D) in zip(axes, [("KEGG", D_kegg), ("COG", D_cog)]):
    coords, vexp = pcoa(D)
    c1, c2 = coords[:, 0], coords[:, 1]

    rp1, _ = stats.pearsonr(PROD_VALS, c1)
    rp2, _ = stats.pearsonr(PROD_VALS, c2)
    sc2    = max(abs(c1).max(), abs(c2).max()) * 0.72

    print(f"\n  PCoA {nome}: PC1 {vexp[0]*100:.1f}%  PC2 {vexp[1]*100:.1f}%")
    print(f"    Productivity × PC1: r = {rp1:.3f} | PC2: r = {rp2:.3f}")

    for i in range(6):
        ax.scatter(c1[i], c2[i], c=[colors_p[i]], s=270, zorder=5,
                   edgecolors="k", linewidths=1.2)

    y_mid = np.median(c2)
    yr    = abs(c2.max() - c2.min())
    xr    = abs(c1.max() - c1.min())
    for i, p in enumerate(SAMPLES):
        ty = c2[i] + yr * (0.22 if c2[i] >= y_mid else -0.26)
        tx = c1[i] + xr * 0.10
        ax.annotate(
            f"{p} ({ISLAND_NAMES[i]})\n{int(PROD_VALS[i])} kg yr⁻¹",
            xy=(c1[i], c2[i]), xytext=(tx, ty),
            fontsize=8.5, fontweight="bold", color="#1a1a2e",
            arrowprops=dict(arrowstyle="->", color="#777", lw=0.7, mutation_scale=9),
            bbox=dict(fc="white", ec="none", alpha=0.8, pad=1.5), zorder=8
        )

    ax.annotate("", xy=(rp1*sc2, rp2*sc2), xytext=(0, 0),
                arrowprops=dict(arrowstyle="->", color="#e63946",
                                lw=2.4, mutation_scale=14))
    ax.text(
        rp1*sc2*1.20, rp2*sc2*1.20,
        f"Productivity\nPC1: r = {rp1:.2f}\nPC2: r = {rp2:.2f}",
        fontsize=8.5, color="#e63946", fontweight="bold", ha="center",
        bbox=dict(fc="white", ec="#e63946", alpha=0.85, lw=1, pad=1.5), zorder=9
    )
    ax.axhline(0, color="grey", lw=0.5, ls="--", alpha=0.5)
    ax.axvline(0, color="grey", lw=0.5, ls="--", alpha=0.5)
    xm = abs(c1.max()-c1.min())*0.50; ym = abs(c2.max()-c2.min())*0.55
    ax.set_xlim(c1.min()-xm, c1.max()+xm)
    ax.set_ylim(c2.min()-ym, c2.max()+ym)
    ax.set_xlabel(f"PC1 ({vexp[0]*100:.1f}% of variance)", fontsize=11)
    ax.set_ylabel(f"PC2 ({vexp[1]*100:.1f}% of variance)", fontsize=11)
    ax.set_title(
        f"PCoA (Bray-Curtis) — {nome} functional profiles\n(n = 6 islands)",
        fontsize=11, fontweight="bold", pad=12
    )
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

sm = plt.cm.ScalarMappable(
    cmap=cmap_p, norm=plt.Normalize(PROD_VALS.min(), PROD_VALS.max())
)
sm.set_array([])
plt.colorbar(sm, ax=axes, shrink=0.48, pad=0.04, label="Productivity (kg yr⁻¹)")
plt.suptitle(
    "Ordination of functional profiles of floodplain islands\n"
    "(KEGG and COG — Bray-Curtis dissimilarity)",
    fontsize=12, fontweight="bold", y=1.01
)
plt.savefig("fig_PCoA.png", dpi=220, bbox_inches="tight")
plt.close()

print("\nFiles saved: mantel_results.csv, bray_curtis_pairwise.csv, fig_PCoA.png")
