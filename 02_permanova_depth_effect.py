"""
================================================================================
Script 02 — Soil Depth Effect: PERMANOVA and Paired Tests
================================================================================
Study: Integration of metagenomic profiles, soil chemical attributes and
       socioenvironmental variables in the determination of floodplain cacao
       productivity in Mocajuba, Pará, Brazil.

Author: [Your name]
Institution: [Your institution]
Date: 2024

Description:
    Tests whether sampling depth (0–10 cm vs. 10–20 cm) significantly affects
    the multivariate soil chemical composition, controlling for island identity
    as a blocking factor (restricted permutation).

    Step 1 — Permutational MANOVA (PERMANOVA):
        Equivalent to adonis2(dist ~ Profundidade, strata = Ponto, permutations = 9999)
        in R/vegan. Labels are permuted ONLY within each island to preserve the
        hierarchical structure of the design (replicates within depth within island).

    Step 2 — Univariate paired t-tests and Wilcoxon signed-rank tests per variable
        (n = 6 pairs; one mean per island per depth), with Benjamini-Hochberg FDR
        correction for multiple comparisons.

Input:
    - solo_quimica_clean.csv   (output of Script 01)

Output:
    - resultado_permanova.txt
    - resultado_teste_profundidade.csv

References:
    Anderson (2001) Austral Ecol. 26:32–46.
    Benjamini & Hochberg (1995) J. R. Stat. Soc. B 57:289–300.
    Oksanen et al. (2022) vegan R package v2.6-4.
================================================================================
"""

import pandas as pd
import numpy as np
from scipy import stats

# ── Configuration ─────────────────────────────────────────────────────────────
RANDOM_SEED  = 42
N_PERM       = 9999
ALPHA        = 0.05

VARS_QUIM = [
    "pH", "Carbono", "MO", "N", "CN", "P", "Al", "Acidez",
    "Na", "K", "Ca", "Mg", "S", "CTC", "V", "Cu", "Zn", "Mn", "Fe"
]

# ── Load data ─────────────────────────────────────────────────────────────────
df = pd.read_csv("solo_quimica_clean.csv")
rng = np.random.default_rng(RANDOM_SEED)

# ── Helper: z-score standardisation ──────────────────────────────────────────
X = df[VARS_QUIM].to_numpy(dtype=float)
Z = (X - X.mean(axis=0)) / X.std(axis=0, ddof=1)

ilhas  = df["Ponto"].to_numpy()
prof   = df["Profundidade"].to_numpy()
ilhas_unicas = np.unique(ilhas)

# ── STEP 1: PERMANOVA ─────────────────────────────────────────────────────────
def ss_among_within(Zmat, grupos):
    """
    Compute SS_total and SS_within using the distance-based partitioning
    of Anderson (2001): SS = (1/n) * sum_{i<j} d_ij^2.
    SS_among = SS_total - SS_within.
    """
    N  = Zmat.shape[0]
    D2 = np.sum((Zmat[:, None, :] - Zmat[None, :, :]) ** 2, axis=-1)
    iu = np.triu_indices(N, k=1)
    ss_total = D2[iu].sum() / N
    ss_within = 0.0
    for g in np.unique(grupos):
        idx = np.where(grupos == g)[0]
        if len(idx) < 2:
            continue
        sub    = D2[np.ix_(idx, idx)]
        iu_sub = np.triu_indices(len(idx), k=1)
        ss_within += sub[iu_sub].sum() / len(idx)
    return ss_total, ss_within


def centre_by_island(Zmat, ilhas_vec):
    """Remove island means to isolate within-island variation (blocking)."""
    Zc = Zmat.copy()
    for g in np.unique(ilhas_vec):
        idx       = np.where(ilhas_vec == g)[0]
        Zc[idx,:] = Zmat[idx,:] - Zmat[idx,:].mean(axis=0)
    return Zc


# Centre by island (remove block effect before testing depth)
Zc = centre_by_island(Z, ilhas)
ss_total, ss_within = ss_among_within(Zc, prof)
ss_among = ss_total - ss_within
df1, df2 = 1, 29   # 35 total - 5 (island) - 1 (depth) = 29 residual

F_obs = (ss_among / df1) / (ss_within / df2)
R2    = ss_among / ss_total

# Restricted permutation: shuffle depth labels ONLY within each island
F_perm = np.empty(N_PERM)
for p in range(N_PERM):
    prof_perm = prof.copy()
    for g in ilhas_unicas:
        idx          = np.where(ilhas == g)[0]
        prof_perm[idx] = rng.permutation(prof[idx])
    _, ss_w_p = ss_among_within(Zc, prof_perm)
    ss_a_p    = ss_total - ss_w_p
    F_perm[p] = (ss_a_p / df1) / (ss_w_p / df2)

p_perm = (np.sum(F_perm >= F_obs) + 1) / (N_PERM + 1)

# Variance explained by island (context)
ss_total_raw, ss_w_ilha = ss_among_within(Z, ilhas)
pct_ilha = 100 * (ss_total_raw - ss_w_ilha) / ss_total_raw
pct_prof  = 100 * ss_among / ss_total

print("=" * 65)
print("PERMANOVA — Effect of depth (stratified by island)")
print("=" * 65)
print(f"  Pseudo-F ({df1}, {df2}) = {F_obs:.3f}")
print(f"  p (permutation, n = {N_PERM}) = {p_perm:.4f}")
print(f"  R² (depth, within island) = {R2:.3f}")
print(f"  Variance explained by island (block): {pct_ilha:.1f}%")
print(f"  Within-island variance explained by depth: {pct_prof:.1f}%")

# Save permanova result
with open("resultado_permanova.txt", "w") as f:
    f.write("PERMANOVA — Effect of sampling depth on soil chemical composition\n")
    f.write("(restricted permutation within island; equivalent to adonis2 strata)\n\n")
    f.write(f"Pseudo-F ({df1}, {df2}) = {F_obs:.4f}\n")
    f.write(f"p-value (permutation, n = {N_PERM}) = {p_perm:.4f}\n")
    f.write(f"R² (depth) = {R2:.4f}\n")
    f.write(f"Variance explained by island (block): {pct_ilha:.2f}%\n")
    f.write(f"Variance explained by depth (within island): {pct_prof:.2f}%\n")
    f.write(f"\nRandom seed: {RANDOM_SEED}\n")

# ── STEP 2: Univariate paired tests per variable ──────────────────────────────
agg = (
    df.groupby(["Ponto", "Profundidade"])[VARS_QUIM]
    .mean()
    .reset_index()
)

resultados = []
for v in VARS_QUIM:
    pivot = agg.pivot(index="Ponto", columns="Profundidade", values=v)
    m010  = pivot["0-10"].to_numpy()
    m1020 = pivot["10-20"].to_numpy()

    t_stat, p_t = stats.ttest_rel(m010, m1020)
    try:
        _, p_w = stats.wilcoxon(m010, m1020)
    except ValueError:
        p_w = np.nan

    resultados.append({
        "variable":     v,
        "mean_0_10":    round(m010.mean(), 3),
        "mean_10_20":   round(m1020.mean(), 3),
        "delta_pct":    round(100 * (m1020.mean() - m010.mean()) / m010.mean(), 1),
        "t_stat":       round(t_stat, 3),
        "p_t":          round(p_t, 4),
        "p_wilcoxon":   round(p_w, 4) if not np.isnan(p_w) else np.nan,
    })

res = pd.DataFrame(resultados).sort_values("p_t")

# Benjamini-Hochberg FDR correction
pvals  = res["p_t"].to_numpy()
order  = np.argsort(pvals)
ranked = pvals[order]
m      = len(pvals)
fdr_sorted = ranked * m / (np.arange(m) + 1)
fdr_sorted = np.minimum.accumulate(fdr_sorted[::-1])[::-1]
fdr    = np.empty(m)
fdr[order] = np.minimum(fdr_sorted, 1.0)
res["p_t_FDR"] = fdr.round(4)

print("\nPaired tests per variable (n = 6 island pairs):")
pd.set_option("display.float_format", lambda x: f"{x:.4f}")
print(res.to_string(index=False))

res.to_csv("resultado_teste_profundidade.csv", index=False)
print("\nFile saved: resultado_teste_profundidade.csv")
print("File saved: resultado_permanova.txt")
