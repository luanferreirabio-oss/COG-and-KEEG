"""
================================================================================
Script 06 — Functional Categories × Productivity: Spearman Correlations
            and Figures (Heatmap, Barplot, Dotplot)
================================================================================
Study: Integration of metagenomic profiles, soil chemical attributes and
       socioenvironmental variables in the determination of floodplain cacao
       productivity in Mocajuba, Pará, Brazil.

Author: [Your name]
Institution: [Your institution]
Date: 2024

Description:
    For each canonical KEGG KO and COG category with non-zero variance across
    the six islands, computes Spearman rank correlation (ρ) with productivity
    (n = 6). Generates:
        - Complete correlation tables (Supplementary Material)
        - Heatmap of Top-10 most abundant categories (Figure 4)
        - Barplot + dotplot of biologically annotated categories (Figure 6)
        - Bar chart of soil chemistry by island and depth (Figure 2)

Input:
    - Join_KEGG_corrigido.csv      (output of Script 04)
    - Join_COG_corrigido.csv       (output of Script 04)
    - solo_quimica_clean.csv       (output of Script 01)

Output:
    - suplementar_KEGG_completo.csv
    - suplementar_COG_completo.csv
    - fig_heatmap.png    (Figure 4)
    - fig_funcional.png  (Figure 6)
    - fig_barras.png     (Figure 2)
================================================================================
"""

import pandas as pd
import numpy as np
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ── Configuration ─────────────────────────────────────────────────────────────
SAMPLES      = ["P1","P2","P3","P4","P5","P6"]
ISLAND_NAMES = ["Santana","Santaninha","Angapijó","Conceição","S. Joaquim","Tauaré"]
PROD_VALS    = np.array([2000, 600, 450, 1035, 1000, 1500])
RHO_THRESH   = 0.886   # |ρ| threshold corresponding to p ≤ 0.019 for n = 6

# Short descriptions for annotated categories
KEGG_DESC = {
    "K00074": "BHBD (C cycling)",
    "K05587": "NifQ (N cycling)",
    "K05816": "FepA (Fe/siderop.)",
    "K11690": "Phos-ABC (P)",
    "K01432": "Allophanate hydrolase",
    "K02193": "CtaG (aerob. resp.)",
    "K13542": "CYP450 (org. ox.)",
    "K15371": "ABC-nutrient",
}
COG_DESC = {
    "COG2373": "Fe-S cluster biosyn.",
    "COG0201": "RpoH (stress resp.)",
    "COG3145": "Fe-receptor OM",
    "COG3745": "Phos-ABC (P)",
    "COG3288": "Glycosyl hydrolase",
    "COG1454": "Laccase (lignin)",
    "COG0367": "Inorg. pyrophosphatase",
    "COG1004": "Thioredoxin (redox)",
    "COG2199": "AraC (C metab.)",
}
KEGG_TOP10_LBL = {
    "K01999":"ABC aa-transporter","K02050":"ABC phosphonate",
    "K01990":"ABC export","K03088":"Sigma RpoS",
    "K00626":"Acetyl-CoA acetyltransf.","K02003":"ABC perm. subunit",
    "K01955":"CPS large subunit","K00059":"3-oxoacyl-ACP reductase",
    "K00817":"PyrAsp aminotransf.","K01920":"D-Ala-D-Ala ligase",
}
COG_TOP10_LBL = {
    "COG0745":"Predicted ATPase","COG1028":"Acyl-CoA dehydrogen.",
    "COG0596":"ABC transporter","COG0683":"ABC perm. protein",
    "COG1116":"ABC ATPase","COG0329":"Glutamine synth. II",
    "COG0531":"ABC aa-binding","COG0466":"RecA recombinase",
    "COG0200":"Ribosomal prot. L2","COG0092":"Ribosomal prot. S3",
}

# ── Load data ─────────────────────────────────────────────────────────────────
kegg = pd.read_csv("Join_KEGG_corrigido.csv", index_col=0)
cog  = pd.read_csv("Join_COG_corrigido.csv",  index_col=0)
solo = pd.read_csv("solo_quimica_clean.csv")

cmap_p   = plt.cm.RdYlGn
norm_p   = (PROD_VALS - PROD_VALS.min()) / (PROD_VALS.max() - PROD_VALS.min())
colors_p = [cmap_p(v) for v in norm_p]


# ── 1. Spearman correlations — all categories ─────────────────────────────────
def spearman_all(df: pd.DataFrame, prefix: str,
                 prod: np.ndarray, samples: list) -> pd.DataFrame:
    """Compute Spearman ρ between each category CPM and productivity."""
    rows = []
    id_col = "KO" if prefix == "K" else "COG"
    for cat_id in df.index:
        if not str(cat_id).startswith(prefix):
            continue
        vals = df.loc[cat_id, samples].values.astype(float)
        if vals.std() < 1e-10:
            continue
        rho, p = stats.spearmanr(vals, prod)
        row = {
            id_col:        cat_id,
            "rho":         round(rho, 4),
            "p":           round(p,   5),
            "sig":         "**" if p < 0.01 else ("*" if p < 0.05 else "ns"),
            "direction":   "Positive" if rho > 0 else "Negative",
        }
        row.update({f"CPM_{s}": round(df.loc[cat_id, s], 4) for s in samples})
        rows.append(row)
    return (
        pd.DataFrame(rows)
        .sort_values("rho", key=abs, ascending=False)
        .reset_index(drop=True)
    )


spk = spearman_all(kegg, "K",   PROD_VALS, SAMPLES)
spc = spearman_all(cog,  "COG", PROD_VALS, SAMPLES)

spk.to_csv("suplementar_KEGG_completo.csv", index=False)
spc.to_csv("suplementar_COG_completo.csv",  index=False)

for name, df_sp in [("KEGG", spk), ("COG", spc)]:
    sig_cnt = (df_sp["p"] < 0.05).sum()
    print(f"{name}: {len(df_sp)} categories | {sig_cnt} with p < 0.05 | "
          f"top 3: {df_sp.iloc[:3, 0].tolist()}")


# ── 2. Figure 4 — Heatmap Top-10 ─────────────────────────────────────────────
from sklearn.preprocessing import StandardScaler

fig, axes = plt.subplots(1, 2, figsize=(18, 7))
fig.subplots_adjust(wspace=0.55)

for ax, (nome, df, lbl_map) in zip(axes, [
    ("KEGG", kegg, KEGG_TOP10_LBL),
    ("COG",  cog,  COG_TOP10_LBL),
]):
    top10_ids = df[SAMPLES].mean(axis=1).nlargest(10).index.tolist()
    sub   = df.loc[top10_ids, SAMPLES].astype(float)
    means = sub.mean(axis=1); stds = sub.std(axis=1).replace(0, 1)
    zmat  = sub.subtract(means, axis=0).divide(stds, axis=0)

    im = ax.imshow(zmat.values, aspect="auto", cmap="RdYlGn",
                   vmin=-2, vmax=2, interpolation="nearest")
    ax.set_xticks(range(6))
    ax.set_xticklabels(
        [f"{s}\n({int(v)} kg yr⁻¹)" for s, v in zip(SAMPLES, PROD_VALS)],
        fontsize=9
    )
    ax.set_yticks(range(10))
    ax.set_yticklabels(
        [f"{t} — {lbl_map.get(t, '—')}" for t in top10_ids], fontsize=8
    )
    ax.set_title(
        f"Top 10 {nome} categories by mean abundance\n(row-standardised z-score)",
        fontsize=11, fontweight="bold", pad=10
    )
    for r in range(10):
        for c in range(6):
            val = zmat.values[r, c]
            tc  = "black" if abs(val) < 1.2 else "white"
            ax.text(c, r, f"{val:.2f}", ha="center", va="center",
                    fontsize=7.5, color=tc)
    plt.colorbar(im, ax=ax, shrink=0.65, pad=0.02, label="Abundance (z-score)")

plt.suptitle(
    "Functional profile of the six floodplain islands\n"
    "Top 10 most abundant KEGG and COG categories",
    fontsize=13, fontweight="bold", y=1.02
)
plt.savefig("fig_heatmap.png", dpi=220, bbox_inches="tight")
plt.close()
print("Saved: fig_heatmap.png")


# ── 3. Figure 6 — Barplot + Dotplot of annotated categories ──────────────────
C_POS = "#e63946"; C_NEG = "#457b9d"

fig, axes2 = plt.subplots(2, 2, figsize=(18, 14))
fig.subplots_adjust(hspace=0.55, wspace=0.42)

for col_idx, (nome, df_sp, desc_map, id_col) in enumerate([
    ("KEGG", spk, KEGG_DESC, "KO"),
    ("COG",  spc, COG_DESC,  "COG"),
]):
    # ── barplot (annotated categories)
    ax_b = axes2[0, col_idx]
    ann  = (
        df_sp[df_sp[id_col].isin(desc_map)]
        .copy()
        .assign(label=lambda d: d[id_col].map(desc_map))
        .sort_values("rho")
    )
    bc = [C_POS if r > 0 else C_NEG for r in ann["rho"]]
    ax_b.barh(range(len(ann)), ann["rho"], color=bc, edgecolor="white", height=0.65)
    for i, (_, row) in enumerate(ann.iterrows()):
        xp  = row["rho"] + (0.02 if row["rho"] >= 0 else -0.02)
        ha  = "left" if row["rho"] >= 0 else "right"
        ax_b.text(xp, i, f"ρ = {row['rho']:.2f}", va="center",
                  ha=ha, fontsize=8, color="#333")
    ax_b.set_yticks(range(len(ann)))
    ax_b.set_yticklabels(ann["label"], fontsize=8.5)
    ax_b.axvline(0, color="grey", lw=0.8, ls="--")
    ax_b.set_xlabel("Spearman's ρ", fontsize=10)
    ax_b.set_title(
        f"{nome} — Annotated functional categories\n"
        "(Spearman correlation with productivity)",
        fontsize=10, fontweight="bold"
    )
    ax_b.spines["top"].set_visible(False)
    ax_b.spines["right"].set_visible(False)
    ax_b.set_xlim(-1.25, 1.25)
    ax_b.grid(axis="x", ls=":", alpha=0.4)

    # ── dotplot (top 10 positive + top 10 negative)
    ax_d   = axes2[1, col_idx]
    top    = pd.concat([
        df_sp[df_sp["rho"] < 0].nsmallest(10, "rho"),
        df_sp[df_sp["rho"] > 0].nlargest(10, "rho"),
    ]).reset_index(drop=True)
    dc     = [C_POS if r > 0 else C_NEG for r in top["rho"]]
    logp   = -np.log10(top["p"].clip(lower=1e-6))
    sizes  = (logp * 40 + 30).clip(upper=300)

    ax_d.scatter(top["rho"], range(len(top)), c=dc, s=sizes,
                 edgecolors="k", lw=0.5, zorder=4)
    ax_d.set_yticks(range(len(top)))
    ax_d.set_yticklabels(top[id_col], fontsize=8)
    ax_d.axvline(0, color="grey", lw=0.8, ls="--")
    ax_d.axhline(9.5, color="grey", lw=0.6, ls=":", alpha=0.6)
    ax_d.set_xlabel("Spearman's ρ", fontsize=10)
    ax_d.set_title(
        f"{nome} — Top 10 positive & negative\n(dot size ∝ −log₁₀(p))",
        fontsize=10, fontweight="bold"
    )
    ax_d.spines["top"].set_visible(False)
    ax_d.spines["right"].set_visible(False)
    ax_d.set_xlim(-1.25, 1.25)
    ax_d.grid(axis="x", ls=":", alpha=0.4)
    for lp_, ls_ in [(0.05, 50), (0.01, 130), (0.001, 210)]:
        ax_d.scatter([], [], c="grey", s=ls_, label=f"p = {lp_}",
                     edgecolors="k", lw=0.4)
    ax_d.legend(title="p-value", fontsize=7.5, title_fontsize=8,
                loc="lower right", framealpha=0.8)

from matplotlib.patches import Patch
leg = [
    Patch(fc=C_POS, label="Positive (↑ in productive islands)"),
    Patch(fc=C_NEG, label="Negative (↓ in productive islands)"),
]
fig.legend(handles=leg, loc="upper center", bbox_to_anchor=(0.5, 1.01),
           ncol=2, fontsize=10, framealpha=0.9)
plt.suptitle(
    "Functional categories correlated with cacao productivity\n"
    "(six Amazonian floodplain islands)",
    fontsize=13, fontweight="bold", y=1.04
)
plt.savefig("fig_funcional.png", dpi=220, bbox_inches="tight")
plt.close()
print("Saved: fig_funcional.png")


# ── 4. Figure 2 — Soil chemistry bar chart ───────────────────────────────────
VARS_BAR = ["pH","Carbono","MO","N","P","Ca","Mg","K","Al","CTC","V","Mn","Fe"]
UNITS    = {
    "pH":"","Carbono":"g kg⁻¹","MO":"g kg⁻¹","N":"g kg⁻¹",
    "P":"mg kg⁻¹","Ca":"cmolc kg⁻¹","Mg":"cmolc kg⁻¹","K":"cmolc kg⁻¹",
    "Al":"cmolc kg⁻¹","CTC":"cmolc kg⁻¹","V":"%","Mn":"mg kg⁻¹","Fe":"mg kg⁻¹"
}
agg = (
    solo
    .groupby(["Ponto","Profundidade"])[VARS_BAR]
    .agg(["mean","std"])
    .reset_index()
)
agg.columns = (
    ["Ponto","Profundidade"] +
    [f"{v}_{s}" for v in VARS_BAR for s in ["mean","std"]]
)
x    = np.arange(6); w = 0.35
cols2 = {"0-10": "#2a9d8f", "10-20": "#e9c46a"}
nrows = int(np.ceil(len(VARS_BAR) / 3))
fig, axes3 = plt.subplots(nrows, 3, figsize=(17, nrows * 3.1))
axes3 = axes3.flatten()

for ai, v in enumerate(VARS_BAR):
    ax = axes3[ai]
    for prof, shift in zip(["0-10","10-20"], [-w/2, w/2]):
        sub = agg[agg["Profundidade"] == prof].set_index("Ponto")
        mn  = [sub.loc[p, f"{v}_mean"] if p in sub.index else 0 for p in SAMPLES]
        sd  = [sub.loc[p, f"{v}_std"]  if p in sub.index else 0 for p in SAMPLES]
        ax.bar(x + shift, mn, w, label=f"{prof} cm", color=cols2[prof],
               edgecolor="white", yerr=sd, capsize=3,
               error_kw=dict(lw=1, capthick=1, ecolor="#444"))
    ax.set_xticks(x)
    ax.set_xticklabels(SAMPLES, fontsize=9)
    u = UNITS.get(v, "")
    ax.set_ylabel(f"{v} ({u})" if u else v, fontsize=9)
    ax.set_title(v, fontsize=10, fontweight="bold")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(labelsize=8)
    if ai == 0:
        ax.legend(title="Depth (cm)", fontsize=8, title_fontsize=8, framealpha=0.8)

for ai in range(len(VARS_BAR), len(axes3)):
    axes3[ai].set_visible(False)

plt.suptitle(
    "Soil chemical attributes by island and depth\n"
    "(mean ± SD, n = 3 per depth)",
    fontsize=12, fontweight="bold", y=1.01
)
plt.tight_layout()
plt.savefig("fig_barras.png", dpi=220, bbox_inches="tight")
plt.close()
print("Saved: fig_barras.png")
print("\nAll outputs saved successfully.")
