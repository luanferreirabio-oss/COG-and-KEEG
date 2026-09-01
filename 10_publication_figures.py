"""
================================================================================
Script 10 - Publication-quality figure generation
================================================================================
Regenerates all manuscript figures to journal submission specification.

Specification
-------------
Resolution      600 dpi (combination art; meets Elsevier, Springer, Wiley,
                Frontiers and MDPI requirements)
Width           190 mm (double-column) or 90 mm (single-column)
Formats         TIFF (LZW compressed) for submission
                PDF  (vector) for typesetting
                PNG  (600 dpi) for preview and preprint
Fonts           Embedded as TrueType (Type 42); Arial/Helvetica family
Minimum type    7 pt at final printed size
Colour          RGB (converted to CMYK by publisher if required)

Author: Luan Daniel Silva Ferreira (ORCID 0000-0001-9187-6988)
Federal University of Para (UFPA), Belem, PA, Brazil
Repository: https://github.com/luandanbio/floodplain-cacao-mocajuba
Last updated: 1 September 2026
================================================================================
"""

import itertools
import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

# ── Publication specification ────────────────────────────────────────────────
DPI = 600
MM = 1 / 25.4
W2 = 190 * MM          # double-column width, inches
W1 = 90 * MM           # single-column width, inches

matplotlib.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["DejaVu Sans", "Arial", "Helvetica"],
    "font.size": 7,
    "axes.labelsize": 7.5,
    "axes.titlesize": 8,
    "xtick.labelsize": 6.5,
    "ytick.labelsize": 6.5,
    "legend.fontsize": 6.5,
    "figure.titlesize": 9,
    "axes.linewidth": 0.6,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "xtick.major.size": 2.5,
    "ytick.major.size": 2.5,
    "lines.linewidth": 1.0,
    "patch.linewidth": 0.5,
    "pdf.fonttype": 42,      # embed TrueType, not Type 3
    "ps.fonttype": 42,
    "svg.fonttype": "none",
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.02,
    "savefig.facecolor": "white",
    "savefig.edgecolor": "none",
    "figure.facecolor": "white",
})

OUT = Path("figures_publication")
OUT.mkdir(exist_ok=True)


def save(fig, stem, max_width_mm=190):
    """Write TIFF (LZW, RGB), PDF (vector) and PNG at publication resolution.

    TIFF output is flattened to RGB on white (journals reject alpha channels)
    and downscaled if the tight bounding box exceeds the column width.
    """
    from PIL import Image
    Image.MAX_IMAGE_PIXELS = None

    tmp = OUT / f"_{stem}_tmp.png"
    fig.savefig(tmp, dpi=DPI, format="png", facecolor="white")
    fig.savefig(OUT / f"{stem}.pdf", format="pdf", facecolor="white")
    plt.close(fig)

    im = Image.open(tmp)
    if im.mode in ("RGBA", "LA", "P"):
        bg = Image.new("RGB", im.size, "white")
        bg.paste(im, mask=im.convert("RGBA").split()[-1])
        im = bg
    else:
        im = im.convert("RGB")

    max_px = int(round(max_width_mm / 25.4 * DPI))
    if im.width > max_px:
        h = int(round(im.height * max_px / im.width))
        im = im.resize((max_px, h), Image.LANCZOS)

    im.save(OUT / f"{stem}.tiff", format="TIFF", compression="tiff_lzw",
            dpi=(DPI, DPI))
    im.save(OUT / f"{stem}.png", format="PNG", dpi=(DPI, DPI))
    tmp.unlink()

    sizes = {e: (OUT / f"{stem}.{e}").stat().st_size / 1e6 for e in ("tiff", "pdf", "png")}
    print(f"  {stem}: {im.width}x{im.height} px, "
          f"{im.width/DPI*25.4:.0f} mm, {im.mode} | "
          f"TIFF {sizes['tiff']:.1f} MB | PDF {sizes['pdf']:.2f} MB | PNG {sizes['png']:.1f} MB")


# ── Data ─────────────────────────────────────────────────────────────────────
VARS = ["pH", "Carbono", "MO", "N", "CN", "P", "Al", "Acidez", "Na", "K",
        "Ca", "Mg", "S", "CTC", "V", "Cu", "Zn", "Mn", "Fe"]
SAMPLES = ["P1", "P2", "P3", "P4", "P5", "P6"]
NAMES = ["Santana", "Santaninha", "Angapijo", "Conceicao", "S. Joaquim", "Tauare"]
PROD = np.array([2000., 600., 450., 1035., 1000., 1500.])


def _find(name):
    for root in [Path("."), Path("data"), Path("..") / "data",
                 Path("/mnt/user-data/outputs")]:
        if (root / name).exists():
            return root / name
    raise FileNotFoundError(name)


raw = pd.read_csv(_find("solo_quimica.csv"))
df = raw[VARS].apply(pd.to_numeric, errors="coerce")
df["ilha"] = ["P" + str(x).replace("S", "") for x in raw["Ponto"]]
df["prof"] = [str(x) for x in raw["Profundidade"]]

med = df.groupby("ilha")[VARS].mean().loc[SAMPLES]
m010 = df[df.prof == "0-10"].groupby("ilha")[VARS].mean().loc[SAMPLES]
m1020 = df[df.prof == "10-20"].groupby("ilha")[VARS].mean().loc[SAMPLES]
s010 = df[df.prof == "0-10"].groupby("ilha")[VARS].std(ddof=1).loc[SAMPLES]
s1020 = df[df.prof == "10-20"].groupby("ilha")[VARS].std(ddof=1).loc[SAMPLES]

kegg = pd.read_csv(_find("Join_KEGG_corrigido.csv"), index_col=0)
egg = pd.read_csv(_find("Join_COG_corrigido.csv"), index_col=0)
cog = egg[[str(i).startswith("COG") for i in egg.index]]

CMAP = plt.cm.RdYlGn
NORM = (PROD - PROD.min()) / (PROD.max() - PROD.min())
COL = [CMAP(v) for v in NORM]

KEGG_LABELS = {
    "K01999": "livK; branched-chain aa transport, substrate-binding",
    "K01997": "livH; branched-chain aa transport, permease",
    "K01996": "livF; branched-chain aa transport, ATP-binding",
    "K02050": "ABC.SN.P; NitT/TauT family transport, permease",
    "K02049": "ABC.SN.A; NitT/TauT family transport, ATP-binding",
    "K01990": "ABC-2.A; ABC-2 type transport, ATP-binding",
    "K03088": "RNA polymerase sigma-70 factor, ECF subfamily",
    "K03704": "cspA; cold shock protein",
    "K00549": "metE; methionine synthase (cobalamin-independent)",
    "K02950": "rpsL; small subunit ribosomal protein S12",
}


def bray_curtis(d):
    a = d[SAMPLES].values.T
    D = np.zeros((6, 6))
    for i in range(6):
        for j in range(i + 1, 6):
            D[i, j] = D[j, i] = np.abs(a[i] - a[j]).sum() / (a[i].sum() + a[j].sum())
    return D


def pcoa(D):
    H = np.eye(6) - np.ones((6, 6)) / 6
    B = -0.5 * H @ (D ** 2) @ H
    ev, evec = np.linalg.eigh(B)
    o = np.argsort(ev)[::-1]
    ev, evec = ev[o], evec[:, o]
    pos = ev > 1e-10
    return evec[:, pos] * np.sqrt(ev[pos]), ev[pos] / ev[pos].sum()


def truncate(s, n):
    s = str(s)
    return s[:n] + "\u2026" if len(s) > n else s


print(f"Generating figures at {DPI} dpi, {190} mm width\n")

# ══ FIGURE 2 — Soil chemistry by island and depth ════════════════════════════
VB = ["pH", "Carbono", "MO", "N", "P", "Ca", "Mg", "K", "Al", "CTC", "V", "Mn", "Fe"]
UN = {"pH": "", "Carbono": "g kg$^{-1}$", "MO": "g kg$^{-1}$", "N": "g kg$^{-1}$",
      "P": "mg kg$^{-1}$", "Ca": "cmol$_c$ kg$^{-1}$", "Mg": "cmol$_c$ kg$^{-1}$",
      "K": "cmol$_c$ kg$^{-1}$", "Al": "cmol$_c$ kg$^{-1}$", "CTC": "cmol$_c$ kg$^{-1}$",
      "V": "%", "Mn": "mg kg$^{-1}$", "Fe": "mg kg$^{-1}$"}

x = np.arange(6)
w = 0.36
nrow = int(np.ceil(len(VB) / 3))
fig, axes = plt.subplots(nrow, 3, figsize=(W2, nrow * 1.05))
axes = axes.flatten()
for i, v in enumerate(VB):
    a = axes[i]
    a.bar(x - w / 2, m010[v], w, yerr=s010[v], label="0\u201310 cm", color="#2a9d8f",
          edgecolor="white", linewidth=0.4, capsize=1.5,
          error_kw=dict(lw=0.5, capthick=0.5, ecolor="#333"))
    a.bar(x + w / 2, m1020[v], w, yerr=s1020[v], label="10\u201320 cm", color="#e9c46a",
          edgecolor="white", linewidth=0.4, capsize=1.5,
          error_kw=dict(lw=0.5, capthick=0.5, ecolor="#333"))
    a.set_xticks(x)
    a.set_xticklabels(SAMPLES)
    u = UN.get(v, "")
    lbl = {"Carbono": "Organic C", "MO": "OM", "Acidez": "Pot. acidity",
           "CTC": "CEC", "V": "V", "CN": "C/N ratio"}.get(v, v)
    a.set_ylabel(f"{lbl} ({u})" if u else lbl)
    a.spines[["top", "right"]].set_visible(False)
    if i == 0:
        a.legend(frameon=False, loc="lower left", bbox_to_anchor=(0.0, 1.02),
                 ncol=2, handlelength=1.0, borderpad=0.2, labelspacing=0.2,
                 columnspacing=1.0)
for i in range(len(VB), len(axes)):
    axes[i].set_visible(False)
fig.tight_layout(pad=0.6, w_pad=1.2, h_pad=0.8)
save(fig, "Figure_2")

# ══ FIGURE 3 — PCA of soil chemistry ═════════════════════════════════════════
Z = StandardScaler().fit_transform(med[VARS].values)
pca = PCA()
sc = pca.fit_transform(Z)
ve = pca.explained_variance_ratio_
pc1, pc2 = sc[:, 0], sc[:, 1]
arrows = {v: (stats.pearsonr(Z[:, i], pc1)[0], stats.pearsonr(Z[:, i], pc2)[0])
          for i, v in enumerate(VARS)}
rp1, pp1 = stats.pearsonr(PROD, pc1)
rp2, pp2 = stats.pearsonr(PROD, pc2)

fig, a = plt.subplots(figsize=(W2 * 0.62, W2 * 0.56))
S = 2.0
for i in range(6):
    a.scatter(pc1[i], pc2[i], c=[COL[i]], s=48, zorder=5,
              edgecolors="k", linewidths=0.5)
ymid = np.median(pc2)
yr, xr = abs(pc2.max() - pc2.min()), abs(pc1.max() - pc1.min())
for i, p in enumerate(SAMPLES):
    a.annotate(f"{p} ({NAMES[i]})\n{int(PROD[i])} kg yr$^{{-1}}$",
               xy=(pc1[i], pc2[i]),
               xytext=(pc1[i] + xr * 0.08, pc2[i] + yr * (0.16 if pc2[i] >= ymid else -0.22)),
               fontsize=5.8, fontweight="bold", color="#1a1a2e",
               arrowprops=dict(arrowstyle="-", color="#888", lw=0.4),
               bbox=dict(fc="white", ec="none", alpha=0.75, pad=0.8), zorder=8)
for v in ["Ca", "Mg", "P", "N", "MO", "pH", "Al", "Fe", "K", "V", "Mn", "Cu"]:
    dx, dy = arrows[v]
    a.annotate("", xy=(dx * S, dy * S), xytext=(0, 0),
               arrowprops=dict(arrowstyle="->", color="#457b9d", lw=0.7, mutation_scale=6))
    a.text(dx * S * 1.08, dy * S * 1.08, v, fontsize=6, color="#457b9d",
           ha="left" if dx >= 0 else "right", fontweight="semibold")
a.annotate("", xy=(rp1 * S, rp2 * S), xytext=(0, 0),
           arrowprops=dict(arrowstyle="->", color="#e63946", lw=1.4, mutation_scale=8))
a.text(rp1 * S * 1.10, rp2 * S * 1.10, "Production", fontsize=6.5, color="#e63946",
       fontweight="bold", bbox=dict(fc="white", ec="#e63946", alpha=0.9, lw=0.5, pad=1.2))
a.axhline(0, color="grey", lw=0.4, ls="--", alpha=0.5)
a.axvline(0, color="grey", lw=0.4, ls="--", alpha=0.5)
a.set_xlabel(f"PC1 ({ve[0]*100:.1f}% of variance)")
a.set_ylabel(f"PC2 ({ve[1]*100:.1f}% of variance)")
a.spines[["top", "right"]].set_visible(False)
a.set_xlim(-3.9, 3.9)
a.set_ylim(-3.9, 3.9)
sm = plt.cm.ScalarMappable(cmap=CMAP, norm=plt.Normalize(PROD.min(), PROD.max()))
sm.set_array([])
cb = fig.colorbar(sm, ax=a, shrink=0.55, pad=0.02, aspect=18)
cb.set_label("Total production (kg yr$^{-1}$)", fontsize=6.5)
cb.ax.tick_params(labelsize=6)
fig.tight_layout(pad=0.4)
save(fig, "Figure_3")

# ══ FIGURE 4 — Heatmap of most abundant identifiers ══════════════════════════
fig = plt.figure(figsize=(W2 * 0.92, W2 * 0.84))
G = gridspec.GridSpec(2, 1, hspace=0.30, figure=fig)
for k, (nm, d, is_kegg) in enumerate([("KEGG", kegg, True), ("COG", cog, False)]):
    a = fig.add_subplot(G[k])
    top = d[SAMPLES].mean(axis=1).nlargest(10).index.tolist()
    sub = d.loc[top, SAMPLES].astype(float)
    z = sub.subtract(sub.mean(axis=1), axis=0).divide(
        sub.std(axis=1, ddof=1).replace(0, 1), axis=0)
    im = a.imshow(z.values, aspect="auto", cmap="RdYlGn", vmin=-2, vmax=2)
    a.set_xticks(range(6))
    a.set_xticklabels([f"{p}\n{int(v)} kg yr$^{{-1}}$" for p, v in zip(SAMPLES, PROD)])
    if is_kegg:
        labs = [f"{t} \u2014 {truncate(KEGG_LABELS.get(t, '\u2014'), 46)}" for t in top]
    else:
        ann = pd.read_csv(_find("COG_exact_correlations.csv")).set_index("ID") \
            if Path("/mnt/user-data/uploads/COG_exact_correlations.csv").exists() \
            else pd.DataFrame()
        labs = [f"{t} \u2014 {truncate(ann.loc[t,'COG_name'] if t in ann.index else '\u2014', 46)}"
                for t in top]
    a.set_yticks(range(10))
    a.set_yticklabels(labs, fontsize=5.8)
    a.set_title(("A. Ten most abundant KEGG orthologues" if is_kegg
                 else "B. Ten most abundant COG orthologous groups"),
                fontsize=7.5, fontweight="bold", loc="left", pad=5)
    for r in range(10):
        for c in range(6):
            val = z.values[r, c]
            a.text(c, r, f"{val:.2f}", ha="center", va="center", fontsize=5.2,
                   color="black" if abs(val) < 1.2 else "white")
    cb = fig.colorbar(im, ax=a, shrink=0.85, pad=0.015, aspect=16)
    cb.set_label("Abundance (z-score)", fontsize=6.5)
    cb.ax.tick_params(labelsize=6)
save(fig, "Figure_4")

# ══ FIGURE 5 — PCoA of functional profiles ═══════════════════════════════════
Dk, De = bray_curtis(kegg), bray_curtis(egg)
fig, axes = plt.subplots(1, 2, figsize=(W2, W2 * 0.44))
fig.subplots_adjust(wspace=0.30)
for a, (nm, D, panel) in zip(axes, [("KEGG", Dk, "A"), ("eggNOG", De, "B")]):
    C, ve_ = pcoa(D)
    c1, c2 = C[:, 0], C[:, 1]
    r1, q1 = stats.pearsonr(PROD, c1)
    r2, q2 = stats.pearsonr(PROD, c2)
    u = D[np.triu_indices(6, 1)]
    for i in range(6):
        a.scatter(c1[i], c2[i], c=[COL[i]], s=48, zorder=5, edgecolors="k", linewidths=0.5)
    ymid = np.median(c2)
    yr, xr = abs(c2.max() - c2.min()), abs(c1.max() - c1.min())
    for i, p in enumerate(SAMPLES):
        a.annotate(f"{p}\n{int(PROD[i])}", xy=(c1[i], c2[i]),
                   xytext=(c1[i] + xr * 0.07, c2[i] + yr * (0.14 if c2[i] >= ymid else -0.20)),
                   fontsize=5.8, fontweight="bold",
                   bbox=dict(fc="white", ec="none", alpha=0.75, pad=0.6), zorder=8)
    a.axhline(0, color="grey", lw=0.4, ls="--", alpha=0.5)
    a.axvline(0, color="grey", lw=0.4, ls="--", alpha=0.5)
    a.set_xlim(c1.min() - xr * 0.42, c1.max() + xr * 0.42)
    a.set_ylim(c2.min() - yr * 0.45, c2.max() + yr * 0.45)
    a.set_xlabel(f"PCo1 ({ve_[0]*100:.1f}% of variance)")
    a.set_ylabel(f"PCo2 ({ve_[1]*100:.1f}% of variance)")
    a.set_title(f"{panel}. {nm} (Bray-Curtis {u.min():.3f}\u2013{u.max():.3f})",
                fontsize=7.5, fontweight="bold", loc="left")
    a.spines[["top", "right"]].set_visible(False)
sm = plt.cm.ScalarMappable(cmap=CMAP, norm=plt.Normalize(PROD.min(), PROD.max()))
sm.set_array([])
cb = fig.colorbar(sm, ax=axes, shrink=0.6, pad=0.02, aspect=18)
cb.set_label("Total production (kg yr$^{-1}$)", fontsize=6.5)
cb.ax.tick_params(labelsize=6)
save(fig, "Figure_5")

# ══ FIGURE 6 — Distribution of correlation coefficients ══════════════════════
kx = pd.read_csv(_find("KEGG_exact_correlations.csv"))
cx = pd.read_csv(_find("COG_exact_correlations.csv"))
fig, axes = plt.subplots(1, 2, figsize=(W2, W2 * 0.38))
fig.subplots_adjust(wspace=0.22)
for a, (nm, d, panel) in zip(axes, [("KEGG orthologues", kx, "A"),
                                     ("COG orthologous groups", cx, "B")]):
    a.hist(d["rho"], bins=45, color="#457b9d", edgecolor="white", linewidth=0.3)
    a.axvspan(-0.886, 0.886, color="grey", alpha=0.16,
              label="$|\\rho|$ < 0.886 ($p$ > 0.05)")
    a.axvline(0, color="#333", lw=0.5, ls="--")
    ns = int((d["p_exact_two_sided"] <= 0.05).sum())
    ex = int(round(len(d) * 0.03333))
    nq = int((d["q_BH"] < 0.05).sum())
    a.set_xlabel("Spearman $\\rho$ (vs. total production)")
    a.set_ylabel("Number of identifiers")
    a.set_title(f"{panel}. {nm} ($n$ = {len(d):,})\n"
                f"$p \\leq$ 0.05: {ns} observed, {ex} expected; $q$ < 0.05: {nq}",
                fontsize=7.5, fontweight="bold", loc="left")
    a.legend(frameon=False, loc="upper right", handlelength=1.2)
    a.spines[["top", "right"]].set_visible(False)
fig.tight_layout(pad=0.4)
save(fig, "Figure_6")

print(f"\nAll figures written to {OUT}/ at {DPI} dpi")
