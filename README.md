# Floodplain Cacao Metagenomics — Analysis Pipeline

[![Python](https://img.shields.io/badge/Python-3.10%2B-blue)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

## Overview

This repository contains the complete statistical analysis pipeline for the study:

> **Integration of metagenomic profiles, soil chemical attributes and socioenvironmental
> variables in the determination of floodplain cacao productivity in Mocajuba, Pará, Brazil**
>
> *[Author names] — [Institution] — [Year]*

The study integrates shotgun soil metagenomics (KEGG and COG functional annotations),
soil chemical attributes (19 variables; 6 islands × 2 depths × 3 replicates), and cacao
productivity data from six *várzea* (Amazonian floodplain) islands in Mocajuba, Pará, Brazil.

---

## Repository Structure

```
.
├── README.md
├── requirements.txt
├── data/                          # Input data (place your files here)
│   ├── solo_quimica.csv
│   ├── Galaxy336-_COG_P1_.tabular
│   ├── Galaxy337-_COG_P2_.tabular
│   ├── ...
│   └── Galaxy347-_KEGG_P6_.tabular
├── scripts/
│   ├── 01_soil_data_loading.py
│   ├── 02_permanova_depth_effect.py
│   ├── 03_spearman_pca_soil_chemistry.py
│   ├── 04_metagenomics_normalisation_qc.py
│   ├── 05_mantel_pcoa.py
│   └── 06_functional_correlations_figures.py
└── outputs/                       # Generated automatically
    ├── tables/
    └── figures/
```

---

## Scripts

| Script | Description | Key outputs |
|--------|-------------|-------------|
| `01_soil_data_loading.py` | Load and clean soil chemistry data; aggregate to island-level mean | `solo_quimica_clean.csv`, `tabela_mestre.csv` |
| `02_permanova_depth_effect.py` | PERMANOVA (restricted permutation by island) + paired t-tests per variable + FDR correction | `resultado_permanova.txt`, `resultado_teste_profundidade.csv` |
| `03_spearman_pca_soil_chemistry.py` | Spearman correlations (soil chemistry × productivity) + PCA biplot | `spearman_quimica_produtividade.csv`, `fig_PCA.png` |
| `04_metagenomics_normalisation_qc.py` | Normalise HUMAnN3 outputs to CPM; build joint matrices; duplicate QC check | `Join_KEGG_corrigido.csv`, `Join_COG_corrigido.csv`, `Top10_*.csv` |
| `05_mantel_pcoa.py` | Mantel tests (3 comparisons) + PCoA (Bray-Curtis) + pairwise dissimilarities | `mantel_results.csv`, `bray_curtis_pairwise.csv`, `fig_PCoA.png` |
| `06_functional_correlations_figures.py` | Spearman ρ (all KOs/COGs × productivity) + heatmap + barplot/dotplot + chemistry bar chart | `suplementar_KEGG_completo.csv`, `suplementar_COG_completo.csv`, `fig_*.png` |

---

## Installation

```bash
git clone https://github.com/[your-username]/floodplain-cacao-metagenomics.git
cd floodplain-cacao-metagenomics
pip install -r requirements.txt
```

---

## Running the Pipeline

Run scripts sequentially from the repository root. Each script reads outputs
from the previous step:

```bash
python scripts/01_soil_data_loading.py
python scripts/02_permanova_depth_effect.py
python scripts/03_spearman_pca_soil_chemistry.py
python scripts/04_metagenomics_normalisation_qc.py
python scripts/05_mantel_pcoa.py
python scripts/06_functional_correlations_figures.py
```

> **Note:** All input files must be placed in the `data/` directory before running.
> Scripts write outputs to the current working directory by default;
> adjust the paths at the top of each script as needed.

---

## Requirements

```
pandas>=2.0
numpy>=1.24
scipy>=1.10
scikit-learn>=1.2
matplotlib>=3.7
```

See `requirements.txt` for exact versions used in the study.

---

## Statistical Methods

| Analysis | Method | Implementation |
|----------|--------|----------------|
| Depth effect (multivariate) | PERMANOVA (restricted permutation) | Custom Python; equivalent to `adonis2(strata=)` in R/vegan |
| Depth effect (univariate) | Paired t-test + Wilcoxon; BH-FDR | `scipy.stats` |
| Soil chemistry × productivity | Spearman ρ | `scipy.stats.spearmanr` |
| PCA | z-score + PCA | `scikit-learn` |
| Functional dissimilarity | Bray-Curtis | Custom Python |
| Matrix congruence | Mantel test (Pearson r, 9,999 perm.) | Custom Python |
| Functional ordination | PCoA (double-centering) | Custom Python |
| Functional × productivity | Spearman ρ (all categories) | `scipy.stats.spearmanr` |

All random seeds are fixed at `42` for reproducibility.

---

## Data Availability

Raw metagenomic sequencing data are available at [NCBI SRA / ENA — accession to be added
upon acceptance]. Processed functional abundance tables are provided in this repository
under `data/`.

---

## Citation

[Citation to be added upon publication]

---

## License

This project is licensed under the MIT License — see the [LICENSE](LICENSE) file for details.
