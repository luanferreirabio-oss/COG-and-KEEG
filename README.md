# Floodplain Cacao — Soil Chemistry and Metagenomics (Mocajuba, Pará)

Reproducible analysis package for the manuscript:

> **Integration of metagenomic profiles and soil chemical attributes in the determination
> of floodplain cacao productivity in Mocajuba, Pará**

Six cacao-cultivating *várzea* islands in the Lower Tocantins region, Pará, Brazil.
Shotgun soil metagenomics (KEGG and eggNOG/COG) integrated with 19 soil chemical
attributes and annual cacao production.

---

## Key results

| Result | Value |
|---|---|
| PERMANOVA — Island | Pseudo-F(5,5) = 3.47; R² = 0.651; p = 0.0062 |
| PERMANOVA — Depth | Pseudo-F(1,5) = 4.30; R² = 0.161; p = 0.060 |
| PERMDISP — Depth | F(1,10) = 0.012; p = 0.903 (homogeneous) |
| Ca × production | ρ = 0.943; exact two-sided p = 0.0167; q(BH) = 0.316 |
| Mg × production | ρ = 0.886; exact two-sided p = 0.0333; q(BH) = 0.316 |
| KEGG orthologues, p ≤ 0.05 | 66 observed vs. 131 expected by chance; **q < 0.05: 0** |
| COG orthologous groups, p ≤ 0.05 | 60 observed vs. 92 expected by chance; **q < 0.05: 0** |
| COG functional categories | all \|ρ\| ≤ 0.714; all q ≥ 0.62 |
| Read mapping rate | 14.54–15.53% (mean 14.83%, CV 2.51%) |
| Depth confounding of richness | total RPK vs. eggNOG identifiers: ρ = 0.943; p = 0.017 |

---

## Experimental unit

Each island × depth combination was characterised by **three analytical determinations**
from equidistant sampling points. The median coefficient of variation among these
determinations is **0.25%**, which reflects analytical precision rather than field-scale
spatial heterogeneity. The island × depth mean was therefore adopted as the experimental
unit: **n = 12** (6 islands × 2 depths) for the PERMANOVA, and **n = 6** (island means)
for all correlations with production.

## Statistical inference at n = 6

The permutation space is finite (6! = 720) and is enumerated exhaustively. All p-values
are **exact two-sided** values. The minimum attainable p-value is **2/720 = 0.0028**;
across ~4,000 simultaneous tests the smallest possible Benjamini-Hochberg q-value exceeds
unity, so **no orthologue-level association can survive FDR correction at this sample
size**, regardless of effect size. Asymptotic approximations (`scipy.stats.spearmanr`)
are not used: they return p = 0 for |ρ| = 1 and propagate that value through any
multiple-testing correction.

---

## Repository structure

```
.
├── README.md
├── requirements.txt
├── LICENSE
├── data/
│   ├── solo_quimica.csv              # 36 rows: 6 islands × 2 depths × 3 determinations
│   ├── Join_KEGG.csv                 # 3,939 KEGG orthologues × 6 islands (CPM)
│   ├── Join_COG_eggNOG.csv           # 18,159 eggNOG identifiers × 6 islands (CPM)
│   ├── RAW_HUMANN_FILES.md       # how to obtain the raw HUMAnN tables
│   └── reference_COG/
│       ├── cog-24.def.tab            # NCBI COG-2024 definitions (retrieved 2026-08-25)
│       └── cog-24.fun.tab            # NCBI COG-2024 functional categories
├── scripts/
│   ├── 07_exact_permutation_and_COG_categories.py
│   ├── 08_corrections_permanova_permdisp_exact.py
│   ├── 09_mapping_rates.py
│   └── kegg_annotations_verified.py
├── results/
│   ├── tables/                       # all output tables (CSV)
│   └── figures/                      # Figures 2–6 (PNG, 200 dpi)
└── manuscript/
    ├── ARTIGO1_FINAL.docx
    └── ARTIGO1_FINAL.pdf
```

---

## Scripts

### `07_exact_permutation_and_COG_categories.py`
Exact permutation inference for KEGG and COG orthologues; Benjamini-Hochberg correction;
re-annotation from the authoritative NCBI COG-2024 release; aggregation into the 25 COG
functional categories.

**Outputs:** `KEGG_exact_correlations.csv`, `COG_exact_correlations.csv`,
`COG_functional_categories_CPM.csv`, `COG_functional_categories_stats.csv`

### `08_corrections_permanova_permdisp_exact.py`
Replicate-structure diagnostic; two-factor PERMANOVA over the twelve island × depth units;
PERMDISP; exact two-sided Spearman correlations for the 19 chemical attributes with BH
correction; depth sensitivity analysis; leave-one-island-out resampling; per-tree
productivity metric.

**Outputs:** `CORR_tabela1_produtividade.csv`, `CORR_tabela2_quimica_exata.csv`,
`CORR_tabela3_permanova.csv`, `CORR_permdisp.csv`,
`CORR_sensibilidade_profundidade.csv`, `CORR_leave_one_out.csv`,
`CORR_testes_profundidade.csv`

### `09_mapping_rates.py`
Extracts UNMAPPED and mapped RPK totals from the raw HUMAnN gene-family tables,
computes per-sample mapping rates, and tests whether detected functional richness
is confounded by sequencing depth. Requires the raw tables (see
`data/RAW_HUMANN_FILES.md`).

**Output:** `Supplementary_Table_S8_mapping_rates.csv`

---

## Running

```bash
git clone https://github.com/luandanbio/floodplain-cacao-mocajuba.git
cd floodplain-cacao-mocajuba
pip install -r requirements.txt

python scripts/07_exact_permutation_and_COG_categories.py
python scripts/08_corrections_permanova_permdisp_exact.py
python scripts/09_mapping_rates.py   # requires the raw HUMAnN tables
```

All random seeds are fixed at `42`. Permutation counts: 9,999 for PERMANOVA and PERMDISP;
complete enumeration (720) for all Spearman and Mantel tests.

---

## COG reference tables

Retrieved from the NCBI COG FTP site on **25 August 2026**:

```bash
wget https://ftp.ncbi.nih.gov/pub/COG/COG2024/data/cog-24.def.tab
wget https://ftp.ncbi.nih.gov/pub/COG/COG2024/data/cog-24.fun.tab
wget https://ftp.ncbi.nih.gov/pub/COG/COG2024/data/Readme.COG2024.txt
```

The COG database receives continuous updates within a release, so the access date is part
of the provenance of the result. Copies of both tables are deposited here under
`data/reference_COG/`.

**Citation:** Galperin, M. Y., Vera Alvarez, R., Karamycheva, S., Makarova, K. S.,
Wolf, Y. I., Landsman, D., Koonin, E. V. (2025). COG database update 2024.
*Nucleic Acids Research* 53(D1), D356–D363. https://doi.org/10.1093/nar/gkae983

COG data are produced by the NCBI / NLM / NIH and are in the public domain in the
United States.

---

## Data availability

Raw sequencing data: NCBI BioProject **PRJNA1224407** (BioSamples SRS24151884–SRS24151889,
corresponding to P1–P6).

---

## Note on the deposited environment specification

The decontamination environment specification deposited at
https://doi.org/10.5281/zenodo.17498295 lists Kraken2 and Bracken among its
dependencies. **These tools were not used in the analyses reported in this study.**
Read decontamination was performed with Bowtie2 against pre-indexed reference
databases, and taxonomic profiling was performed with MetaPhlAn as part of the
HUMAnN pipeline. The Kraken2 and Bracken entries are residual dependencies of the
environment template and produced no output used here.

## Software

| Package | Version |
|---|---|
| Python | 3.10 |
| pandas | 2.0.3 |
| numpy | 1.24.4 |
| scipy | 1.11.4 |
| scikit-learn | 1.3.2 |
| matplotlib | 3.8.2 |

Bioinformatic pipeline: FastQC v0.11.9, Fastp v0.23.2, Bowtie2 v2.4.5, HUMAnN3,
MetaPhlAn, DIAMOND.

---

## Archived releases

Research data are deposited in two Zenodo records:

| Record | Content | DOI |
|---|---|---|
| Decontamination | Reference indices, parameter files and per-sample read filtering outputs | https://doi.org/10.5281/zenodo.17498295 |
| Functional analysis | KEGG and COG abundance matrices, analysis scripts, derived tables, COG-2024 reference tables | https://doi.org/10.5281/zenodo.21345125 |

This repository mirrors the functional analysis record. Cite the corresponding DOI
according to which materials are being referenced.

## Citation

> Ferreira, L. D. S., da Silva, J. F. B. R., Vilhena, M. P. S. P., Alegria, O. C.,
> Ramos, R. T. J., de Sousa, M. P. A. Integration of metagenomic profiles and soil
> chemical attributes in the determination of floodplain cacao productivity in
> Mocajuba, Pará. [Journal, year — to be completed upon publication].

## Contact

Luan Daniel Silva Ferreira — luan.ferreirabio@gmail.com
Federal University of Pará (UFPA), Belém, PA, Brazil
ORCID: https://orcid.org/0000-0001-9187-6988

## License

MIT License — see `LICENSE`.
