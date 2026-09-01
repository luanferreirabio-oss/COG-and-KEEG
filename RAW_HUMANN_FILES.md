# Raw HUMAnN gene-family tables

The six raw HUMAnN gene-family output files (one per island, ~40 MB each,
221 MB total) are **not included in this repository** because of GitHub file
size constraints. They are deposited at:

**https://doi.org/10.5281/zenodo.21345125**

## Files and island correspondence

Galaxy dataset numbers 28–33 correspond to islands P1–P6 in order.

| Island | Name | Galaxy file |
|---|---|---|
| P1 | Santana | `Galaxy181-_HUMAnN_on_dataset_28_and_113__Gene_families_and_their_abundance_.tabular` |
| P2 | Santaninha | `Galaxy190-_HUMAnN_on_dataset_29_and_157__Gene_families_and_their_abundance_.tabular` |
| P3 | Angapijó | `Galaxy316-_HUMAnN_on_dataset_30_and_161__Gene_families_and_their_abundance_.tabular` |
| P4 | Conceição | `Galaxy325-_HUMAnN_on_dataset_31_and_165__Gene_families_and_their_abundance_.tabular` |
| P5 | São Joaquim | `Galaxy222-_HUMAnN_on_dataset_32_and_169__Gene_families_and_their_abundance_.tabular` |
| P6 | Tauaré | `Galaxy231-_HUMAnN_on_dataset_33_and_173__Gene_families_and_their_abundance_.tabular` |

The P3 file carries the column header `humann-P3_Abundance-RPKs`, which
independently confirms the dataset-to-island correspondence.

## Format

Two columns, tab-separated: gene family identifier and abundance in RPK
(reads per kilobase). The first data row is `UNMAPPED`, giving the RPK total
for reads not assigned to any UniRef90 family. Rows containing `|` are
species-stratified duplicates of the corresponding unstratified row and must
be excluded to avoid double counting.

## Reproducing Supplementary Table S8

Download the six files into `data/raw_humann/` and run:

```bash
python scripts/09_mapping_rates.py
```

This regenerates `Supplementary_Table_S8_mapping_rates.csv` and prints the
sequencing-depth confounding tests reported in Section 3.4.
