# UBE2H analysis: data download and placement

This repository does not include the large public input datasets. Download the required files from NCBI GEO and place them under the project-level `data/` directory exactly as shown below.

All analysis scripts use `here::here()`, so open the project through `UBE2H.Rproj` (or run R from anywhere inside the project) before executing the notebooks.

## Expected directory structure

```text
UBE2H/
├── code/
├── data/
│   ├── arm_coords.csv
│   ├── GSE111672/
│   ├── GSE84465/
│   ├── GSE188644/
│   └── GSE313191/
├── figures/
├── R/
├── results/
├── renv/
├── renv.lock
└── UBE2H.Rproj
```

Create the dataset directories from the project root if they are absent:

```bash
mkdir -p data/GSE111672 data/GSE84465 data/GSE188644 data/GSE313191
```

---

## GSE111672: PDAC single-cell RNA-seq

1. Open the NCBI GEO record for accession **GSE111672**.
2. Download the PDAC-A inDrop filtered expression matrix from the supplementary files.
3. Place it at:

```text
data/GSE111672/GSE111672_PDAC-A-indrop-filtered-expMat.txt
```

The downstream PDAC-A analysis also loads a generated CopyKAT object from:

```text
data/GSE111672/GSE111672_PDAC_A_inDrop_copykat.rds
```

This `.rds` file is not a GEO download. It is included with submission. Optionally generate it by running the PDAC-A CopyKAT script included in `code/`. The script saves the object directly to the path above.

Expected folder contents:

```text
data/GSE111672/
├── GSE111672_PDAC-A-indrop-filtered-expMat.txt
└── GSE111672_PDAC_A_inDrop_copykat.rds
```

---

## GSE84465: four-patient GBM single-cell RNA-seq

1. Open the NCBI GEO record for accession **GSE84465**.
2. Download the all-cell expression matrix and the GEO series-matrix metadata file.
3. Place the files at:

```text
data/GSE84465/GSE84465_GBM_All_data.csv.gz
data/GSE84465/GSE84465_series_matrix.txt
```

The downstream analysis loads a generated combined CopyKAT object from:

```text
data/GSE84465/GSE84465_all_CopyKAT_results.rds
```

This `.rds` file is not a GEO download. It contains the CopyKAT results for `BT_S1`, `BT_S2`, `BT_S4`, and `BT_S6`.Optionally generate it by running the GSE84465 notebook with:

```yaml
params:
  run_copykat: true
```

Else, use the default:

```yaml
params:
  run_copykat: false
```

Expected folder contents:

```text
data/GSE84465/
├── GSE84465_GBM_All_data.csv.gz
├── GSE84465_series_matrix.txt
└── GSE84465_all_CopyKAT_results.rds
```

---

## GSE188644: RPE-1 RNA-seq

1. Open the NCBI GEO record for accession **GSE188644**.
2. Download the processed supplementary data files used by the GSE188644 analysis script.
3. Place the downloaded files in:

```text
data/GSE188644/
```

Preserve the original GEO filenames. The analysis script resolves its inputs relative to this directory using `here::here("data", "GSE188644", ...)`.

Expected placement:

```text
data/GSE188644/
└── <GSE188644 supplementary files downloaded from GEO>
```

Do not place these files directly in `data/` or rename them unless the corresponding filenames in the analysis script are updated.

---

## GSE313191: ETiX embryoid single-cell RNA-seq

1. Open the NCBI GEO record for accession **GSE313191**.
2. Download the processed supplementary single-cell data used by the GSE313191 analysis script.
3. Place the downloaded files in:

```text
data/GSE313191/
```

Preserve the original GEO filenames. A processed Seurat object used in this project is named:

```text
data/GSE313191/GSE313191_seurat_list_9samples.rds
```

A downstream quality-controlled object may be generated as:

```text
data/GSE313191/GSE313191_seurat_list_DD_AA_QCd.rds
```

The quality-controlled `.rds` is a generated intermediate rather than a raw GEO download; regenerate it with the corresponding preprocessing script when necessary.

Expected folder contents:

```text
data/GSE313191/
├── GSE313191_seurat_list_9samples.rds
└── GSE313191_seurat_list_DD_AA_QCd.rds   # generated intermediate, when required
```

---

## Shared chromosome-arm coordinates

The CopyKAT arm-level scoring notebooks also require:

```text
data/arm_coords.csv
```

This small project input should remain included with the code archive.

---

## Restoring the R environment

From the project root, open `UBE2H.Rproj` and run:

```r
renv::restore()
```

Then verify the required input paths before running an analysis:

```r
file.exists(here::here(
  "data",
  "GSE111672",
  "GSE111672_PDAC-A-indrop-filtered-expMat.txt"
))

file.exists(here::here(
  "data",
  "GSE84465",
  "GSE84465_GBM_All_data.csv.gz"
))
```

A script will stop with a missing-file message when an expected input has not been placed in the correct dataset directory.
