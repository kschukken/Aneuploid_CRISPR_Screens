# Commands to retrieve data

```shell
mkdir -p data/{GSE188644,GSE313191,GSE111672,GSE84465}

## GSE84465 
cd data/GSE84465
curl -O "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE84nnn/GSE84465/suppl/GSE84465_GBM_All_data.csv.gz"
curl -O "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE84nnn/GSE84465/matrix/GSE84465_series_matrix.txt.gz"
gunzip -k GSE84465_series_matrix.txt.gz 
cd ../..

## GSE111672 
cd data/GSE111672
curl -O "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE111nnn/GSE111672/suppl/GSE111672_PDAC-A-indrop-filtered-expMat.txt.gz"
gunzip -k GSE111672_PDAC-A-indrop-filtered-expMat.txt.gz
cd ../..

## GSE188644 
cd data/GSE188644
curl -O "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE188nnn/GSE188644/suppl/GSE188644_RAW.tar"
tar -xvf GSE188644_RAW.tar
cd ../..

## GSE313191 
cd data/GSE313191
curl -O "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE313nnn/GSE313191/suppl/GSE313191_RAW.tar"
tar -xvf GSE313191_RAW.tar
cd ../..
```

## metadata for 313191:

````R
library(GEOquery)
library(here)
gse <- getGEO("GSE313191")
meta <- pData(gse[[1]])
write.csv(meta, here("data", "GSE313191", "GSE313191_metadata.csv"), row.names = FALSE)

# title matching check
samples <- c("ESC_ETiX_DD","ESC_ETiX_DA","ESC_ETiX_AA",
             "iGata4_ETiX_DD","iGata4_ETiX_DA","iGata4_ETiX_AA",
             "TSC_ETiX_DD","TSC_ETiX_DA","TSC_ETiX_AA")
stopifnot(length(setdiff(samples, meta$title)) == 0)
```
