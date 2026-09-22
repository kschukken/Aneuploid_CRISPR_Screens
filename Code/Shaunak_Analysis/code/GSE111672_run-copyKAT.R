suppressPackageStartupMessages({
  library(here)
  library(data.table)
  library(copykat)
})

patient <- "PDAC-A"
n_cores <- 12L

counts_file <- here(
  "data",
  "GSE111672",
  "GSE111672_PDAC-A-indrop-filtered-expMat.txt"
)

copykat_rds <- here(
  "data",
  "GSE111672",
  "GSE111672_PDAC_A_inDrop_copykat.rds"
)

predictions_file <- here(
  "data",
  "GSE111672",
  "GSE111672_PDAC_A_inDrop_copykat_predictions.txt"
)

if (!file.exists(counts_file)) {
  stop("Counts file not found: ", counts_file)
}

message("Reading ", patient, " expression matrix...")

df <- data.table::fread(
  counts_file,
  sep = "\t",
  header = TRUE,
  data.table = FALSE,
  check.names = FALSE
)

# Repeated column headers are cell-type labels.
cell_type <- trimws(
  gsub("\r", "", colnames(df)[-1])
)

# Unique IDs within PDAC-A.
cell_ids <- paste0(
  "PDAC_A_cell_",
  seq_along(cell_type)
)

# Gene symbols are stored in the first column.
genes <- trimws(
  gsub("\r", "", as.character(df[[1]]))
)

keep_gene <- (
  !is.na(genes) &
    nzchar(genes) &
    !duplicated(genes)
)

mat <- as.matrix(
  df[keep_gene, -1, drop = FALSE]
)

storage.mode(mat) <- "integer"

rownames(mat) <- genes[keep_gene]
colnames(mat) <- cell_ids

rm(df)

# Remove RBC cells.
drop_rbc <- grepl(
  "^RBC",
  cell_type,
  ignore.case = TRUE
)

mat <- mat[, !drop_rbc, drop = FALSE]
cell_type <- cell_type[!drop_rbc]
cell_ids <- cell_ids[!drop_rbc]

stopifnot(
  identical(colnames(mat), cell_ids)
)

metadata <- data.frame(
  cell_name = cell_ids,
  cell_type = cell_type,
  patient = patient,
  row.names = cell_ids,
  stringsAsFactors = FALSE
)

cat(
  patient, ":",
  nrow(mat), "genes x",
  ncol(mat), "cells after removing RBCs\n"
)

print(
  sort(
    table(metadata$cell_type),
    decreasing = TRUE
  )
)

# Ductal and Cancer groups are observations.
# Other cell types are normal references.
obs_pattern <- "Ductal|Cancer"

is_observation <- grepl(
  obs_pattern,
  metadata$cell_type,
  ignore.case = TRUE
)

norm_cells <- rownames(metadata)[!is_observation]
obs_cells <- rownames(metadata)[is_observation]

cat("\n", patient, " reference groups:\n", sep = "")

print(
  sort(
    table(metadata[norm_cells, "cell_type"]),
    decreasing = TRUE
  )
)

cat("\n", patient, " observation groups:\n", sep = "")

print(
  sort(
    table(metadata[obs_cells, "cell_type"]),
    decreasing = TRUE
  )
)

stopifnot(
  length(norm_cells) >= 5,
  length(obs_cells) >= 1,
  all(norm_cells %in% colnames(mat))
)

message("\nRunning CopyKAT for ", patient, "...")

copykat_result <- copykat::copykat(
  rawmat = mat,
  id.type = "S",
  norm.cell.names = norm_cells,
  ngene.chr = 5,
  win.size = 25,
  KS.cut = 0.05,
  sam.name = "GSE111672_PDAC_A_inDrop",
  distance = "euclidean",
  genome = "hg20",
  n.cores = n_cores
)

pred <- copykat_result$prediction

stopifnot(
  all(pred$cell.names %in% rownames(metadata))
)

pred$patient <- patient
pred$cell_type <- metadata[
  pred$cell.names,
  "cell_type"
]

pred$input_group <- ifelse(
  pred$cell.names %in% norm_cells,
  "reference",
  "observation"
)

cat("\n", patient, " CopyKAT calls:\n", sep = "")

print(
  table(
    pred$copykat.pred,
    useNA = "ifany"
  )
)

cat("\n", patient, " calls by input group:\n", sep = "")

print(
  table(
    pred$input_group,
    pred$copykat.pred,
    useNA = "ifany"
  )
)

cat("\n", patient, " calls by cell type:\n", sep = "")

print(
  table(
    pred$cell_type,
    pred$copykat.pred,
    useNA = "ifany"
  )
)

saveRDS(
  copykat_result,
  copykat_rds
)

write.table(
  pred,
  predictions_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
