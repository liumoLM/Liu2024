# Build a wide-format CSV joining sample metadata (Table S17) to the
# per-sample signature assignment matrix (Table S10).
#
# Steps:
#   1. Read Table S17 metadata (.xlsx) and Table S10 assignments (.tsv).
#   2. Force both column names (sample IDs) and row names (signature names)
#      to syntactically valid R names with make.names(unique = TRUE), so that
#      the same transformation can be applied to the metadata Patient column
#      and the assignment column names for a clean join.
#   3. Transpose the assignment matrix so samples are rows and signatures are
#      columns.
#   4. Join to the metadata by the (transformed) sample / Patient identifier.
#   5. Write a CSV with columns: cancer.type, sample, <signature1>,
#      <signature2>, ... (one row per sample).
#
# Output (written to the same folder):
#   Table_S10_S17_joined_wide.csv

library(readxl)

this_dir <- local({
  f <- NULL
  try(f <- sys.frame(1)$ofile, silent = TRUE)
  if (is.null(f)) {
    a <- commandArgs(trailingOnly = FALSE)
    m <- a[grep("^--file=", a)]
    if (length(m) > 0) f <- sub("^--file=", "", m)
  }
  if (is.null(f)) return(getwd())
  f <- gsub("~\\+~", " ", f)
  dirname(normalizePath(f))
})

sup_dir <- file.path(this_dir, "..", "Sup Tables")
assignment_file <- file.path(
  sup_dir,
  "Table S10 83-type and 89-type signature assignment.tsv"
)
metadata_file <- file.path(
  sup_dir,
  "Table S17 metadata of 6975 samples.xlsx"
)
output_file <- file.path(this_dir, "Table_S10_S17_joined_wide.csv")

message("Loading metadata: ", metadata_file)
meta <- as.data.frame(read_excel(metadata_file, sheet = "Sheet1"))
stopifnot(all(c("Patient", "Cancer Type") %in% colnames(meta)))

message("Loading assignment matrix: ", assignment_file)
assign_mat <- read.delim(assignment_file, row.names = 1, check.names = FALSE)

# Force valid R names on sample column names and signature row names.
colnames(assign_mat) <- make.names(colnames(assign_mat), unique = TRUE)
rownames(assign_mat) <- make.names(rownames(assign_mat), unique = TRUE)

# Transpose: samples become rows, signatures become columns.
t_assign <- as.data.frame(t(assign_mat), check.names = FALSE)
t_assign$sample <- rownames(t_assign)
rownames(t_assign) <- NULL

# Apply the same name transformation to the metadata Patient column so the
# join key is consistent with the transformed assignment column names.
meta$sample <- make.names(meta$Patient, unique = TRUE)
meta$cancer.type <- meta[["Cancer Type"]]

joined <- merge(
  meta[, c("sample", "cancer.type")],
  t_assign,
  by = "sample",
  all = FALSE
)

# Order columns: cancer.type, sample, then signatures in the original
# assignment-matrix order.
sig_cols <- rownames(assign_mat)
joined <- joined[, c("cancer.type", "sample", sig_cols)]

message(
  "Joined ",
  nrow(joined),
  " samples across ",
  length(sig_cols),
  " signatures"
)

write.csv(joined, output_file, row.names = FALSE)
message("Wrote: ", output_file)
