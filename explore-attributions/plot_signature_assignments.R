# Generate per-signature assignment hamburger plots from the supplementary
# tables in this folder.
#
# Inputs (in the same folder as this script):
#   Table S10 83-type and 89-type signature assignment.tsv
#       Tab-separated matrix; first column is the signature name (e.g.
#       "ID_K/InsDel_K"), remaining columns are bare sample IDs with mutation
#       counts attributed to that signature.
#   Table S17 metadata of 6975 samples.xlsx
#       Sheet "Sheet1" with columns including "Patient", "Cancer Type", and
#       "MSIseq" ("MSI-H"/"Non-MSI-H").
#
# Output (written to the same folder):
#   Table_S10_hamburger_plots_mf<min_fraction>.pdf
#       One page per signature; each dot is a sample, red dashed line is the
#       per-cancer-type median, MSI-H samples are red triangles, others are
#       black circles.
#
# Usage (run from the repo root or from this folder):
#   Rscript "Sup Tables/plot_signature_assignments.R"
#   Rscript "Sup Tables/plot_signature_assignments.R" --min-fraction 0.05

library(ggplot2)
library(scales)
library(argparser)
library(readxl)

p <- arg_parser("Hamburger plots from Sup Tables S10 + S17")
p <- add_argument(
  p,
  "--min-fraction",
  type = "double",
  default = 0.0,
  help = "Only plot a sample's point if the signature accounts for >= this fraction of the sample's total mutations"
)
args <- parse_args(p)

# Locate this script's directory so inputs/outputs resolve regardless of cwd
this_dir <- local({
  f <- NULL
  try(f <- sys.frame(1)$ofile, silent = TRUE)
  if (is.null(f)) {
    a <- commandArgs(trailingOnly = FALSE)
    m <- a[grep("^--file=", a)]
    if (length(m) > 0) f <- sub("^--file=", "", m)
  }
  if (is.null(f)) return(getwd())
  # Rscript replaces spaces in argv with "~+~"; undo that
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
output_file <- file.path(
  this_dir,
  sprintf("Table_S10_hamburger_plots_mf%s.pdf", args$min_fraction)
)

# Load sample metadata from S17
message("Loading metadata from: ", metadata_file)
sample_info <- as.data.frame(read_excel(metadata_file, sheet = "Sheet1"))
# Normalize column names we depend on
stopifnot(all(c("Patient", "Cancer Type", "MSIseq") %in% colnames(sample_info)))
sample_info$Cancer_Type <- sample_info[["Cancer Type"]]
sample_info$MSIseq_MSI.H <- sample_info[["MSIseq"]] == "MSI-H"

# Plotting parameters
params <- list(
  width = 12,
  height = 6,
  base_size = 12,
  point_size = 1.2,
  point_alpha = 0.6,
  min_samples = 1,
  genome_size_mb = NULL
)

plot_signature_by_cancer_type <- function(
  signature_values,
  signature_name = "Signature",
  genome_size_mb = NULL,
  min_samples = 1,
  log_scale = TRUE,
  exclude_zero = TRUE,
  order_by = c("median", "name", "count"),
  bg_colors = c("#EDF8B1", "#2C7FB8"),
  point_color = "black",
  median_color = "red",
  point_size = 1.5,
  point_alpha = 0.7,
  base_size = 14,
  sample_info = NULL
) {
  order_by <- match.arg(order_by)

  sample_names <- names(signature_values)
  cancer_types <- sub("::.*", "", sample_names)

  df_all <- data.frame(
    sample = sample_names,
    cancer_type = cancer_types,
    mutations = as.numeric(signature_values),
    stringsAsFactors = FALSE
  )

  if (!is.null(sample_info)) {
    df_all$patient_id <- sub("^.*::", "", df_all$sample)
    df_all$msi_status <- sample_info$MSIseq_MSI.H[
      match(df_all$patient_id, sample_info$Patient)
    ]
    df_all$dot_color <- ifelse(
      !is.na(df_all$msi_status) & df_all$msi_status == TRUE,
      "red",
      point_color
    )
    df_all$dot_shape <- ifelse(
      !is.na(df_all$msi_status) & df_all$msi_status == TRUE,
      17,
      16
    )
  } else {
    df_all$dot_color <- point_color
    df_all$dot_shape <- 16
  }

  total_by_type <- as.data.frame(table(df_all$cancer_type))
  colnames(total_by_type) <- c("cancer_type", "total_samples")

  nonzero_by_type <- as.data.frame(table(df_all$cancer_type[
    df_all$mutations > 0
  ]))
  colnames(nonzero_by_type) <- c("cancer_type", "nonzero_samples")

  sample_counts <- merge(
    total_by_type,
    nonzero_by_type,
    by = "cancer_type",
    all.x = TRUE
  )
  sample_counts$nonzero_samples[is.na(sample_counts$nonzero_samples)] <- 0

  df <- df_all

  if (!is.null(genome_size_mb)) {
    df$mutations <- df$mutations / genome_size_mb
    y_label <- "Number of Mutations per Megabase"
  } else {
    y_label <- "Number of Mutations"
  }

  if (exclude_zero) {
    df <- df[df$mutations > 0, ]
  }

  type_counts <- table(df$cancer_type)
  valid_types <- names(type_counts)[type_counts >= min_samples]
  df <- df[df$cancer_type %in% valid_types, ]

  if (nrow(df) == 0) {
    warning("No data remaining after filtering for ", signature_name)
    return(list(plot = NULL, sample_counts = sample_counts))
  }

  medians <- aggregate(mutations ~ cancer_type, data = df, FUN = median)
  colnames(medians) <- c("cancer_type", "median_mutations")

  counts <- as.data.frame(table(df$cancer_type))
  colnames(counts) <- c("cancer_type", "n_samples")

  summary_df <- merge(medians, counts, by = "cancer_type")

  if (order_by == "median") {
    summary_df <- summary_df[order(summary_df$median_mutations), ]
  } else if (order_by == "name") {
    summary_df <- summary_df[order(summary_df$cancer_type), ]
  } else if (order_by == "count") {
    summary_df <- summary_df[order(-summary_df$n_samples), ]
  }

  type_order <- summary_df$cancer_type
  df$cancer_type <- factor(df$cancer_type, levels = type_order)
  summary_df$cancer_type <- factor(summary_df$cancer_type, levels = type_order)

  df$x_pos <- as.numeric(df$cancer_type)
  summary_df$x_pos <- as.numeric(summary_df$cancer_type)

  df <- df[order(df$cancer_type, df$mutations), ]

  df$within_rank <- ave(
    seq_len(nrow(df)),
    df$cancer_type,
    FUN = function(x) seq_along(x)
  )
  df$group_size <- ave(
    seq_len(nrow(df)),
    df$cancer_type,
    FUN = length
  )

  df$x_plot <- df$x_pos +
    0.8 * (df$within_rank - 1) / pmax(df$group_size - 1, 1) -
    0.4
  df$x_plot[df$group_size == 1] <- df$x_pos[df$group_size == 1]

  n_types <- length(type_order)
  bg_df <- data.frame(
    xmin = seq(0.5, n_types - 0.5, by = 1),
    xmax = seq(1.5, n_types + 0.5, by = 1),
    fill = rep(bg_colors, length.out = n_types)
  )

  counts_for_plot <- sample_counts[sample_counts$cancer_type %in% type_order, ]
  counts_for_plot$cancer_type <- factor(
    counts_for_plot$cancer_type,
    levels = type_order
  )
  counts_for_plot$x_pos <- as.numeric(counts_for_plot$cancer_type)
  y_min <- min(df$mutations)
  y_label_pos <- y_min / 10

  p <- ggplot2::ggplot(df, ggplot2::aes(x = x_plot, y = mutations)) +
    ggplot2::geom_rect(
      data = bg_df,
      ggplot2::aes(
        xmin = xmin,
        xmax = xmax,
        ymin = -Inf,
        ymax = Inf,
        fill = fill
      ),
      inherit.aes = FALSE,
      alpha = 0.2
    ) +
    ggplot2::scale_fill_identity() +
    ggplot2::geom_point(
      ggplot2::aes(color = dot_color, shape = dot_shape),
      size = point_size,
      alpha = point_alpha
    ) +
    ggplot2::scale_color_identity() +
    ggplot2::scale_shape_identity() +
    ggplot2::geom_segment(
      data = summary_df,
      ggplot2::aes(
        x = x_pos - 0.4,
        xend = x_pos + 0.4,
        y = median_mutations,
        yend = median_mutations
      ),
      color = median_color,
      linewidth = 1,
      linetype = "dashed",
      inherit.aes = FALSE
    ) +
    ggplot2::geom_text(
      data = counts_for_plot,
      ggplot2::aes(
        x = x_pos,
        y = y_label_pos,
        label = paste0(
          nonzero_samples,
          "\n",
          total_samples,
          "\n",
          round(nonzero_samples / total_samples, 2)
        )
      ),
      vjust = 2.0,
      size = 3.5,
      lineheight = 0.9,
      inherit.aes = FALSE
    ) +
    ggplot2::scale_x_continuous(
      breaks = seq_along(type_order),
      labels = type_order,
      position = "top"
    ) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::labs(
      title = paste0(signature_name, " mutations"),
      x = NULL,
      y = y_label
    ) +
    ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        angle = 60,
        hjust = 0.0,
        vjust = -0.5
      ),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor.x = ggplot2::element_blank(),
      legend.position = "none",
      plot.margin = ggplot2::margin(t = 5, r = 5, b = 100, l = 5, unit = "pt")
    )

  if (log_scale) {
    p <- p +
      ggplot2::scale_y_log10(
        labels = scales::label_number(drop0trailing = TRUE)
      )
  }

  list(plot = p, sample_counts = sample_counts)
}

plot_signature_hamburger <- function(
  signature_name,
  assignment_matrix,
  sample_info = NULL,
  min_fraction = 0.0,
  ...
) {
  if (!signature_name %in% rownames(assignment_matrix)) {
    stop("Signature '", signature_name, "' not found in assignment matrix")
  }

  sig_values <- as.numeric(assignment_matrix[signature_name, ])
  names(sig_values) <- colnames(assignment_matrix)

  if (min_fraction > 0) {
    sample_totals <- colSums(assignment_matrix)
    fractions <- sig_values / sample_totals
    fractions[is.na(fractions)] <- 0
    sig_values[fractions < min_fraction] <- 0
  }

  plot_signature_by_cancer_type(
    signature_values = sig_values,
    signature_name = signature_name,
    sample_info = sample_info,
    ...
  )
}

message("Loading assignment matrix from: ", assignment_file)
assignment_matrix <- read.delim(
  assignment_file,
  row.names = 1,
  check.names = FALSE
)

# Prefix column names with "CancerType::" using the S17 lookup
ctype <- sample_info$Cancer_Type[
  match(colnames(assignment_matrix), sample_info$Patient)
]
ctype[is.na(ctype)] <- "Unknown"
colnames(assignment_matrix) <- paste0(ctype, "::", colnames(assignment_matrix))

signature_names <- rownames(assignment_matrix)
message(
  "Found ",
  length(signature_names),
  " signatures in assignment matrix"
)

message("Generating PDF: ", output_file)
pdf(output_file, width = params$width, height = params$height)

for (sig_name in signature_names) {
  message("  Plotting: ", sig_name)
  result <- tryCatch(
    {
      plot_signature_hamburger(
        signature_name = sig_name,
        assignment_matrix = assignment_matrix,
        sample_info = sample_info,
        min_fraction = args$min_fraction,
        genome_size_mb = params$genome_size_mb,
        min_samples = params$min_samples,
        log_scale = TRUE,
        exclude_zero = TRUE,
        order_by = "median",
        point_size = params$point_size,
        point_alpha = params$point_alpha,
        base_size = params$base_size
      )
    },
    error = function(e) {
      warning("Failed to plot ", sig_name, ": ", e$message)
      list(plot = NULL)
    }
  )

  if (!is.null(result$plot)) {
    suppressWarnings(print(result$plot))
  }
}

dev.off()
message("Generated: ", output_file)
