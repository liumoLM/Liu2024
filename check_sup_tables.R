# Check that the supplementary tables in "Sup Tables/" match their sources in
# the code repository (~/github/Liu2026_code_and_data).
#
# The mapping from paper file name to source file is in
# sup_tables_manifest.tsv, because the two repos name the files differently.
# Tables are compared by cell contents, sheet by sheet, not byte by byte,
# because two .xlsx files with identical cells can differ in zip metadata.
#
# Usage:
#   Rscript check_sup_tables.R [--code-repo <path>] [--manifest <path>]
#
# Exits with status 1 if any table with a source is out of date, if a file
# in "Sup Tables/" is missing from the manifest, or if a file named in the
# manifest is missing. ms.qmd also calls check_sup_tables() at render time
# and warns about any problem.

#' Read every sheet of an .xlsx file, or a delimited text file, as text.
#'
#' @param path File to read.
#' @return A named list of data frames, one per sheet.
read_table_cells <- function(path) {
  if (grepl("\\.xlsx$", path, ignore.case = TRUE)) {
    sheets <- readxl::excel_sheets(path)
    out <- lapply(sheets, function(s) {
      as.data.frame(suppressMessages(readxl::read_excel(
        path, sheet = s, col_names = FALSE, col_types = "text"
      )))
    })
    names(out) <- sheets
    out
  } else {
    list(data = utils::read.delim(
      path, header = FALSE, colClasses = "character", check.names = FALSE
    ))
  }
}

#' Compare the supplementary tables with their sources in the code repo.
#'
#' @param code_repo Path to the code repository.
#' @param manifest Path to the manifest that maps paper files to source files.
#' @param sup_dir Directory holding the paper's supplementary tables.
#' @return A data frame with one row per problem found (zero rows if none).
check_sup_tables <- function(
    code_repo = "~/github/Liu2026_code_and_data",
    manifest = "sup_tables_manifest.tsv",
    sup_dir = "Sup Tables") {
  man <- utils::read.delim(manifest, colClasses = "character")
  problems <- data.frame(file = character(), problem = character())
  add <- function(file, problem) {
    problems[nrow(problems) + 1, ] <<- c(file, problem)
  }

  # Files in sup_dir (top level only) that the manifest does not mention.
  present <- list.files(sup_dir)
  present <- present[!dir.exists(file.path(sup_dir, present))]
  for (f in setdiff(present, man$paper_file)) {
    add(f, "in Sup Tables/ but not in the manifest")
  }

  for (i in seq_len(nrow(man))) {
    paper <- file.path(sup_dir, man$paper_file[i])
    if (!file.exists(paper)) {
      add(man$paper_file[i], "in the manifest but missing from Sup Tables/")
    } else if (nzchar(man$source_file[i])) {
      src <- file.path(path.expand(code_repo), man$source_file[i])
      if (!file.exists(src)) {
        add(man$paper_file[i], paste("source missing:", man$source_file[i]))
      } else if (!identical(read_table_cells(paper), read_table_cells(src))) {
        add(man$paper_file[i], paste("differs from", man$source_file[i]))
      }
    }
  }
  problems
}

if (sys.nframe() == 0) {
  p <- argparser::arg_parser(
    "Check Sup Tables/ against their sources in the code repo"
  )
  p <- argparser::add_argument(
    p, "--code-repo", help = "Path to the code repository",
    default = "~/github/Liu2026_code_and_data"
  )
  p <- argparser::add_argument(
    p, "--manifest", help = "Manifest mapping paper files to source files",
    default = "sup_tables_manifest.tsv"
  )
  args <- argparser::parse_args(p)

  problems <- check_sup_tables(args$code_repo, args$manifest)
  man <- utils::read.delim(args$manifest, colClasses = "character")
  n_checked <- sum(nzchar(man$source_file))
  if (nrow(problems) == 0) {
    message(
      "All ", n_checked, " supplementary tables with a source match it. ",
      sum(!nzchar(man$source_file)), " files have no source in the code repo."
    )
  } else {
    message(nrow(problems), " problem(s):")
    message(paste0("  ", problems$file, ": ", problems$problem, collapse = "\n"))
    quit(status = 1)
  }
}
