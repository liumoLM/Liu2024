# Check that the main figures in main_figures/ are not older than their
# sources in the code repository (~/github/Liu2026_code_and_data).
#
# This uses the same comparison as check_sup_figures.R, with the mapping from
# paper figure to source plots in main_figures_manifest.tsv. Subdirectories of
# main_figures/ (cropped/, older/, and so on) are not checked.
#
# Usage:
#   Rscript check_main_figures.R [--code-repo <path>] [--manifest <path>]
#
# Exits with status 1 if any figure is possibly stale, if a file in
# main_figures/ is missing from the manifest, or if a file named in the
# manifest is missing. ms.qmd also calls check_main_figures() at render time
# and warns about any problem.

source("check_sup_figures.R")

#' Compare the main figures with their sources in the code repo.
#'
#' @param code_repo Path to the code repository.
#' @param manifest Path to the manifest that maps paper figures to sources.
#' @param fig_dir Directory holding the paper's main figures.
#' @return A data frame with one row per problem found (zero rows if none).
check_main_figures <- function(
    code_repo = "~/github/Liu2026_code_and_data",
    manifest = "main_figures_manifest.tsv",
    fig_dir = "main_figures") {
  check_sup_figures(code_repo, manifest, sup_dir = fig_dir)
}

if (sys.nframe() == 0) {
  p <- argparser::arg_parser(
    "Check main_figures/ against their sources in the code repo"
  )
  p <- argparser::add_argument(
    p, "--code-repo", help = "Path to the code repository",
    default = "~/github/Liu2026_code_and_data"
  )
  p <- argparser::add_argument(
    p, "--manifest", help = "Manifest mapping paper figures to source plots",
    default = "main_figures_manifest.tsv"
  )
  args <- argparser::parse_args(p)

  problems <- check_main_figures(args$code_repo, args$manifest)
  if (nrow(problems) == 0) {
    message("No main figure is older than its sources.")
  } else {
    message(nrow(problems), " problem(s):")
    message(paste0("  ", problems$file, ": ", problems$problem, collapse = "\n"))
    quit(status = 1)
  }
}
