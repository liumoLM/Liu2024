# Check that the supplementary figures in "Sup Figures/" are not older than
# their sources in the code repository (~/github/Liu2026_code_and_data).
#
# The paper's figures are assembled or edited by hand from the plots the code
# repo produces, so they cannot be compared with them pixel by pixel.
# Instead, a figure is flagged as possibly stale when one of its sources was
# changed after the paper's copy was last changed. A file's change time is its
# last git commit time, or its modification time if it has uncommitted
# changes. The mapping from paper figure to source plots is in
# sup_figures_manifest.tsv, one row per source. Its reviewed_through column
# (YYYY-MM-DD) records that source changes up to and including that date have
# been checked and need no change to the paper's figure.
#
# Usage:
#   Rscript check_sup_figures.R [--code-repo <path>] [--manifest <path>]
#
# Exits with status 1 if any figure is possibly stale, if a file in
# "Sup Figures/" is missing from the manifest, or if a file named in the
# manifest is missing. ms.qmd also calls check_sup_figures() at render time
# and warns about any problem.

#' Time a file's contents last changed, from git or, if uncommitted, from the
#' file system.
#'
#' Commits that only rename or move the file (git status R100) are skipped,
#' because they do not change its contents.
#'
#' @param repo Git repository holding the file.
#' @param path File path relative to `repo`.
#' @return A POSIXct time.
file_change_time <- function(repo, path) {
  dirty <- system2(
    "git", c("-C", shQuote(repo), "status", "--porcelain", "--", shQuote(path)),
    stdout = TRUE
  )
  if (length(dirty) > 0) {
    return(file.mtime(file.path(repo, path)))
  }
  # Output alternates a "@<commit time>" line with name-status lines.
  log <- system2(
    "git", c("-C", shQuote(repo), "log", "--follow", "--format=@%ct",
             "--name-status", "--", shQuote(path)),
    stdout = TRUE
  )
  log <- log[nzchar(log)]
  time <- NA
  for (line in log) {
    if (startsWith(line, "@")) {
      time <- as.numeric(substring(line, 2))
    } else if (!startsWith(line, "R100")) {
      return(as.POSIXct(time, origin = "1970-01-01"))
    }
  }
  file.mtime(file.path(repo, path))
}

#' Compare the supplementary figures with their sources in the code repo.
#'
#' @param code_repo Path to the code repository.
#' @param manifest Path to the manifest that maps paper figures to sources.
#' @param sup_dir Directory holding the paper's supplementary figures.
#' @return A data frame with one row per problem found (zero rows if none).
check_sup_figures <- function(
    code_repo = "~/github/Liu2026_code_and_data",
    manifest = "sup_figures_manifest.tsv",
    sup_dir = "Sup Figures") {
  code_repo <- path.expand(code_repo)
  paper_repo <- system2("git", c("rev-parse", "--show-toplevel"), stdout = TRUE)
  sup_rel <- file.path(
    sub(paste0("^", paper_repo, "/?"), "", normalizePath(sup_dir)), ""
  )
  sup_rel <- sub("^/", "", sup_rel)
  man <- utils::read.delim(manifest, colClasses = "character")
  problems <- data.frame(file = character(), problem = character())
  add <- function(file, problem) {
    problems[nrow(problems) + 1, ] <<- c(file, problem)
  }

  # Files in sup_dir (top level only) that the manifest does not mention,
  # ignoring Office lock and temporary files whose names start with "~".
  present <- list.files(sup_dir)
  present <- present[!dir.exists(file.path(sup_dir, present)) &
                       !startsWith(present, "~")]
  for (f in setdiff(present, man$paper_file)) {
    add(f, "in Sup Figures/ but not in the manifest")
  }

  for (i in seq_len(nrow(man))) {
    paper <- man$paper_file[i]
    src <- man$source_file[i]
    if (!file.exists(file.path(sup_dir, paper))) {
      add(paper, "in the manifest but missing from Sup Figures/")
    } else if (nzchar(src)) {
      if (!file.exists(file.path(code_repo, src))) {
        add(paper, paste("source missing:", src))
      } else {
        paper_time <- file_change_time(paper_repo, paste0(sup_rel, paper))
        src_time <- file_change_time(code_repo, src)
        reviewed <- man$reviewed_through[i]
        if (nzchar(reviewed)) {
          # Changes on or before the reviewed_through date count as reflected.
          paper_time <- max(
            paper_time, as.POSIXct(as.Date(reviewed) + 1) - 1
          )
        }
        if (src_time > paper_time) {
          add(paper, sprintf(
            "source %s changed %s, after the paper's copy (%s)",
            src, format(src_time, "%Y-%m-%d"), format(paper_time, "%Y-%m-%d")
          ))
        }
      }
    }
  }
  problems
}

if (sys.nframe() == 0) {
  p <- argparser::arg_parser(
    "Check Sup Figures/ against their sources in the code repo"
  )
  p <- argparser::add_argument(
    p, "--code-repo", help = "Path to the code repository",
    default = "~/github/Liu2026_code_and_data"
  )
  p <- argparser::add_argument(
    p, "--manifest", help = "Manifest mapping paper figures to source plots",
    default = "sup_figures_manifest.tsv"
  )
  args <- argparser::parse_args(p)

  problems <- check_sup_figures(args$code_repo, args$manifest)
  if (nrow(problems) == 0) {
    message("No supplementary figure is older than its sources.")
  } else {
    message(nrow(problems), " problem(s):")
    message(paste0("  ", problems$file, ": ", problems$problem, collapse = "\n"))
    quit(status = 1)
  }
}
