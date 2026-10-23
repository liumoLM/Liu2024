## Build the interactive data-flow viewer from flow_nodes.tsv and
## flow_edges.tsv.
##
## The two TSVs are the record of how the paper's figures, tables, and
## numbers are made: which script reads and writes which file, across the
## data repo (unified_indels_data), the code repo (Liu2026_code_and_data),
## and this repo. Edit the TSVs, then rerun this script.
##
## Usage, from the repository root:
##   Rscript data_flow/build_flow_viewer.R
##
## Writes data_flow/data_and_analysis_flow.html (self-contained). In the
## viewer, click a node to highlight everything upstream and downstream of
## it, pick a node from the menu at the top left, or filter by kind.

suppressPackageStartupMessages({
  library(dplyr)
  library(visNetwork)
  library(argparser)
})

parser <- arg_parser("Build the data-flow viewer")
parser <- add_argument(
  parser, "--out",
  help = "HTML file to write",
  default = here::here("data_flow", "data_and_analysis_flow.html")
)
args <- parse_args(parser)

nodes <- read.delim(here::here("data_flow", "flow_nodes.tsv"),
                    na.strings = "", colClasses = "character")
edges <- read.delim(here::here("data_flow", "flow_edges.tsv"),
                    na.strings = "", colClasses = "character")
stopifnot(!anyDuplicated(nodes$id),
          all(edges$from %in% nodes$id), all(edges$to %in% nodes$id))

## find_data.R is sourced by most scripts only to locate the data repo. Its
## edges would cross the whole graph, so the viewer leaves it out.
hidden <- "code:find_data.R"
nodes <- nodes %>% dplyr::filter(!id %in% hidden)
edges <- edges %>% dplyr::filter(!from %in% hidden, !to %in% hidden)

## Colors match the earlier PowerPoint deck: green for data, blue for
## outputs, amber for intermediate files.
kind_style <- tibble::tribble(
  ~kind,            ~shape,     ~color,
  "data repo file", "database", "#2E7D4F",
  "script",         "box",      "#3B4A54",
  "code repo file", "ellipse",  "#C98A1B",
  "paper file",     "ellipse",  "#1F5F8B",
  "external",       "ellipse",  "#7A8A94"
)

esc <- function(x) htmltools::htmlEscape(ifelse(is.na(x), "", x))

vis_nodes <- nodes %>%
  left_join(kind_style, by = "kind") %>%
  mutate(
    group = kind,
    font.color = ifelse(kind == "script", "white", "black"),
    ## Files flagged as missing, without a writer, or unread get a red
    ## border.
    borderWidth = ifelse(is.na(flags), 1, 3),
    color.border = ifelse(is.na(flags), color, "#C0392B"),
    color.background = ifelse(kind == "script", color,
                              paste0(color, "33")),
    title = paste0(
      "<b>", esc(path), "</b><br>",
      "repo: ", esc(repo), "<br>kind: ", esc(kind),
      ifelse(is.na(exists), "", paste0("<br>on disk: ", esc(exists))),
      ifelse(is.na(flags), "", paste0("<br><span style='color:#C0392B'>",
                                      esc(flags), "</span>")),
      ifelse(is.na(note), "", paste0("<br>note: ", esc(note)))
    )
  ) %>%
  ## Font color of non-script nodes is set per row above; boxes use white.
  mutate(font.color = ifelse(kind == "script", "white", "#1F2A30")) %>%
  select(id, label, group, shape, title, borderWidth, color.border,
         color.background, font.color)

vis_edges <- edges %>%
  mutate(
    arrows = "to",
    dashes = relation %in% c("copied", "edited by hand"),
    color = case_when(
      relation == "edited by hand" ~ "#C0392B",
      relation == "copied" ~ "#1F5F8B",
      TRUE ~ "#9AA5AD"
    ),
    title = paste0(esc(relation), ifelse(is.na(note), "", paste0(": ", esc(note))))
  ) %>%
  select(from, to, arrows, dashes, color, title)

## Layout: Sugiyama layers (igraph), left to right. Each weakly connected
## part of the graph is laid out on its own and the parts are stacked
## vertically, largest first, so separate pipelines do not interleave.
g <- igraph::graph_from_data_frame(vis_edges[, c("from", "to")],
                                   vertices = vis_nodes$id)
parts <- igraph::components(g, mode = "weak")
part_order <- order(-parts$csize)
coords <- matrix(NA_real_, nrow(vis_nodes), 2,
                 dimnames = list(vis_nodes$id, NULL))
y_offset <- 0
for (k in part_order) {
  ids <- names(parts$membership)[parts$membership == k]
  sub <- igraph::induced_subgraph(g, ids)
  lay <- if (length(ids) == 1) {
    matrix(c(0, 0), 1)
  } else {
    igraph::layout_with_sugiyama(sub, hgap = 1, vgap = 1, maxiter = 200)$layout
  }
  ## Sugiyama puts layers on y; turn them into columns on x.
  x <- lay[, 2] * -1
  y <- lay[, 1]
  x <- x - min(x)
  y <- y - min(y)
  coords[igraph::V(sub)$name, ] <- cbind(x, y + y_offset)
  y_offset <- y_offset + max(y) + 2
}
vis_nodes$x <- coords[vis_nodes$id, 1] * 800
vis_nodes$y <- coords[vis_nodes$id, 2] * 45

net <- visNetwork(vis_nodes, vis_edges, width = "100%", height = "92vh",
                  main = "Liu et al. 2026: data and analysis flow",
                  submain = paste(
                    "Click a node to highlight what it depends on and what",
                    "depends on it. Dashed edges are copies; red dashed",
                    "edges are hand edits. A red border marks a file that is",
                    "missing, has no writer, or is never read."
                  )) %>%
  visNodes(physics = FALSE) %>%
  visPhysics(enabled = FALSE) %>%
  visEdges(smooth = list(type = "cubicBezier", forceDirection = "horizontal")) %>%
  visOptions(
    highlightNearest = list(enabled = TRUE, algorithm = "hierarchical",
                            degree = list(from = 100, to = 100),
                            hover = FALSE),
    ## Full ids (code:..., data:..., paper:...) in the menu, because
    ## short labels such as fig4.pdf repeat across repos.
    nodesIdSelection = list(enabled = TRUE, useLabels = FALSE,
                            main = "Select a file or script"),
    selectedBy = list(variable = "group", main = "Select by kind")
  ) %>%
  visInteraction(navigationButtons = TRUE, tooltipDelay = 100) %>%
  visLegend(
    addNodes = kind_style %>%
      transmute(label = kind, shape, color = ifelse(kind == "script", color,
                                                    paste0(color, "33")),
                font.color = "black"),
    useGroups = FALSE, position = "right", width = 0.12
  )

visSave(net, file = args$out, selfcontained = TRUE)
## The HTML embeds its libraries; this folder is left over from saving.
unlink(sub("\\.html$", "_files", args$out), recursive = TRUE)
message("Wrote ", args$out, " (", nrow(vis_nodes), " nodes, ",
        nrow(vis_edges), " edges)")
