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
## Writes data_flow/data_and_analysis_flow.html. The page draws the graph
## with Cytoscape.js and its dagre layout, loaded from cdnjs and jsdelivr, so
## it needs an internet connection.
## The viewer's controls are described on the page itself.

suppressPackageStartupMessages({
  library(dplyr)
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

## ---- Initial layout ---------------------------------------------------------
## Sugiyama layers (igraph), left to right. Each weakly connected part of the
## graph is laid out on its own and the parts are stacked vertically, largest
## first, so separate pipelines do not interleave. Nodes can be dragged in
## the viewer.
g <- igraph::graph_from_data_frame(edges[, c("from", "to")], vertices = nodes$id)
parts <- igraph::components(g, mode = "weak")
coords <- matrix(NA_real_, nrow(nodes), 2, dimnames = list(nodes$id, NULL))
y_offset <- 0
for (k in order(-parts$csize)) {
  ids <- names(parts$membership)[parts$membership == k]
  sub <- igraph::induced_subgraph(g, ids)
  lay <- if (length(ids) == 1) {
    matrix(c(0, 0), 1)
  } else {
    igraph::layout_with_sugiyama(sub, hgap = 1, vgap = 1, maxiter = 200)$layout
  }
  ## Sugiyama puts layers on y; turn them into columns on x.
  x <- -lay[, 2]
  y <- lay[, 1]
  coords[igraph::V(sub)$name, ] <- cbind(x - min(x), y - min(y) + y_offset)
  y_offset <- y_offset + max(y) - min(y) + 2
}

## ---- Elements as Cytoscape.js JSON -----------------------------------------
cy_nodes <- lapply(seq_len(nrow(nodes)), function(i) {
  n <- nodes[i, ]
  list(
    data = Filter(Negate(is.na), list(
      id = n$id, label = n$label, kind = n$kind, repo = n$repo,
      path = n$path, exists = n$exists, flags = n$flags, note = n$note
    )),
    position = list(x = coords[n$id, 1] * 420, y = coords[n$id, 2] * 40),
    classes = paste(gsub(" ", "-", n$kind), if (!is.na(n$flags)) "flagged")
  )
})
cy_edges <- lapply(seq_len(nrow(edges)), function(i) {
  e <- edges[i, ]
  list(
    data = Filter(Negate(is.na), list(
      id = paste0("e", i), source = e$from, target = e$to,
      relation = e$relation, note = e$note
    )),
    classes = gsub(" ", "-", e$relation)
  )
})
elements_json <- jsonlite::toJSON(c(cy_nodes, cy_edges), auto_unbox = TRUE)

## ---- Page -------------------------------------------------------------------
template <- r"--(<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Data and analysis flow</title>
<script src="https://cdnjs.cloudflare.com/ajax/libs/cytoscape/3.30.2/cytoscape.min.js"></script>
<script src="https://cdnjs.cloudflare.com/ajax/libs/dagre/0.8.5/dagre.min.js"></script>
<script src="https://cdn.jsdelivr.net/npm/cytoscape-dagre@2.5.0/cytoscape-dagre.min.js"></script>
<style>
  :root { --bg: #ffffff; --fg: #1a1a1a; --muted: #5f6b73; --line: #d5dbdf; --panel: #f6f8f9; }
  * { box-sizing: border-box; }
  html, body { margin: 0; height: 100%; background: var(--bg); color: var(--fg);
    font: 14px/1.4 system-ui, -apple-system, "Segoe UI", sans-serif; }
  #app { display: grid; grid-template-columns: 1fr 340px; grid-template-rows: auto 1fr; height: 100vh; }
  header { grid-column: 1 / 3; padding: 8px 12px; border-bottom: 1px solid var(--line);
    display: flex; flex-wrap: wrap; gap: 8px; align-items: center; }
  header h1 { font-size: 16px; margin: 0 12px 0 0; }
  header input[type=search] { width: 340px; padding: 4px 6px; }
  button { padding: 4px 9px; border: 1px solid #aab4ba; background: #fff; border-radius: 4px; cursor: pointer; }
  button:hover { background: #eef2f4; }
  #cy { width: 100%; height: 100%; }
  aside { border-left: 1px solid var(--line); background: var(--panel); padding: 10px 12px; overflow: auto; }
  aside h2 { font-size: 14px; margin: 14px 0 6px; }
  .legend-row { display: flex; align-items: center; gap: 8px; margin: 4px 0; }
  .sw { width: 26px; height: 14px; border: 1.5px solid #333; display: inline-block; }
  .kind-toggle { display: flex; align-items: center; gap: 6px; }
  #info { font-size: 13px; word-break: break-all; }
  #info .flag { color: #c0392b; }
  #help li { margin: 2px 0; }
  .muted { color: var(--muted); }
</style>
</head>
<body>
<div id="app">
<header>
  <h1>Liu et al. 2026: data and analysis flow</h1>
  <input type="search" id="search" list="ids" placeholder="Find a file or script (type part of a path)">
  <datalist id="ids"></datalist>
  <button id="btn-lineage" title="Highlight everything upstream and downstream of the selected nodes">Lineage of selection</button>
  <button id="btn-upstream" title="Only what the selection depends on">Upstream</button>
  <button id="btn-downstream" title="Only what depends on the selection">Downstream</button>
  <button id="btn-isolate" title="Hide everything not highlighted">Hide the rest</button>
  <button id="btn-zoom" title="Zoom to the selected nodes">Zoom to selection</button>
  <button id="btn-relayout" title="Lay out only the nodes now shown, left to right">Re-layout shown</button>
  <button id="btn-fit">Fit all</button>
  <button id="btn-reset">Reset</button>
</header>
<div id="cy"></div>
<aside>
  <h2>Selected</h2>
  <div id="info" class="muted">Click a node or an edge.</div>
  <h2>Show kinds</h2>
  <div id="kinds"></div>
  <h2>Edges</h2>
  <div class="legend-row"><svg width="40" height="10"><line x1="0" y1="5" x2="40" y2="5" stroke="#8c979e" stroke-width="2"/></svg> reads or writes</div>
  <div class="legend-row"><svg width="40" height="10"><line x1="0" y1="5" x2="40" y2="5" stroke="#1f63b5" stroke-width="2" stroke-dasharray="6 4"/></svg> copied</div>
  <div class="legend-row"><svg width="40" height="10"><line x1="0" y1="5" x2="40" y2="5" stroke="#c0392b" stroke-width="2" stroke-dasharray="6 4"/></svg> edited by hand</div>
  <div class="legend-row"><span class="sw" style="border:3px solid #c0392b;background:#fff"></span> flagged: missing, no writer, or never read</div>
  <h2>How to use</h2>
  <ul id="help" class="muted">
    <li>Scroll to zoom, drag the background to pan.</li>
    <li>Drag a node to move it.</li>
    <li>Shift-drag (or Ctrl-drag) on the background to select a box of nodes, then <b>Zoom to selection</b>.</li>
    <li>Click a node to select it; Shift-click or Ctrl-click to add more.</li>
    <li><b>Lineage of selection</b> highlights what the selection depends on and what depends on it; <b>Hide the rest</b> then leaves only that and lays it out compactly.</li>
    <li><b>Re-layout shown</b> lays out the nodes now shown, for example after unchecking some kinds. <b>Reset</b> restores the full graph and its original layout.</li>
  </ul>
</aside>
</div>
<script>
const elements = __ELEMENTS__;
const KINDS = [
  ["data repo file", "data-repo-file", "#8fd3b6", "#1b7f5a", "barrel"],
  ["script", "script", "#e7dcff", "#5b3fa8", "round-rectangle"],
  ["code repo file", "code-repo-file", "#ffd59a", "#b86a00", "ellipse"],
  ["paper file", "paper-file", "#a9d1f5", "#1f63b5", "rectangle"],
  ["external", "external", "#e6e6e6", "#6b6b6b", "diamond"]
];
const style = [
  { selector: "node", style: {
      "label": "data(label)", "font-size": 11, "color": "#111", "text-valign": "center",
      "text-halign": "center", "text-wrap": "ellipsis", "text-max-width": 170,
      "width": "label", "height": 22, "padding": "6px", "border-width": 1.5 } },
  ...KINDS.map(([k, cls, fill, border, shape]) => ({ selector: "node." + cls,
      style: { "background-color": fill, "border-color": border, "shape": shape } })),
  { selector: "node.external", style: { "height": 26, "text-max-width": 150 } },
  { selector: "node.flagged", style: { "border-color": "#c0392b", "border-width": 3 } },
  { selector: "edge", style: { "width": 1.2, "line-color": "#b8c1c7", "target-arrow-color": "#b8c1c7",
      "target-arrow-shape": "triangle", "arrow-scale": 0.8, "curve-style": "bezier" } },
  { selector: "edge.copied", style: { "line-style": "dashed", "line-color": "#1f63b5", "target-arrow-color": "#1f63b5" } },
  { selector: "edge.edited-by-hand", style: { "line-style": "dashed", "line-color": "#c0392b", "target-arrow-color": "#c0392b" } },
  { selector: ".faded", style: { "opacity": 0.12 } },
  { selector: "node.hl", style: { "border-width": 3 } },
  { selector: "edge.hl", style: { "width": 2.2, "line-color": "#47535a", "target-arrow-color": "#47535a" } },
  { selector: "edge.copied.hl", style: { "line-color": "#1f63b5", "target-arrow-color": "#1f63b5" } },
  { selector: "edge.edited-by-hand.hl", style: { "line-color": "#c0392b", "target-arrow-color": "#c0392b" } },
  { selector: ":selected", style: { "overlay-color": "#ffbf00", "overlay-opacity": 0.35, "overlay-padding": 5 } }
];
const cy = cytoscape({
  container: document.getElementById("cy"), elements, style,
  layout: { name: "preset" }, wheelSensitivity: 0.25,
  boxSelectionEnabled: true, selectionType: "additive", minZoom: 0.02, maxZoom: 4
});
cy.fit(undefined, 30);
const original = {};
cy.nodes().forEach(n => { original[n.id()] = Object.assign({}, n.position()); });
cytoscape.use(cytoscapeDagre);
function relayoutShown() {
  const shown = cy.elements(":visible");
  shown.layout({ name: "dagre", rankDir: "LR", nodeSep: 12, rankSep: 90, edgeSep: 6,
    animate: true, animationDuration: 300, fit: true, padding: 30 }).run();
}

// Search box: list every node id.
const dl = document.getElementById("ids");
cy.nodes().map(n => n.id()).sort().forEach(id => {
  const o = document.createElement("option"); o.value = id; dl.appendChild(o); });
document.getElementById("search").addEventListener("change", e => {
  const n = cy.getElementById(e.target.value);
  if (n.nonempty()) { cy.elements().unselect(); n.select(); showInfo(n);
    cy.animate({ center: { eles: n }, zoom: 1.2 }, { duration: 300 }); }
});

// Kind toggles.
const kindsDiv = document.getElementById("kinds");
KINDS.forEach(([k, cls, fill, border]) => {
  const lab = document.createElement("label"); lab.className = "kind-toggle legend-row";
  lab.innerHTML = `<input type="checkbox" checked data-cls="${cls}"> <span class="sw" style="background:${fill};border-color:${border}"></span> ${k} (${cy.nodes("." + cls).length})`;
  kindsDiv.appendChild(lab);
});
kindsDiv.addEventListener("change", e => {
  const cls = e.target.dataset.cls;
  cy.nodes("." + cls).style("display", e.target.checked ? "element" : "none");
});

function esc(s) { return String(s).replace(/[&<>"]/g, c => ({"&":"&amp;","<":"&lt;",">":"&gt;","\"":"&quot;"}[c])); }
function showInfo(el) {
  const d = el.data(), info = document.getElementById("info");
  info.classList.remove("muted");
  if (el.isNode()) {
    const ins = el.incomers("node").length, outs = el.outgoers("node").length;
    info.innerHTML = `<b>${esc(d.id)}</b><br>kind: ${esc(d.kind)}` +
      (d.exists ? `<br>on disk: ${esc(d.exists)}` : "") +
      (d.flags ? `<br><span class="flag">${esc(d.flags)}</span>` : "") +
      (d.note ? `<br>note: ${esc(d.note)}` : "") +
      `<br>${ins} input(s), ${outs} output(s)`;
  } else {
    info.innerHTML = `<b>${esc(d.source)}</b><br>&rarr; ${esc(d.relation)} &rarr;<br><b>${esc(d.target)}</b>` +
      (d.note ? `<br>note: ${esc(d.note)}` : "");
  }
}
cy.on("tap", "node, edge", e => showInfo(e.target));

// Messages go to the side panel; browser alert() pop-ups would block the page.
function say(msg) { const i = document.getElementById("info"); i.classList.add("muted"); i.textContent = msg; }
function selection() {
  const s = cy.nodes(":selected");
  if (s.empty()) say("Select one or more nodes first (click, or Shift-drag a box).");
  return s;
}
function highlight(eles) {
  cy.elements().removeClass("hl").addClass("faded");
  eles.removeClass("faded").addClass("hl");
}
document.getElementById("btn-lineage").onclick = () => { const s = selection(); if (s.nonempty())
  highlight(s.union(s.predecessors()).union(s.successors())); };
document.getElementById("btn-upstream").onclick = () => { const s = selection(); if (s.nonempty())
  highlight(s.union(s.predecessors())); };
document.getElementById("btn-downstream").onclick = () => { const s = selection(); if (s.nonempty())
  highlight(s.union(s.successors())); };
document.getElementById("btn-isolate").onclick = () => {
  const shown = cy.elements(".hl");
  if (shown.empty()) { say("Highlight a lineage first (Lineage, Upstream, or Downstream)."); return; }
  cy.elements().not(shown).style("display", "none");
  relayoutShown();
};
document.getElementById("btn-relayout").onclick = relayoutShown;
document.getElementById("btn-zoom").onclick = () => { const s = selection(); if (s.nonempty())
  cy.animate({ fit: { eles: s, padding: 60 } }, { duration: 300 }); };
document.getElementById("btn-fit").onclick = () =>
  cy.animate({ fit: { eles: cy.elements(":visible"), padding: 30 } }, { duration: 300 });
document.getElementById("btn-reset").onclick = () => {
  cy.elements().removeClass("faded hl").style("display", "element");
  document.querySelectorAll("#kinds input").forEach(i => i.checked = true);
  cy.elements().unselect();
  cy.nodes().forEach(n => n.position(original[n.id()]));
  cy.fit(undefined, 30);
};
</script>
</body>
</html>
)--"

html <- sub("__ELEMENTS__", elements_json, template, fixed = TRUE)
writeLines(html, args$out, useBytes = TRUE)
message("Wrote ", args$out, " (", nrow(nodes), " nodes, ", nrow(edges), " edges)")
