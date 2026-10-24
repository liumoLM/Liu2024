# Data and analysis flow

A file-level graph of how the paper is made: which script reads and writes
which file, across the data repo (`unified_indels_data`), the code repo
(`Liu2026_code_and_data`), and this repo. It replaced the PowerPoint deck
`data_and_analysis_flow.pptx` on 2026-10-23.

- `flow_nodes.tsv`: one row per file or script.
- `flow_edges.tsv`: one row per connection.
- `build_flow_viewer.R`: builds `data_and_analysis_flow.html` from the two
  TSVs. Run `Rscript data_flow/build_flow_viewer.R` from the repository root.
- `data_and_analysis_flow.html`: the interactive viewer, built on
  Cytoscape.js (loaded from a CDN, so it needs an internet connection). Drag
  nodes to move them, Shift-drag a box to select nodes and zoom to them,
  find a node by path, highlight the lineage (upstream, downstream, or both)
  of the selection, hide everything else and lay out what is left, and show
  or hide each kind of node. Instructions are on the page.

## Node ids

Each id names the repo and the path in it:

- `data:` a file in the data repo, for example `data:liu_et_al_83_assignment.tsv`
- `code:` a file in the code repo, for example `code:build_fig4/fig4.pdf`
- `paper:` a file in this repo, for example `paper:main_figures/fig4.pdf`
- `ext:` anything outside the three repos, such as files on the Huang-lab
  server (`/public/data/...`) or local VCF folders (`~/MEGA/...`)

`flow_nodes.tsv` columns:

| Column | Meaning |
| --- | --- |
| `id` | as above |
| `kind` | `script`, `data repo file`, `code repo file`, `paper file`, or `external` |
| `repo`, `path`, `label` | the parts of the id, and the file name shown in the viewer |
| `exists` | `yes` or `no` on disk when the TSVs were made; blank for `ext:` |
| `flags` | `missing on disk`, `no script writes it`, `nothing reads it` |
| `note` | for example "made by hand" |

`flow_edges.tsv` columns: `from`, `to`, `relation`, `note`. The relation is
`reads` (file to script), `writes` (script to file), `copied` (a file copied
unchanged, such as a figure copied from the code repo to this repo), or
`edited by hand` (for example `sup_fig_s2_msighdp_signatures_draft.pdf` to
the paper's Figure S2), or `adapted from` (a script whose code was adapted
from another script, such as `read_annotated_vcf()` in
`~/github/zz-liu_2025_draft_release/code_for_internal_exploration/`, which
`build_fig6/make_fig6_logo_counts.R` follows to find the annotated indel
VCFs in `~/MEGA/important_mut_sig_data/`). A script that sources another script `reads` it.

## How the TSVs were made

On 2026-10-23 every script in the code repo (except `not-used/`) and the two
scripts in the data repo were read, and each file a script reads or writes
was listed, including names built with `paste0()` or `file.path()`, files
read by sourced helpers, and argparser defaults. The copies into this repo
come from `main_figures_manifest.tsv`, `sup_figures_manifest.tsv`, and
`sup_tables_manifest.tsv`. Copies made by hand inside the code repo (for
example the topography step 05 and 06 tables copied into `build_fig_s3/`,
which match the pipeline's output to within 4e-13) were added as `copied`
edges.

When a script changes what it reads or writes, edit the two TSVs and rerun
`build_flow_viewer.R`. The viewer leaves out `code:find_data.R`, which most
scripts source only to locate the data repo.
