# Data and analysis flow

A file-level graph of how the paper is made: which script reads and writes
which file, across the data repo (`unified_indels_data`), the code repo
(`Liu2026_code_and_data`), and this repo. It replaced the PowerPoint deck
`data_and_analysis_flow.pptx` on 2026-10-23.

- `flow_nodes.tsv`: one row per file or script.
- `flow_edges.tsv`: one row per connection.
- `cancer_types_nodes.tsv`, `cancer_types_edges.tsv`: how the cancer type
  columns of the data repo's `sample_info.tsv` are made
  (`build_cancer_types/build_cancer_types.R` in the code repo and its three
  sources), kept apart since 2026-11-23. The main graph starts from
  `sample_info.tsv` with its finished cancer types.

The ID89 topography scripts are not in the TSVs. An interactive viewer of
the TSVs (`build_flow_viewer.R` and its HTML pages) was dropped on
2026-11-23, because no layout made the graph readable. It is in the git
history.

## Node ids

Each id names the repo and the path in it:

- `data:` a file in the data repo, for example `data:liu_et_al_83_assignment.tsv`
- `code:` a file in the code repo, for example `code:build_fig4/fig4.pdf`
- `paper:` a file in this repo, for example `paper:main_figures/fig4.pdf`
- `ext:` anything outside the three repos, such as files of the original
  topography run (`/public/data/...`) or local VCF folders (`~/MEGA/...`)

`flow_nodes.tsv` columns:

| Column | Meaning |
| --- | --- |
| `id` | as above |
| `kind` | `script`, `data repo file`, `code repo file`, `paper file`, or `external` |
| `repo`, `path`, `label` | the parts of the id, and the file name |
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
VCFs in `~/MEGA/important_mut_sig_data/`). The annotated indel VCFs are not
raw VCFs. They are made from the raw VCFs (Hartwig `*.purple.somatic.vcf.gz`,
PCAWG `*.consensus.indel.vcf.gz`) by the indel annotation code of the
mSigSpectra package (`read_vcf()` then `annotate_id_vcf()`,
https://github.com/steverozen/mSigSpectra), which adds `seq.context`,
`COSMIC_83`, `Koh_89`, and `Koh_476`. That step is shown in the graph as an
`ext:` node, because no script in the three repos runs it. A script that sources another script `reads` it.

## How the TSVs were made

On 2026-10-23 every script in the code repo (except `not-used/`) and the two
scripts in the data repo were read, and each file a script reads or writes
was listed, including names built with `paste0()` or `file.path()`, files
read by sourced helpers, and argparser defaults. The copies into this repo
come from `main_figures_manifest.tsv`, `sup_figures_manifest.tsv`, and
`sup_tables_manifest.tsv`. Copies made by hand inside the code repo (for
example the topography step 05 and 06 tables once copied into
`build_fig_s3/`) were added as `copied`
edges.

On 2026-11-08 the steps that make the topography inputs were added
(`partial.credit.ID83.R`, `add_strand_bias.R`, SigProfilerSimulator, and
the ICAMS annotation). Their scripts are not yet in the code repo, so they
are `ext:` nodes labeled "(not yet in the code repo)".

The graph is a draft of the end-user documentation. Node and edge notes
describe what a file holds or what a script does. Who sent a file, and
when, belongs in `paper_to_do.md`, not in the TSVs.

On 2026-11-23 the streamed route of the topography analysis was added
(`ID83_topography_analysis_code/streaming/`, see that folder's README,
"How to run"): `rewrite_to_hmf_ids.R` and the copies keyed by HMF IDs,
the local run of `normalize.sh`, `run_streaming_counts.sh` and its awk
scripts and outputs, and `counts_to_matrices.R`, which now writes the r.04
tables in the code repo. The edges from `r.04_create4InfMatrix.R` to those
tables were removed, because they no longer come from it. Also added:
`cancer_subtypes.R` with the scripts that source it. The making of the
cancer types (`build_cancer_types.R`) is in the two cancer type TSVs.

When a script changes what it reads or writes, edit the TSVs.
