# Indel Paper

**Status:** In progress — text work this morning (2026-04-22)
**Relative priority:** Secondary to the mutational signatures grant; can continue into late May or later per `contexts/current-priorities.md`.

## Overview

Paper on indel mutational signatures / new indel classification work.

## Open tasks (from current-priorities.md)

- Main Idea: run through paper in order and carry out 
- To do list for Mo
- 
- Clean up text in shiny app

- Move Signature F to the vignette; release / ask for citation (include lovelace_luquette@hms.harvard.edu)
- Look at code base and think more about what goes in it
- Add mSigPlot to CRAN
- Finish mSigSpectra and submit to CRAN
    - Added code to translate?
- Update MTIS slides and publish; include "takeaways"

## Waiting
Updates for shiny app from xueming
mSigPlot CI results, then CRAN


## Paper deliverables

- Vignette

## Related / potential grant deliverables

(Seed list from current-priorities.md — for the mutational-signatures grant, not necessarily the paper)
- Generic mutation simulator in C++ (any classification; opportunities are classification-dependent)
- Plotting for new classifications — how to make these generic?
- Cluster samples by activity?
- SVD or other dimension reduction, then look for signatures in reduced dimension
- Synthetic data (VAE / GAN / latent Dirichlet allocation)
- Lack of benchmarking

## Repo

Leave the to do list here for now

## Memory / Note

Paper shiny app at: 
Review shiny app: https://indelsignaturebrowser.shinyapps.io/indel-signature-browser/?nav=Home



## Literature notes

- [`2026-08-11-replication-slippage-poly-T-tracts.md`](2026-08-11-replication-slippage-poly-T-tracts.md)
  — replication slippage at homopolymers. Two results bear directly on this paper:
  (1) ID1 (+T at poly-T) vs ID2 (-T at poly-T) are ascribed to nascent- vs
  template-strand slippage, and Koh et al. 2025 separates them experimentally by
  tract length (Pol-dys insertions peak at poly-T5-7, MMRd deletions grow at poly-T6+);
  (2) there is no evidence that A:T slips more than G:C at matched tract length. The
  poly-A/T dominance of human indel spectra is a target-abundance and tract-length
  effect, which matters for how opportunity normalization is done in any new
  classification.
- [`comparison-with-koh-2025.md`](comparison-with-koh-2025.md) — how our extraction
  differs from Koh et al. 2025.
