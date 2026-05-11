# To do list for shiny app

## Critical

- On two test files, did not complete building 476-type matrix before server disconnected; not sure if this because there was a problem processing the 
  the test files or if thing just timed out.
  
  - test file 1 (hg38): https://raw.githubusercontent.com/steverozen/mSigSpectra/refs/heads/main/inst/reijns_cell_data/reijns_cell_AKO1_4.vcf
  - test file 2 (hg37): https://raw.githubusercontent.com/steverozen/mSigSpectra/refs/heads/main/inst/extdata/Strelka-ID-GRCh37/Strelka.ID.GRCh37.s1.vcf

Make copies of these files in your repo and use them for testing. They are cell line data so no privacy concerns.

Please use mSigSpectra for this -- it should be faster, and we need to publish mSigSpectra as part of this paper.

## Presentation etc.

- On the front page, put order 476, 89, 83.

- On the front page, for 476, change the text to:
 "This new indel classification (Koh et all, 2025) offers a highly granular classification of 476 types
  of indel. It provides the highest resolution for analyzing similarities and differences
  between signatures.

- Please regenerate the 89-type and 476-type plots using the latest mSigPlot from Github; after thinking
  about it a lot I decide we should open "righ-open" intervals e.g. Ins2:U2:R(5,), not Ins2:U2:R(5,). 
  The reason is that if use  Ins2:U2:R(5,). and then start to analyze e.g. the Hartwig data
  in isolation, we would need to add a new category, Ins2:U2:R(10,). I don't think we should go into
  the area of new classificaiton systems in this paper.

- Put a link to the shiny app on the zenodo front page

- On the shiny app front page, remove the "Algorithmic tranlsation ... " block and change 
"This web site presents signatures extracted..." to "This web site presents corrsponding signatures extracted..."



