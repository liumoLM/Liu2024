#!/bin/bash
# Crop the white page margins from the supplementary figure PDFs into
# sup_figs/cropped/ for sup_figs/sup_figs.qmd. A figure is re-cropped only if
# its source is newer than the cropped copy, or the cropped copy is missing.
# An optional third argument to crop() is a bounding box to cut out.
set -e
cd "$(dirname "$0")/sup_figs"
mkdir -p cropped
crop() {
  if [ ! -e "cropped/$2" ] || [ "$1" -nt "cropped/$2" ]; then
    if [ -n "$3" ]; then
      pdfcrop --margins 4 --bbox "$3" "$1" "cropped/$2" > /dev/null
    else
      pdfcrop --margins 4 "$1" "cropped/$2" > /dev/null
    fi
    echo "cropped $1"
  fi
}
crop Figure_S1a.pdf fig_s1a.pdf
crop Figure_S1b.pdf fig_s1b.pdf
# Figure S2 is shown on two pages, panels A and B on the first and panel C
# on the second. The bounding boxes are in PDF points, origin at the lower
# left.
crop figure_s2.pdf fig_s2ab.pdf "8 369 596 827"
crop figure_s2.pdf fig_s2c.pdf "8 115 596 358"
crop "Figure S3a.pdf" fig_s3a.pdf
crop "Figure S3b.pdf" fig_s3b.pdf
crop "Figure S3c.pdf" fig_s3c.pdf
crop "Figure S4.pdf" fig_s4.pdf
crop figure_s5.pdf fig_s5.pdf
crop "Figure S6.pdf" fig_s6.pdf
