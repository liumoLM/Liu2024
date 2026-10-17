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
# figure_s2_to_s4.pdf holds three panels stacked top to bottom. Each becomes
# its own figure. The bounding boxes (PDF points, origin at the lower left)
# leave out the panel letters A to C at the left edge.
crop figure_s2_to_s4.pdf fig_s2.pdf "22 602 596 804"
crop figure_s2_to_s4.pdf fig_s3.pdf "22 369 596 580"
crop figure_s2_to_s4.pdf fig_s4.pdf "22 115 596 342"
crop "Figure S5.pdf" fig_s5.pdf
crop "Figure S6.pdf" fig_s6.pdf
crop "Figure S7.pdf" fig_s7.pdf
crop "Figure S8.pdf" fig_s8.pdf
crop figure_s9.pdf fig_s9.pdf
crop "Figure S10.pdf" fig_s10.pdf
