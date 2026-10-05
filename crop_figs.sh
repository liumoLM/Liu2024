#!/bin/bash
# Crop the white page margins from the main figure PDFs into
# main_figures/cropped/ for ms.qmd. A figure is re-cropped only if its
# source is newer than the cropped copy, or the cropped copy is missing.
set -e
cd "$(dirname "$0")/main_figures"
mkdir -p cropped
crop() {
  if [ ! -e "cropped/$2" ] || [ "$1" -nt "cropped/$2" ]; then
    pdfcrop --margins 4 "$1" "cropped/$2" > /dev/null
    echo "cropped $1"
  fi
}
crop "Graphical abstract.v5.1200px.pdf" graphical_abstract.pdf
crop fig1.pdf fig1.pdf
crop fig2_Cap9_2026.pdf fig2.pdf
crop fig3_Cap9_2026.pdf fig3.pdf
crop fig4_combined.pdf fig4.pdf
crop Figure5.pdf fig5.pdf
crop Figure6-TOP1-TAM.pdf fig6.pdf
