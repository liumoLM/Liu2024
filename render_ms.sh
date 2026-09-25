#!/bin/bash
# Re-crop the main figures, then render ms.qmd to ms.pdf. Pass --open to open the result in the default viewer.
set -e
cd "$(dirname "$0")"
./crop_figs.sh
quarto render ms.qmd --to pdf
if [ "$1" = "--open" ]; then
  xdg-open ms.pdf > /dev/null 2>&1 &
fi
