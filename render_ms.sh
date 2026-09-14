#!/bin/bash
# Render ms.qmd to ms.pdf. Pass --open to open the result in the default viewer.
set -e
cd "$(dirname "$0")"
quarto render ms.qmd --to pdf
if [ "$1" = "--open" ]; then
  xdg-open ms.pdf > /dev/null 2>&1 &
fi
