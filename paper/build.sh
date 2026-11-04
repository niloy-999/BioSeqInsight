#!/usr/bin/env bash
# Rebuild the manuscript from the validation study's generated numbers.
#
#   ./build.sh [path/to/BioSeqInsight/validation/results/manuscript]
#
# With an argument, copies numbers.json and numbers.tex from there first. Then:
#   1. check_numbers.py   fails if the text uses a number key that does not exist
#                         (LaTeX would silently drop it)
#   2. gen_validation_figures.py   draws the validation figures from numbers.json
#   3. pdflatex + bibtex  (needs a TeX installation and matplotlib/numpy)
set -euo pipefail
cd "$(dirname "$0")"
if [ "${1:-}" != "" ]; then
    mkdir -p results
    cp "$1/numbers.json" "$1/numbers.tex" results/
    echo "copied numbers from $1"
fi
python3 check_numbers.py
python3 gen_validation_figures.py
pdflatex -interaction=nonstopmode main.tex >/dev/null
bibtex main >/dev/null || true
pdflatex -interaction=nonstopmode main.tex >/dev/null
pdflatex -interaction=nonstopmode main.tex >/dev/null
echo "built main.pdf ($(pdfinfo main.pdf | awk '/Pages/{print $2}') pages)"
grep -q "PROVISIONAL" <(pdftotext -f 1 -l 1 main.pdf -) && echo "NOTE: this build still carries the PROVISIONAL banner (numbers.tex is provisional)" || echo "final numbers in use"
