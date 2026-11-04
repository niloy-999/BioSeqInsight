#!/usr/bin/env python3
"""Fail if any \\N{key} used in the manuscript is missing from results/numbers.json.
LaTeX silently expands an undefined key to nothing, which drops the number from the text."""
import json, re, sys, glob
from pathlib import Path
here = Path(__file__).resolve().parent
flat = json.loads((here / "results" / "numbers.json").read_text())["flat"]
used = {}
for f in sorted(glob.glob(str(here / "sections" / "*.tex"))) + [str(here / "main.tex")]:
    for m in re.finditer(r"\\N\{([^}]*)\}", Path(f).read_text()):
        used.setdefault(m.group(1), set()).add(Path(f).name)
missing = {k: sorted(v) for k, v in used.items() if k not in flat}
unused = sorted(set(flat) - set(used))
print(f"{len(used)} distinct \\N keys used, {len(flat)} defined, {len(unused)} defined but unused")
for k, v in sorted(missing.items()): print("MISSING:", k, "in", ", ".join(v))
sys.exit(1 if missing else 0)
