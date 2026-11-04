#!/usr/bin/env python3
"""Draw the validation-study figures from results/numbers.json (written by
validation/scripts/07_manuscript_numbers.py). No number is typed here."""
import json, sys
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np

HERE = Path(__file__).resolve().parent
data = json.loads((HERE / "results" / "numbers.json").read_text())
raw = data["raw"]
plt.rcParams.update({"font.size": 8, "axes.spines.top": False, "axes.spines.right": False, "font.family": "DejaVu Sans"})
C = {"true": "#2a7f62", "same": "#e0b252", "near": "#e08a3c", "diff": "#c0392b", "unm": "#8e8e8e", "none": "#e9e9e9"}

# ---- Figure A: outcomes of the four acceptance policies over all 1,000 queries
policies = [("method1_naive", "Naive: accept top RCSB hit"), ("method2_naive", "Naive: top hit with RCSB 0.9 identity cutoff"),
            ("method1_validated_permissive", "Validated, permissive (M3 or better)"), ("method1_validated_strict", "Validated, strict (M4)")]
fig, ax = plt.subplots(figsize=(6.8, 2.9))
for i, (m, label) in enumerate(policies):
    b = raw["report_rcsb"][m]["all"]; t = b["false_association_tiers"]; n = 1000
    parts = [("true", b["true_associations"], "Correct structure"), ("same", t["same_sequence"], "Wrong accession, same sequence"),
             ("near", t["near_identical"], "Wrong accession, 90-99% identical"), ("diff", t["different"], "Different protein (<90%)"),
             ("unm", t["unmeasured"], "Wrong accession, not measurable")]
    used = sum(v for _, v, _ in parts); parts.append(("none", n - used, "No structure accepted"))
    left = 0
    for key, v, lab in parts:
        ax.barh(i, v, left=left, color=C[key], edgecolor="white", linewidth=0.4, label=lab if i == 0 else None)
        if v >= 45: ax.text(left + v / 2, i, str(v), ha="center", va="center", fontsize=7, color="white" if key in ("true", "diff") else "black")
        left += v
ax.set_yticks(range(len(policies))); ax.set_yticklabels([l for _, l in policies]); ax.invert_yaxis()
ax.set_xlim(0, 1000); ax.set_xlabel("Queries (of 1,000 random Swiss-Prot proteins)")
ax.legend(ncol=3, fontsize=6.3, loc="upper center", bbox_to_anchor=(0.32, -0.30), frameon=False, columnspacing=1.0)
fig.tight_layout(); fig.savefig(HERE / "figures" / "val_outcomes.pdf", bbox_inches="tight"); plt.close(fig)

# ---- Figure B: what the strict gate did with the correct naive hits
g = raw["gate_outcomes"]
labels = {"accepted": "Accepted (M4)", "partial": "Rejected: partial structure of the protein", "unmeasured": "Rejected: no PDB-format file to validate",
          "M3": "Rejected: M3 (accession not confirmed)", "other": "Rejected: other (M1/M2)"}
vals = {"accepted": 0, "partial": 0, "unmeasured": 0, "M3": 0, "other": 0}
for k, v in g.items():
    key = "accepted" if k.startswith("accepted") else "partial" if k.startswith("partial") else "unmeasured" if k.startswith("unmeasured") else "M3" if k.startswith("M3") else "other"
    vals[key] += v
fig, ax = plt.subplots(figsize=(5.4, 2.0))
cols = [C["true"], "#4c78a8", C["unm"], C["same"], "#b9a7c9"]
for i, (k, lab) in enumerate(labels.items()):
    ax.barh(i, vals[k], color=cols[i]); ax.text(vals[k] + 1.5, i, str(vals[k]), va="center", fontsize=7)
ax.set_yticks(range(len(labels))); ax.set_yticklabels(list(labels.values())); ax.invert_yaxis()
ax.set_xlabel(f"Correct naive top hits (n = {sum(vals.values())})"); ax.set_xlim(0, max(vals.values()) * 1.15)
fig.tight_layout(); fig.savefig(HERE / "figures" / "val_gate_cost.pdf"); plt.close(fig)

# ---- Figure C: adversarial classes by mapping level
adv = raw["adversarial"]["classes"]
order = [("B1", "B1\nwrong paralogue"), ("B2", "B2\nfragment of\nthe protein"), ("B4", "B4\nunrelated"), ("B5", "B5\nnear-identical,\nother entry"), ("B3", "B3\nisoform\n(labels\nunverified)")]
lev_col = {"M0": "#dddddd", "M1": "#9ecae1", "M2": "#6baed6", "M3": "#e0b252", "M4": "#c0392b"}
fig, ax = plt.subplots(figsize=(6.0, 2.8))
for i, (cls, lab) in enumerate(order):
    blk = adv[cls]; n = max(blk["n"], 1); bottom = 0
    for lv in ("M1", "M2", "M3", "M4"):
        h = 100 * blk["levels"][lv] / n
        ax.bar(i, h, bottom=bottom, color=lev_col[lv], edgecolor="white", linewidth=0.4, hatch="//" if cls == "B3" else None, label=lv if i == 0 else None)
        if h >= 7: ax.text(i, bottom + h / 2, str(blk["levels"][lv]), ha="center", va="center", fontsize=6.5)
        bottom += h
    ax.text(i, 103, f"n = {blk['n']}", ha="center", fontsize=6.5)
ax.set_xticks(range(len(order))); ax.set_xticklabels([l for _, l in order], fontsize=6.5); ax.set_ylim(0, 112)
ax.set_ylabel("Cases (%)"); ax.legend(title="Mapping level", fontsize=6.5, title_fontsize=7, ncol=4, loc="lower center", bbox_to_anchor=(0.5, 1.0), frameon=False)
fig.tight_layout(); fig.savefig(HERE / "figures" / "val_adversarial.pdf"); plt.close(fig)

# ---- Figure D: identity x coverage grids, ground-truth positives vs negatives
gp, gn = np.array(raw["grid_positive"]), np.array(raw["grid_negative"])
fig, axes = plt.subplots(1, 2, figsize=(6.4, 2.9), sharey=True)
vmax = max(gp.max(), gn.max())
for ax, grid, title in ((axes[0], gp, "Correct hits (ground truth)"), (axes[1], gn, "Wrong hits and constructed negatives")):
    im = ax.imshow(np.where(grid == 0, np.nan, grid), origin="lower", aspect="auto", cmap="viridis", norm=LogNorm(1, vmax), extent=[0, 100, 0, 100])
    ax.set_title(title, fontsize=7.5); ax.set_xlabel("Coverage of query (%)")
    th = raw["thresholds"]
    ax.axvline(th["exact_coverage"], color="white", ls=":", lw=0.8); ax.axhline(th["exact_identity"], color="white", ls=":", lw=0.8)
axes[0].set_ylabel("Identity over aligned region (%)")
cb = fig.colorbar(im, ax=axes, fraction=0.04, pad=0.02); cb.set_label("Hits (log scale)", fontsize=7)
fig.savefig(HERE / "figures" / "val_grid.pdf", bbox_inches="tight"); plt.close(fig)
print("figures written")
