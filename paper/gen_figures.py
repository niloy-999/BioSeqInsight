import json, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

plt.rcParams.update({
    "font.size": 9, "font.family": "sans-serif",
    "axes.spines.top": False, "axes.spines.right": False,
    "figure.dpi": 300,
})
REPO = "/home/claude/v2/BioSeqInsight/benchmarks/results"

# ---------- Figure: fault injection recovery ----------
d = json.load(open(f"{REPO}/fault_injection.json"))
scen = d["scenarios"]
names = list(scen.keys())
labels = [n.replace("_", " ") for n in names]
success = [scen[n]["successful"] for n in names]
failed = [scen[n]["requests"] - scen[n]["successful"] for n in names]

fig, ax = plt.subplots(figsize=(7.0, 3.6))
y = range(len(names))
ax.barh(y, success, color="#2c7a4b", label="Succeeded")
ax.barh(y, failed, left=success, color="#b23b3b", label="Failed")
ax.set_yticks(list(y))
ax.set_yticklabels(labels, fontsize=7.5)
ax.invert_yaxis()
ax.set_xlabel(f"Requests out of {d['repeats_per_scenario']} (max {d['max_attempts']} attempts)")
ax.set_title("Fault-injection outcomes across 13 scripted failure scenarios")
ax.legend(loc="lower right", fontsize=8, frameon=False)
fig.tight_layout()
fig.savefig("figures/fault_injection.pdf")
plt.close(fig)

# ---------- Figure: batch throughput scaling ----------
p = json.load(open(f"{REPO}/performance.json"))
sizes = [row["sequences"] for row in p["batch_scalability"]]
tps = [row["sequences_per_second"] for row in p["batch_scalability"]]
per_seq = [row["per_sequence_ms"] for row in p["batch_scalability"]]

fig, ax1 = plt.subplots(figsize=(4.6, 3.4))
ax1.plot(sizes, tps, marker="o", color="#1f5fa8")
ax1.set_xscale("log")
ax1.set_xlabel("Batch size (sequences)")
ax1.set_ylabel("Throughput (sequences / s)", color="#1f5fa8")
ax1.tick_params(axis="y", labelcolor="#1f5fa8")
ax1.set_xticks(sizes)
ax1.get_xaxis().set_major_formatter(mticker.ScalarFormatter())
ax2 = ax1.twinx()
ax2.plot(sizes, per_seq, marker="s", color="#b23b3b", linestyle="--")
ax2.set_ylabel("Per-sequence cost (ms)", color="#b23b3b")
ax2.tick_params(axis="y", labelcolor="#b23b3b")
ax2.set_ylim(0, max(per_seq) * 1.4)
ax1.set_title("Offline batch throughput vs. batch size")
fig.tight_layout()
fig.savefig("figures/performance.pdf")
plt.close(fig)

# ---------- Figure: live benchmark per-provider outcomes ----------
s = json.load(open(f"{REPO}/structure_benchmark_2026-09-23c_summary.json"))
providers = ["rcsb", "alphafold", "esmatlas"]
disp = {"rcsb": "RCSB PDB", "alphafold": "AlphaFold DB", "esmatlas": "ESM Atlas"}
cats = ["retrieved", "not_found", "service_errors", "invalid_input"]
cat_colors = {"retrieved": "#2c7a4b", "not_found": "#c9a227",
              "service_errors": "#b23b3b", "invalid_input": "#7f7f7f"}
cat_labels = {"retrieved": "Retrieved (usable structure)",
              "not_found": "Coverage gap (not found)",
              "service_errors": "Service error (unrecoverable)",
              "invalid_input": "Refused (over length limit)"}

fig, ax = plt.subplots(figsize=(5.6, 3.6))
bottom = [0]*len(providers)
for cat in cats:
    vals = [s["per_provider"][p][cat] for p in providers]
    ax.bar([disp[p] for p in providers], vals, bottom=bottom,
           color=cat_colors[cat], label=cat_labels[cat])
    bottom = [b+v for b, v in zip(bottom, vals)]
ax.set_ylabel("Candidate proteins (of 69)")
ax.set_title("Live structure-resource benchmark outcomes\n(69 curated proteins, incl. paralogues, orthologues, disordered regions)")
ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.12), ncol=1, fontsize=7.5, frameon=False)
fig.tight_layout()
fig.savefig("figures/live_benchmark.pdf")
plt.close(fig)

# ---------- Figure: mapping-level distribution ----------
dist = s["best_mapping_distribution"]
levels = ["M0", "M1", "M2", "M3", "M4"]
vals = [dist[l] for l in levels]
colors = ["#6e7781", "#cf222e", "#9a6700", "#0969da", "#1a7f37"]

fig, ax = plt.subplots(figsize=(4.2, 3.4))
bars = ax.bar(levels, vals, color=colors)
for b, v in zip(bars, vals):
    if v > 0:
        ax.text(b.get_x()+b.get_width()/2, v+0.6, str(v), ha="center", fontsize=8)
ax.set_ylabel("Proteins (of 69 resolved)")
ax.set_title("Best mapping level achieved per protein\n(across all three structural resources)")
fig.tight_layout()
fig.savefig("figures/mapping_distribution.pdf")
plt.close(fig)

# ---------- Figure: workflow comparison ----------
w = json.load(open(f"{REPO}/workflow_comparison.json"))
tasks = w["tasks"]
task_labels = ["T1\nCharacterise\na sequence", "T2\nRetrieve +\nvalidate structure",
               "T3\nBatch of\n50 sequences", "T4\nReproduce\nan analysis",
               "T5\nShare with a\ncollaborator"]
manual_vals = [t.get("manual_expanded_steps", t["manual"]["steps"]) for t in tasks]
integrated_vals = [t["integrated"]["steps"] for t in tasks]

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.2, 3.2), gridspec_kw={"width_ratios": [1, 1]})
x = list(range(len(tasks)))
ax1.bar([i-0.2 for i in x], manual_vals, width=0.4, color="#b23b3b", label="Manual (multi-tool)")
ax1.bar([i+0.2 for i in x], integrated_vals, width=0.4, color="#2c7a4b", label="BioSeqInsight")
ax1.set_yscale("log")
ax1.set_xticks(x); ax1.set_xticklabels(task_labels, fontsize=6.5)
ax1.set_ylabel("Interaction steps (log scale)")
ax1.legend(fontsize=7, frameon=False)
ax1.set_title("Steps required, log scale\n(T3 needs log scale: 804 vs 3)", fontsize=8.5)

keep = [0, 1, 3, 4]
ax2.bar([i-0.2 for i in range(len(keep))], [manual_vals[k] for k in keep],
        width=0.4, color="#b23b3b")
ax2.bar([i+0.2 for i in range(len(keep))], [integrated_vals[k] for k in keep],
        width=0.4, color="#2c7a4b")
ax2.set_xticks(range(len(keep)))
ax2.set_xticklabels([task_labels[k] for k in keep], fontsize=6.5)
ax2.set_ylabel("Interaction steps (linear)")
ax2.set_title("Same data, linear scale\n(T3 omitted: 804 vs 3, off-scale)", fontsize=8.5)
fig.tight_layout()
fig.savefig("figures/workflow_comparison.pdf")
plt.close(fig)

print("All data-driven figures generated.")
