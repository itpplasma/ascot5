"""Plot synthetic orbit coverage: python plot.py (matplotlib required)."""
from pathlib import Path
import json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

root = Path(__file__).resolve().parent
data = json.loads((root / "coverage.json").read_text())
ids = data["marker_id"]
fig, axes = plt.subplots(1, 3, figsize=(10.5, 3.6), sharey=True,
                         layout="constrained")
for ax, ranks in zip(axes, (1, 2, 4)):
    counts = data["mpi_ranks"][str(ranks)]
    ax.axhline(22, color="0.2", linestyle=":", linewidth=2,
               label="Expected (fixed sampling)")
    ax.plot(ids, counts["before"], "s--", color="#D55E00", markersize=7,
            markerfacecolor="none", label="Before")
    ax.plot(ids, counts["after"], "o-", color="#0072B2", markersize=4,
            label="After")
    ax.set(title="%d MPI rank%s" % (ranks, "" if ranks == 1 else "s"),
           xlabel="Marker ID [1]", xticks=ids, ylim=(-1.5, 29),
           yticks=[0, 11, 22])
    ax.grid(axis="y", alpha=.2)
axes[0].set_ylabel("Orbit samples per marker [1]")
handles, labels = axes[0].get_legend_handles_labels()
fig.legend(handles, labels, loc="outside upper center", ncol=3,
           fontsize=9, frameon=False)
fig.savefig(root / "orbit-coverage.png", dpi=180)
