#!/usr/bin/env python3
"""Figures from the ranked neoantigen table -> example_results/figures/."""
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

F = Path("example_results/figures"); F.mkdir(parents=True, exist_ok=True)
d = pd.read_csv(next(Path("example_results/predictions").glob("*.ranked_neoantigen_candidates.tsv")), sep="\t")

fig, axs = plt.subplots(1, 3, figsize=(13, 3.8))
axs[0].hist(np.log10(d.affinity_nm), bins=40, color="#4c72b0")
axs[0].axvline(np.log10(50), c="r", ls="--", lw=1); axs[0].axvline(np.log10(500), c="r", ls=":", lw=1)
axs[0].set_xlabel("log10 predicted affinity (nM)"); axs[0].set_title("Binding affinity (dashed 50 nM, dotted 500 nM)")
c = d.binder_class.value_counts(); axs[1].bar(c.index, c.values, color="#55a868"); axs[1].set_title("Binder class")
a = d.groupby("allele").size().sort_values(); axs[2].barh(a.index, a.values, color="#dd8452"); axs[2].set_title("Candidates per HLA allele")
fig.tight_layout(); fig.savefig(F / "affinity_overview.png", dpi=130); plt.close(fig)

fig, ax = plt.subplots(figsize=(6, 5))
ok = d.wildtype_affinity_nm.notna()
ax.scatter(d.wildtype_affinity_nm[ok], d.affinity_nm[ok], s=5, alpha=.4, c="#4c72b0")
lim = [1, max(d.affinity_nm.max(), d.wildtype_affinity_nm.max())]
ax.plot(lim, lim, c="k", lw=.5); ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xlabel("wild-type peptide affinity (nM)"); ax.set_ylabel("mutant peptide affinity (nM)")
ax.set_title("Mutant vs wild-type binding (below line = gain)"); fig.tight_layout()
fig.savefig(F / "mutant_vs_wildtype.png", dpi=130); plt.close(fig)

top = d.nsmallest(15, "affinity_nm")
fig, ax = plt.subplots(figsize=(7, 4.5))
lab = top.gene_symbol + " " + top.hgvsp.str.split(":").str[-1] + "  " + top.peptide + " (" + top.allele.str.replace("HLA-", "") + ")"
ax.barh(lab[::-1], top.affinity_nm[::-1], color="#c44e52"); ax.set_xlabel("predicted affinity (nM)")
ax.tick_params(axis="y", labelsize=7); ax.set_title("Top 15 candidates"); fig.tight_layout()
fig.savefig(F / "top_candidates.png", dpi=130); plt.close(fig)
print(len(d), "candidate peptide-allele pairs")
