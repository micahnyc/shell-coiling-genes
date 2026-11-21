"""
Step 2 (development): are the candidate genes switched on when the shell forms?

Uses NCBI RefSeq gene-level TPM tables (same source as mantle_expression.py) for
bivalves whose annotation includes larval / juvenile RNA-seq:
  Crassostrea virginica  1-day larvae (first shell, prodissoconch I), 22-day juveniles,
                         adult mantle, gonad, gill, hemocytes, whole body
  Ruditapes philippinarum larvae, umbo larvae, mantle, gonad, gill, foot, ...
Sample types are grouped by the "Sample Name" field (tissue/stage) in *_rnaseq_runs.txt.

Genes: candidates from the gene trees (analysis/results/gene_tree_*_assignments.csv):
the mollusc-specific dispatched-like lineage, DISP1/2, DISP3 and diaphanous.
Controls: shell-forming gene names (should be high in larvae and mantle).

Usage:  python analysis/dev_expression.py
Writes analysis/results/dev_expression_long.csv and dev_expression.png
"""

import os
import re
import urllib.request

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..")
RES = os.path.join(HERE, "results")
DATA = os.path.join(ROOT, "data", "dev_expression")

SPECIES = {  # species: (RefSeq assembly dir, annotation release)
    "Crassostrea virginica": ("GCF_053477285.1_ASM5347728v1", "GCF_053477285.1-RS_2025_12"),
    "Ruditapes philippinarum": ("GCF_026571515.1_ASM2657151v2", "GCF_026571515.1-RS_2023_10"),
}
# Mercenaria mercenaria (RS_2023_01) ships raw counts only, no TPM table, so it is left out.
CONTROLS = r"tyrosinase|chitin synthase|carbonic anhydrase|perlucin|nacrein|shell matrix"


def ftp_dir(asm):
    a = asm.split("_")[1].split(".")[0]
    return f"https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/{a[:3]}/{a[3:6]}/{a[6:9]}/{asm}/"


def fetch(url, path):
    if not os.path.exists(path):
        os.makedirs(os.path.dirname(path), exist_ok=True)
        urllib.request.urlretrieve(url, path)
    return path


def tissue_label(name):
    """Sample Name -> short tissue/stage label, e.g. 'mantle', 'larvae (1 day)'."""
    tissue = name.split(" (")[0].split(",")[0].strip().lower()
    m = re.search(r"\(([^)]*)\)\s*$", name)
    extra = m.group(1) if m else ""
    age = re.search(r"(\d+(?:\.\d+)?\s*(?:day|days|year|years))", extra)
    if any(w in tissue for w in ("larva", "juvenile", "embryo")) and age:
        return f"{tissue} ({age.group(1)})"
    return tissue


def main():
    assign = pd.concat([pd.read_csv(os.path.join(RES, f"gene_tree_{f}_assignments.csv"))
                        for f in ("dispatched", "formin")])
    assign = assign[assign.subfamily.isin(["mollusc-specific dispatched-like lineage",
                                           "DISP1/2 (Hedgehog-releasing orthologs)", "DISP3", "diaphanous"])]
    assign["gene_id"] = assign.symbol.str.replace("LOC", "", regex=False)

    rows = []
    for sp, (asm, rel) in SPECIES.items():
        base = ftp_dir(asm)
        d = os.path.join(DATA, rel)
        runs = pd.read_csv(fetch(base + f"{rel}_rnaseq_runs.txt", os.path.join(d, "runs.txt")), sep="\t")
        tpm = pd.read_csv(fetch(base + f"{rel}_normalized_gene_expression_counts.txt.gz",
                                os.path.join(d, "tpm.txt.gz")), sep="\t")
        tpm["gene_id"] = tpm.GeneID.str.replace("GeneID:", "", regex=False)
        runs["label"] = runs["Sample Name"].map(tissue_label)
        cols = [r for r in runs["SRA run accession"] if r in tpm.columns]
        lab = runs.set_index("SRA run accession").label[cols]
        by = tpm[cols].T.groupby(lab.values).mean().T  # genes x sample types
        by.index = tpm.gene_id
        n_runs = lab.value_counts()

        # percentile of each gene within each sample type (expressed genes only)
        pct = by.where(by > 0).rank(pct=True) * 100

        sel = assign[assign.species == sp]
        for g in sel.itertuples():
            if g.gene_id not in by.index:
                continue
            for t in by.columns:
                rows.append(dict(species=sp, group=g.subfamily, gene=g.symbol, name=g.name,
                                 sample_type=t, runs=int(n_runs[t]), TPM=by.at[g.gene_id, t],
                                 pct_within_sample=pct.at[g.gene_id, t]))
        # controls: gene names come from the assembly's feature table (TPM table has only LOC symbols)
        ft = fetch(base + f"{asm}_feature_table.txt.gz", os.path.join(d, "feature_table.txt.gz"))
        f = pd.read_csv(ft, sep="\t", usecols=["# feature", "name", "GeneID"], low_memory=False)
        f = f[(f["# feature"] == "CDS") & f.name.str.contains(CONTROLS, case=False, na=False)]
        ctl_ids = set(f.GeneID.astype(str))
        c = by[by.index.isin(ctl_ids)]
        for t in by.columns:
            rows.append(dict(species=sp, group="shell-gene controls (median)", gene="controls",
                             name=f"{len(c)} genes", sample_type=t, runs=int(n_runs[t]),
                             TPM=c[t].median(), pct_within_sample=pct.loc[c.index, t].median()))
        print(f"{sp}: {len(cols)} runs, {by.shape[1]} sample types, {len(sel)} candidate genes, {len(c)} controls")

    long = pd.DataFrame(rows)
    long.to_csv(os.path.join(RES, "dev_expression_long.csv"), index=False, float_format="%.3f")

    # summary: mean TPM per group x sample type
    summ = long.groupby(["species", "group", "sample_type", "runs"]).agg(
        genes=("gene", "nunique"), mean_TPM=("TPM", "mean"), max_TPM=("TPM", "max"),
        median_pct=("pct_within_sample", "median")).reset_index()
    pd.set_option("display.width", 220)
    for sp in SPECIES:
        s = summ[summ.species == sp].pivot_table(index="group", columns="sample_type", values="max_TPM")
        print(f"\n{sp}: max TPM per group (columns = sample types)\n", s.round(1).to_string())

    # figure: one heatmap per species, rows = candidate genes, columns = sample types
    fig, axes = plt.subplots(len(SPECIES), 1, figsize=(13, 4.2 * len(SPECIES)))
    order = ["shell-gene controls (median)", "mollusc-specific dispatched-like lineage",
             "DISP1/2 (Hedgehog-releasing orthologs)", "DISP3", "diaphanous"]
    for ax, sp in zip(axes, SPECIES):
        s = long[(long.species == sp) & (long.gene != "controls")].copy()
        s["row"] = s.group.str.split(" \\(").str[0] + " | " + s.gene
        s["o"] = s.group.map({g: i for i, g in enumerate(order)})
        m = s.pivot_table(index=["o", "row"], columns="sample_type", values="TPM").sort_index()
        im = ax.imshow(np.log10(m.values + 1), aspect="auto", cmap="YlGnBu", vmin=0, vmax=2.5)
        ax.set_yticks(range(len(m)))
        ax.set_yticklabels([r[1] for r in m.index], fontsize=7)
        ax.set_xticks(range(m.shape[1]))
        ax.set_xticklabels(m.columns, rotation=30, ha="right", fontsize=8)
        for i in range(m.shape[0]):
            for j in range(m.shape[1]):
                ax.text(j, i, f"{m.values[i, j]:.0f}", ha="center", va="center", fontsize=6)
        ax.set_title(f"{sp}: TPM (colour = log10 TPM+1)", fontsize=10, loc="left")
    fig.tight_layout()
    fig.savefig(os.path.join(RES, "dev_expression.png"), dpi=110)
    print("saved dev_expression.png")


if __name__ == "__main__":
    main()
