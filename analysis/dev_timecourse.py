"""
Candidate genes across embryo and larval development, in two species with time courses.

  Mytilus galloprovincialis (mussel), Miglioli et al. 2024, PRJNA996031:
      15 stages, 0-72 h post fertilisation (hpf), 2 replicates. Shell field forms from
      ~16 h, trochophore 20-28 h, D-veliger (first full shell) 44-72 h.
      Counts: NCBI RefSeq GCF_965363235.1-RS_2025_08 gene_expression_counts (raw counts).
  Pomacea canaliculata (apple snail), Accorsi et al. 2025, GEO GSE240084:
      Stages 2-10 and hatchling, 3-4 replicates; STAR counts on RefSeq gene IDs.
      The stage-to-morphology mapping is not given in the deposit (stage N ~ N days).

Expression is counts per million (CPM). It compares one gene across stages, not between genes.
Reference curve: the mean scaled profile (each gene's max = 1) of shell-forming gene
names (tyrosinase, chitin synthase, carbonic anhydrase, ...) that reach >= 5 CPM.
Candidate groups come from the dispatched and formin gene trees.

Usage:  python analysis/dev_timecourse.py
Writes analysis/results/dev_timecourse.csv and dev_timecourse.png
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
CONTROLS = r"tyrosinase|chitin synthase|carbonic anhydrase|perlucin|nacrein|shell matrix"
MYT = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/965/363/235/GCF_965363235.1_xbMytGall1.hap1.1/"
MYT_REL = "GCF_965363235.1-RS_2025_08"
POM = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/073/045/GCF_003073045.1_ASM307304v1/"
POM_COUNTS = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE240nnn/GSE240084/suppl/GSE240084_devel_star_count.csv.gz"
GROUPS = ["mollusc-specific dispatched-like lineage", "DISP1/2 (Hedgehog-releasing orthologs)",
          "DISP3", "diaphanous"]


def fetch(url, name):
    path = os.path.join(DATA, name)
    if not os.path.exists(path):
        os.makedirs(DATA, exist_ok=True)
        urllib.request.urlretrieve(url, path)
    return path


def control_ids(base, asm):
    ft = pd.read_csv(fetch(base + f"{asm}_feature_table.txt.gz", f"{asm}_feature_table.txt.gz"),
                     sep="\t", usecols=["# feature", "name", "GeneID"], low_memory=False)
    ft = ft[(ft["# feature"] == "CDS") & ft.name.str.contains(CONTROLS, case=False, na=False)]
    return set(ft.GeneID.dropna().astype("int64").astype(str))  # older tables read GeneID as float


def cpm(counts):
    return counts / counts.sum() * 1e6


def mytilus():
    runs = pd.read_csv(fetch(MYT + f"{MYT_REL}_rnaseq_runs.txt", "mytilus_runs.txt"), sep="\t")
    runs = runs[runs["Bioproject Accession"] == "PRJNA996031"]
    stage = {r: float(re.search(r"(\d+) hpf", s).group(1))
             for r, s in zip(runs["SRA run accession"], runs["Sample Name"])}
    c = pd.read_csv(fetch(MYT + f"{MYT_REL}_gene_expression_counts.txt.gz", "mytilus_counts.txt.gz"), sep="\t")
    c.index = c.GeneID.str.replace("GeneID:", "", regex=False)
    x = cpm(c[list(stage)])
    by = x.T.groupby(pd.Series(stage)).mean().T  # genes x hpf
    return by, control_ids(MYT, "GCF_965363235.1_xbMytGall1.hap1.1"), "hours post fertilisation"


def pomacea():
    c = pd.read_csv(fetch(POM_COUNTS, "pomacea_GSE240084_counts.csv.gz"), index_col=0)
    c.index = c.index.str.replace("LOC", "", regex=False)
    x = cpm(c)
    st = [s.split("_")[0] for s in x.columns]
    order = {f"Stage-{i}": i for i in range(2, 11)} | {"Hatchling": 11}
    by = x.T.groupby([order[s] for s in st]).mean().T
    return by, control_ids(POM, "GCF_003073045.1_ASM307304v1"), "stage (11 = hatchling)"


def main():
    assign = pd.concat([pd.read_csv(os.path.join(RES, f"gene_tree_{f}_assignments.csv"))
                        for f in ("dispatched", "formin")])
    assign = assign[assign.subfamily.isin(GROUPS)]
    assign["gene_id"] = assign.symbol.str.replace("LOC", "", regex=False)

    rows, fig, axes = [], *plt.subplots(1, 2, figsize=(14, 6.4))
    colours = dict(zip(GROUPS, ["#c0392b", "#3e6f8e", "#8e8e8e", "#b8741a"]))
    for ax, (sp, loader) in zip(axes, [("Mytilus galloprovincialis", mytilus),
                                       ("Pomacea canaliculata", pomacea)]):
        by, ctl, xlabel = loader()
        ctl = by[by.index.isin(ctl) & (by.max(axis=1) >= 5)]
        ref = (ctl.div(ctl.max(axis=1), axis=0)).mean()
        ax2 = ax.twinx()
        ax2.fill_between(ref.index, 0, ref.values, color="#4d7c4f", alpha=0.15, lw=0)
        ax2.set_ylim(0, 1.05)
        ax2.set_ylabel(f"shell-gene reference (mean scaled profile, {len(ctl)} genes)", color="#4d7c4f", fontsize=8)
        sel = assign[(assign.species == sp) & assign.gene_id.isin(by.index)]
        for g in sel.itertuples():
            y = by.loc[g.gene_id]
            ax.plot(y.index, y.values + 0.1, "-o", ms=3, color=colours[g.subfamily], lw=1.4,
                    label=f"{g.subfamily.split(' (')[0]} {g.symbol}")
            for s, v in y.items():
                rows.append(dict(species=sp, group=g.subfamily, gene=g.symbol, stage=s, CPM=v))
        for s, v in ref.items():
            rows.append(dict(species=sp, group="shell-gene reference (scaled)", gene="controls", stage=s, CPM=v))
        ax.set_yscale("log")
        ax.set_xlabel(xlabel)
        ax.set_ylabel("CPM + 0.1 (log scale)")
        ax.set_title(sp, loc="left", fontsize=10)
        ax.set_zorder(ax2.get_zorder() + 1)
        ax.patch.set_visible(False)
        ax.legend(fontsize=7, loc="upper center", bbox_to_anchor=(0.5, -0.14), ncol=2, frameon=False)
        print(f"{sp}: {by.shape[1]} stages, {len(sel)} candidate genes, {len(ctl)} reference genes; "
              f"reference peaks at {ref.idxmax()}")
    axes[0].axvspan(16, 28, color="#999999", alpha=0.12)
    axes[0].text(22, 0.45, "shell field /\ntrochophore", ha="center", fontsize=7)
    fig.tight_layout()
    fig.savefig(os.path.join(RES, "dev_timecourse.png"), dpi=120)
    out = pd.DataFrame(rows)
    out.to_csv(os.path.join(RES, "dev_timecourse.csv"), index=False, float_format="%.3f")
    pd.set_option("display.width", 250)
    print(out[out.gene != "controls"].pivot_table(index=["species", "group", "gene"], columns="stage",
                                                    values="CPM").round(1).to_string())
    print(out[out.gene == "controls"].pivot_table(index="species", columns="stage", values="CPM").round(2).to_string())


if __name__ == "__main__":
    main()
