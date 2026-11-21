"""
Are the QTL-region genes switched on in the mantle, the tissue that builds the shell?

Data: NCBI's RefSeq annotation of the 2024 Littorina saxatilis genome
(GCF_037325665.1-RS_2024_12) ships gene-level TPMs for every RNA-seq run it used.
Ten of those runs are a tissue panel (BioProject PRJNA550990): mantle x2, foot x2,
hepatopancreas x2, reproductive organs x2, penis, head.

For each gene:
  mantle_TPM     mean TPM of the two mantle samples
  other_TPM      mean TPM of the tissue means for the five other tissues
  mantle_log2FC  log2((mantle_TPM + 1) / (other_TPM + 1))
  mantle_pct     percentile of mantle_log2FC among expressed genes genome-wide
                 (100 = most mantle-biased gene in the genome)
With one or two samples per tissue there is no replication worth a statistical
test; this is a ranking, not a differential-expression analysis.

Usage:  python analysis/mantle_expression.py
Writes analysis/results/mantle_expression_genome.csv.gz,
       mantle_expression_qtl_genes.csv, mantle_expression.png
"""

import os
import re

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..")
RES = os.path.join(HERE, "results")
EXP = os.path.join(ROOT, "data", "littorina_v2", "expression")
TPM = os.path.join(EXP, "GCF_037325665.1-RS_2024_12_normalized_gene_expression_counts.txt.gz")
RUNS = os.path.join(EXP, "GCF_037325665.1-RS_2024_12_rnaseq_runs.txt")
PANEL = "PRJNA550990"
EXPRESSED = 1.0  # TPM in at least one panel tissue

# Known shell-forming gene names: if the mantle data behave, these should rank high.
CONTROLS = r"tyrosinase|chitin synthase|carbonic anhydrase|perlucin|shell matrix|nacrein|pif|collagen alpha-1\(XII\)"


def tissue_of(sample_name):
    return sample_name.split(" (")[0]


def main():
    runs = pd.read_csv(RUNS, sep="\t")
    runs = runs[runs["Bioproject Accession"] == PANEL]
    tissue = dict(zip(runs["SRA run accession"], runs["Sample Name"].map(tissue_of)))

    t = pd.read_csv(TPM, sep="\t")
    t["gene_id"] = t["GeneID"].str.replace("GeneID:", "", regex=False)
    by_tissue = t[list(tissue)].T.groupby(pd.Series(tissue)).mean().T  # genes x tissues
    by_tissue.index = t.gene_id

    others = [c for c in by_tissue.columns if c != "mantle"]
    g = pd.DataFrame(index=by_tissue.index)
    g["gene"] = t.GeneSym.values
    g["chrom"], g["start"] = t.Chr.values, t.GeneStart.values
    g["mantle_TPM"] = by_tissue["mantle"]
    g["other_TPM"] = by_tissue[others].mean(axis=1)
    g["top_tissue"] = by_tissue.idxmax(axis=1)
    g["mantle_log2FC"] = np.log2((g.mantle_TPM + 1) / (g.other_TPM + 1))
    expressed = by_tissue.max(axis=1) >= EXPRESSED
    g["expressed"] = expressed
    g.loc[expressed, "mantle_pct"] = 100 * g.loc[expressed, "mantle_log2FC"].rank(pct=True)
    for c in by_tissue.columns:
        g[f"TPM_{c}"] = by_tissue[c]

    # gene descriptions from the GFF already used in map_qtl_genes.py
    desc = {}
    with open(os.path.join(ROOT, "data", "littorina_v2", "genomic.gff"), encoding="utf-8") as f:
        for line in f:
            if "\tgene\t" not in line:
                continue
            a = dict(kv.split("=", 1) for kv in line.rstrip("\n").split("\t")[8].split(";") if "=" in kv)
            m = re.search(r"GeneID:(\d+)", a.get("Dbxref", ""))
            if m:
                desc[m.group(1)] = a.get("description", "").replace("%2C", ",")
    g["description"] = g.index.map(desc)
    g.to_csv(os.path.join(RES, "mantle_expression_genome.csv.gz"), float_format="%.3f")

    ex = g[g.expressed]
    print(f"{len(ex)} of {len(g)} genes expressed (>= {EXPRESSED} TPM in some panel tissue)")
    print("Tissues:", list(by_tissue.columns))

    ctl = ex[ex.description.fillna("").str.contains(CONTROLS, case=False, regex=True)]
    print(f"\nPositive controls (shell-forming gene names, n={len(ctl)}): median mantle percentile "
          f"{ctl.mantle_pct.median():.0f} (genome median 50)")
    print(ctl.sort_values("mantle_pct", ascending=False)
          [["gene", "description", "mantle_TPM", "other_TPM", "mantle_pct"]].head(12).to_string())

    # QTL regions
    q = pd.read_csv(os.path.join(RES, "qtl_region_genes.csv"), dtype={"gene_id": str})
    q = q[q.biotype == "protein_coding"]
    m = q.merge(g.drop(columns=["gene", "description"]), left_on="gene_id", right_index=True, how="left")
    m.to_csv(os.path.join(RES, "mantle_expression_qtl_genes.csv"), index=False, float_format="%.3f")

    print("\nQTL regions: share of expressed genes in the genome's top 5% / 10% most mantle-biased")
    rows = []
    for reg, r in m.groupby("region"):
        e = r[r.expressed == True]
        rows.append(dict(region=reg, genes=len(r), expressed=len(e),
                         top5=int((e.mantle_pct >= 95).sum()), top10=int((e.mantle_pct >= 90).sum()),
                         median_pct=round(e.mantle_pct.median(), 1)))
    print(pd.DataFrame(rows).to_string(index=False))

    print("\nNamed candidates:")
    pat = r"diaphanous|\bformin\b|\bformin-|dispatched|perlucin"
    cand = m[m.description.str.contains(pat, case=False, regex=True)].drop_duplicates("gene_id")
    print(cand[["region", "gene", "start_x", "description", "mantle_TPM", "other_TPM",
                "top_tissue", "mantle_pct"]].sort_values("start_x").to_string(index=False))

    for reg in ["W_LG17", "S_LG6"]:
        r = m[(m.region == reg) & (m.expressed == True)].sort_values("mantle_log2FC", ascending=False)
        print(f"\nMost mantle-biased genes in {reg}:")
        print(r[["gene", "start_x", "description", "mantle_TPM", "other_TPM", "mantle_pct"]]
              .head(15).to_string(index=False))

    # figure: genome distribution of mantle bias, with candidates marked
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.scatter(np.log10(ex.other_TPM + 1), np.log10(ex.mantle_TPM + 1), s=2, c="#cccccc", lw=0,
               label="all expressed genes")
    ax.scatter(np.log10(ctl.other_TPM + 1), np.log10(ctl.mantle_TPM + 1), s=12, c="#3a6ea5",
               label="shell-forming gene names (controls)")
    cd = cand[cand.expressed == True]
    ax.scatter(np.log10(cd.other_TPM + 1), np.log10(cd.mantle_TPM + 1), s=30, c="#c0392b",
               label="QTL candidates")
    short = lambda d: (d.replace("protein ", "").replace(" homolog", "").replace("-like", "")
                       .replace("dispatched 1", "dispatched"))
    for r in cd.itertuples():
        ax.annotate(f"{short(r.description)} ({r.region.split('_')[-1]})",(np.log10(r.other_TPM + 1), np.log10(r.mantle_TPM + 1)), fontsize=7,
                    xytext=(3, 3), textcoords="offset points")
    lim = max(ax.get_xlim()[1], ax.get_ylim()[1])
    ax.plot([0, lim], [0, lim], "--", color="grey", lw=0.8)
    ax.set_xlabel("log10(TPM + 1), mean of other tissues")
    ax.set_ylabel("log10(TPM + 1), mantle")
    ax.set_title("Littorina saxatilis tissue panel (PRJNA550990, via NCBI RefSeq annotation)")
    ax.legend(fontsize=8, loc="upper left")
    fig.tight_layout()
    fig.savefig(os.path.join(RES, "mantle_expression.png"), dpi=130)
    print("\nsaved mantle_expression.png")


if __name__ == "__main__":
    main()
