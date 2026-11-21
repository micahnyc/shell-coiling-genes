"""
Which periwinkle genes correspond to proteins found in mollusc shells, and are the QTL
inversions enriched for them?

Shell-matrix proteomes (proteins identified by mass spectrometry in shell):
  Lottia gigantea, Mann et al. 2012 (consensus lists, Additional files 1-2; 600 proteins)
  Crepidula fornicata, Batzel et al. 2022 (185); Crepidula atrasolea, 2024 (18)
  Haliotis laevigata, Mann et al. 2018 (80); Pinctada fucata, Liu et al. 2015 (43)
  Cepaea nemoralis, Mann & Jackson 2014 (59 contigs, translated)
Sequences were gathered from GenBank / the papers' supplements into data/shell_proteomes/.
Search: DIAMOND (analysis/orthologs/run_shellprot.sh), e <= 1e-10; each shell protein is
assigned to its best periwinkle hit (gene level). Some listed proteins are cellular
carry-over (actin, tubulin, histones), so the gene set is "shell-proteome homologs", not a
clean list of secreted shell proteins.

Enrichment: hypergeometric test of homolog genes among protein-coding genes in
  LG6 inversion  NC_090250.1 1.35-17.72 Mb
  LG17 inversion NC_090261.1 63.0-70.3 Mb (Koch 2021 Table S1 bounds via our marker placement)
Also lists which homologs are mantle-specific (top 5% mantle bias).

Usage:  python analysis/shell_proteome_homologs.py
Writes analysis/results/shell_proteome_homologs.csv
"""

import os

import pandas as pd
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..")
RES = os.path.join(HERE, "results")
D = os.path.join(ROOT, "data", "shell_proteomes")
SETS = ["lottia_consensus", "crepidula_fornicata", "crepidula_atrasolea", "haliotis_laevigata",
        "pinctada_fucata", "cepaea_orfs"]
REGIONS = {"LG6 inversion": ("NC_090250.1", 1.35e6, 17.72e6), "LG17 inversion": ("NC_090261.1", 63.0e6, 70.3e6)}


def main():
    ft = pd.read_csv(os.path.join(ROOT, "data", "orthologs", "GCF_037325665.1_US_GU_Lsax_2.0_feature_table.txt.gz"),
                     sep="\t", low_memory=False)
    cds = ft[ft["# feature"] == "CDS"].dropna(subset=["product_accession"])
    p2g = dict(zip(cds.product_accession, cds.GeneID.astype("int64")))
    genes = ft[(ft["# feature"] == "gene") & (ft["class"] == "protein_coding")][
        ["GeneID", "genomic_accession", "start", "end", "name", "symbol"]].copy()
    genes["GeneID"] = genes.GeneID.astype("int64")
    names = dict(zip(cds.GeneID.astype("int64"), cds.name))

    hits = []
    for s in SETS:
        h = pd.read_csv(os.path.join(D, f"{s}.vs_lit.tsv"), sep="\t", header=None,
                        names=["q", "s", "pid", "len", "e", "bits"])
        h = h.sort_values("bits", ascending=False).drop_duplicates("q")
        h["gene"] = h.s.map(p2g)
        h["set"] = s
        hits.append(h)
    h = pd.concat(hits)
    hom = h.groupby("gene").agg(n_shell_proteins=("q", "size"), sets=("set", lambda x: ",".join(sorted(set(x)))),
                                best_pid=("pid", "max")).reset_index()
    hom["name"] = hom.gene.map(names)

    mx = pd.read_csv(os.path.join(RES, "mantle_expression_genome.csv.gz"), index_col=0)
    hom["mantle_pct"] = hom.gene.map(mx.mantle_pct)
    hom = hom.merge(genes[["GeneID", "genomic_accession", "start"]], left_on="gene", right_on="GeneID", how="left")
    hom["region"] = ""
    N, K = len(genes), hom.gene.nunique()
    print(f"{K} periwinkle genes are best hits of shell-proteome proteins "
          f"({len(h)} shell proteins with hits); {N} protein-coding genes")
    for reg, (chrom, a, b) in REGIONS.items():
        inreg = genes[(genes.genomic_accession == chrom) & (genes.start.between(a, b))]
        k = hom.gene.isin(inreg.GeneID).sum()
        n = len(inreg)
        p = stats.hypergeom.sf(k - 1, N, K, n)
        print(f"{reg}: {k} of {n} genes are shell-proteome homologs (expected {n * K / N:.1f}); "
              f"enrichment p = {p:.3g}")
        hom.loc[hom.gene.isin(inreg.GeneID), "region"] = reg
    hom = hom.sort_values(["region", "start"], ascending=[False, True])
    hom.drop(columns=["GeneID"]).to_csv(os.path.join(RES, "shell_proteome_homologs.csv"), index=False, float_format="%.1f")
    pd.set_option("display.width", 220)
    print(hom[hom.region != ""][["region", "start", "name", "sets", "n_shell_proteins", "best_pid", "mantle_pct"]]
          .round(1).to_string(index=False))
    top = hom[hom.mantle_pct >= 95].sort_values("mantle_pct", ascending=False)
    print(f"\nMantle-specific (top 5%) shell-proteome homologs genome-wide: {len(top)} of {hom.mantle_pct.notna().sum()}")
    print(top[["genomic_accession", "start", "name", "sets", "mantle_pct"]].head(25).to_string(index=False))


if __name__ == "__main__":
    main()
