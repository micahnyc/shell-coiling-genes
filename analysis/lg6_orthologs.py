"""
Do the LG6 candidate genes behave like shell genes in other molluscs?

Candidates: the 14 genes in the LG6 shape region that are among the periwinkle genome's 5%
most mantle-biased (mantle_expression.py), 8 inside the inversion and 6 outside, plus the
LG6 diaphanous gene. Controls: 9 known shell-forming genes that are strongly mantle-biased
in periwinkle (carbonic anhydrases, tyrosinase, chitin synthase, nacrein, perlucin).

Orthologs: reciprocal best hits with DIAMOND (analysis/orthologs/run_rbh.sh) in
  Crassostrea virginica (oyster)   adult tissues + 1-day larvae (NCBI RefSeq TPMs)
  Mytilus galloprovincialis (mussel) 0-72 h development (NCBI RefSeq counts, CPM)
  Dreissena rostriformis           trochophore single cells (GSE192624)
If the best hit is not reciprocal, it is still reported but flagged, since gene families
such as carbonic anhydrases have many close paralogs.

Shell-like behaviour, per species:
  oyster  mantle_pct      percentile of mantle bias (mantle vs mean of other adult tissues)
          larva_pct       percentile of expression level in 1-day larvae (first shell)
  mussel  sf_fold         log2 CPM fold, shell-field window (16-28 h) vs 0-8 h
          veliger_fold    log2 CPM fold, D-veliger (44-72 h) vs 0-8 h
  Dreissena rho_shell / rho_cilia   per-cell Spearman with the shell-field and cilia scores

Usage:  python analysis/lg6_orthologs.py
Writes analysis/results/lg6_orthologs.csv and lg6_orthologs.png
"""

import os
import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import dev_expression as de  # noqa: E402
import dev_timecourse as dt  # noqa: E402
import dreissena_singlecell as ds  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..")
RES = os.path.join(HERE, "results")
O = os.path.join(ROOT, "data", "orthologs")
COLS = ["qseqid", "sseqid", "pident", "length", "evalue", "bitscore"]


def rbh(sp):
    f = pd.read_csv(os.path.join(O, f"fwd_{sp}.tsv"), sep="\t", names=COLS)
    r = pd.read_csv(os.path.join(O, f"rev_{sp}.tsv"), sep="\t", names=COLS)
    fb = f.sort_values("bitscore", ascending=False).drop_duplicates("qseqid").set_index("qseqid")
    rb = r.sort_values("bitscore", ascending=False).drop_duplicates("qseqid").set_index("qseqid")
    out = {}
    for q, row in fb.iterrows():
        back = rb.sseqid.get(row.sseqid)
        out[q] = (row.sseqid, row.pident, back == q)
    return out


def protein_to_gene(feature_table):
    ft = pd.read_csv(feature_table, sep="\t", usecols=["# feature", "product_accession", "GeneID"], low_memory=False)
    ft = ft[(ft["# feature"] == "CDS") & ft.product_accession.notna()]
    return dict(zip(ft.product_accession, ft.GeneID.astype("int64").astype(str)))


def oyster_table():
    asm, rel = de.SPECIES["Crassostrea virginica"]
    d = os.path.join(de.DATA, rel)
    runs = pd.read_csv(os.path.join(d, "runs.txt"), sep="\t")
    tpm = pd.read_csv(os.path.join(d, "tpm.txt.gz"), sep="\t")
    tpm.index = tpm.GeneID.str.replace("GeneID:", "", regex=False)
    runs["label"] = runs["Sample Name"].map(de.tissue_label)
    cols = [r for r in runs["SRA run accession"] if r in tpm.columns]
    lab = runs.set_index("SRA run accession").label[cols]
    by = tpm[cols].T.groupby(lab.values).mean().T
    adults = ["gill", "gonad", "hemocytes", "mantle"]
    bias = np.log2((by["mantle"] + 1) / (by[[a for a in adults if a != "mantle"]].mean(1) + 1))
    expressed = by.max(1) >= 1
    out = pd.DataFrame(index=by.index)
    out["mantle_TPM"] = by["mantle"]
    out["mantle_pct"] = bias[expressed].rank(pct=True) * 100
    lv = [c for c in by.columns if c.startswith("larvae")][0]
    out["larva_TPM"] = by[lv]
    out["larva_pct"] = by[lv][expressed].rank(pct=True) * 100
    return out, os.path.join(d, "feature_table.txt.gz")


def mussel_table():
    by, _, _ = dt.mytilus()
    early = by[[c for c in by.columns if c <= 8]].mean(1)
    sf = by[[c for c in by.columns if 16 <= c <= 28]].mean(1)
    vel = by[[c for c in by.columns if c >= 44]].mean(1)
    out = pd.DataFrame(index=by.index)
    out["early_CPM"], out["sf_CPM"], out["veliger_CPM"] = early, sf, vel
    out["sf_fold"] = np.log2((sf + 1) / (early + 1))
    out["veliger_fold"] = np.log2((vel + 1) / (early + 1))
    return out, os.path.join(dt.DATA, "GCF_965363235.1_xbMytGall1.hap1.1_feature_table.txt.gz")


def main():
    q = pd.read_csv(os.path.join(O, "littorina_queries.csv"))
    oy, oy_ft = oyster_table()
    mu, mu_ft = mussel_table()
    X, gi, pcs, lab, scores = ds.load_cells(verbose=False)
    maps = {"cvi": protein_to_gene(oy_ft), "myt": protein_to_gene(mu_ft)}
    hits = {sp: rbh(sp) for sp in ("cvi", "myt", "dre")}

    rows = []
    for r in q.itertuples():
        row = dict(set=r.set, description=r.description, littorina_gene=r.gene_id)
        for sp, label in (("cvi", "oyster"), ("myt", "mussel"), ("dre", "dreissena")):
            h = hits[sp].get(r.protein)
            if h is None:
                row[f"{label}_ortholog"] = "none"
                continue
            prot, pid, recip = h
            row[f"{label}_ortholog"] = f"{prot} ({pid:.0f}%{'' if recip else ', not reciprocal'})"
            row[f"{label}_reciprocal"] = recip
            if sp == "cvi":
                g = maps["cvi"].get(prot)
                if g in oy.index:
                    row.update(oyster_mantle_pct=oy.at[g, "mantle_pct"], oyster_larva_pct=oy.at[g, "larva_pct"],
                               oyster_mantle_TPM=oy.at[g, "mantle_TPM"])
            elif sp == "myt":
                g = maps["myt"].get(prot)
                if g in mu.index:
                    row.update(mussel_sf_fold=mu.at[g, "sf_fold"], mussel_veliger_fold=mu.at[g, "veliger_fold"],
                               mussel_veliger_CPM=mu.at[g, "veliger_CPM"])
            else:
                if prot in gi:
                    y = X[gi[prot]]
                    row["dreissena_cells"] = int((y > 0).sum())
                    if (y > 0).sum() >= 10:
                        row["dreissena_rho_shell"] = stats.spearmanr(y, scores["shellfield"])[0]
                        row["dreissena_rho_cilia"] = stats.spearmanr(y, scores["cilia"])[0]
        rows.append(row)
    t = pd.DataFrame(rows)
    t.to_csv(os.path.join(RES, "lg6_orthologs.csv"), index=False, float_format="%.3f")
    pd.set_option("display.width", 260)
    show = ["set", "description", "oyster_mantle_pct", "oyster_larva_pct", "mussel_sf_fold", "mussel_veliger_fold",
            "dreissena_cells", "dreissena_rho_shell", "dreissena_rho_cilia", "oyster_reciprocal", "mussel_reciprocal",
            "dreissena_reciprocal"]
    print(t[[c for c in show if c in t]].round(2).to_string(index=False))

    # figure: metric heatmap, rows grouped by set
    metrics = [("oyster_mantle_pct", "oyster: mantle bias\n(percentile)", 0, 100),
               ("oyster_larva_pct", "oyster: level in\n1-day larvae (pct)", 0, 100),
               ("mussel_sf_fold", "mussel: shell field vs\nearly (log2 fold)", -4, 4),
               ("mussel_veliger_fold", "mussel: D-veliger vs\nearly (log2 fold)", -4, 4),
               ("dreissena_rho_shell", "Dreissena: tracks\nshell-field cells (rho)", -0.3, 0.3),
               ("dreissena_rho_cilia", "Dreissena: tracks\nciliated cells (rho)", -0.3, 0.3)]
    order = ["control: known shell gene", "LG6 inversion", "LG6 inversion (diaphanous)", "LG6 S region, outside inversion"]
    t["o"] = t.set.map({s: i for i, s in enumerate(order)})
    t = t.sort_values(["o", "description"])
    fig, ax = plt.subplots(figsize=(12, 0.36 * len(t) + 2))
    for j, (c, lab_, lo, hi) in enumerate(metrics):
        v = t[c].values.astype(float) if c in t else np.full(len(t), np.nan)
        norm = (np.clip(v, lo, hi) - lo) / (hi - lo)
        for i, (x, n) in enumerate(zip(v, norm)):
            colour = plt.cm.RdYlBu_r(n) if not np.isnan(x) else (0.93, 0.93, 0.93, 1)
            ax.add_patch(plt.Rectangle((j, i), 1, 1, color=colour))
            ax.text(j + 0.5, i + 0.5, "" if np.isnan(x) else (f"{x:.0f}" if "pct" in c else f"{x:+.2f}"),
                    ha="center", va="center", fontsize=7)
    ax.set_xlim(0, len(metrics))
    ax.set_ylim(len(t), 0)
    ax.set_xticks(np.arange(len(metrics)) + 0.5)
    ax.set_xticklabels([m[1] for m in metrics], fontsize=8)
    ax.xaxis.tick_top()
    ax.set_yticks(np.arange(len(t)) + 0.5)
    ax.set_yticklabels([f"{s.replace('control: known shell gene', 'control')[:22]} | {d[:45]}"
                        for s, d in zip(t.set, t.description)], fontsize=7.5)
    for i, (a, b) in enumerate(zip(t.o.values[:-1], t.o.values[1:])):
        if a != b:
            ax.axhline(i + 1, color="black", lw=1)
    ax.set_title("LG6 candidates vs known shell genes in three other molluscs (grey = no ortholog or too few cells)",
                 fontsize=9, pad=40)
    fig.tight_layout()
    fig.savefig(os.path.join(RES, "lg6_orthologs.png"), dpi=120)


if __name__ == "__main__":
    main()
