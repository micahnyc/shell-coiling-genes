"""
Is the mollusc-specific dispatched-like lineage expressed in shell-field cells or in
ciliated cells? Single-cell RNA-seq of the Dreissena rostriformis trochophore larva
(Salamanca-Diaz et al. 2022; GEO GSE192624, 10x Genomics, 13-15 h post fertilisation).
The authors describe shell-field (Hox1), prototroch (Caveolin) and other ciliated (FoxJ)
cell populations.

Steps
1. Translate the deposited gene models (GSM5753523_all_genes.fasta) to proteins and search
   them with DIAMOND (run_diamond.sh in data/dreissena_sc, WSL) using:
   - our dispatched proteins (data/gene_trees/dispatched.faa): each Dreissena gene takes the
     group (mollusc-specific lineage, DISP1/2, DISP3) of its best-scoring query;
   - human marker proteins: cilia (FOXJ1, tektins, radial spoke, axonemal dyneins, SPAG6,
     CFAP52), shell field (EN1/engrailed, HOXA1/Hox1, BMP2/4, GATA2/3, SMAD1, CHRDL1),
     prototroch (CAV1/2). Each marker takes its best Dreissena hit.
2. Cells with >= 300 UMIs; log CPM; 2,000 most variable genes; PCA (20 PCs); k-means (k=7,
   as in the paper; 20 restarts).
3. Per cluster: gene-set scores (mean z-scored log CPM) and the share of cells expressing
   each dispatched gene. Per cell: correlation of the main lineage gene with each score.

Usage:  python analysis/dreissena_singlecell.py   (after the DIAMOND searches)
Writes analysis/results/dreissena_singlecell_clusters.csv and dreissena_singlecell.png
"""

import gzip
import os

import numpy as np
import pandas as pd
import scipy.io
from scipy import stats
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..")
RES = os.path.join(HERE, "results")
D = os.path.join(ROOT, "data", "dreissena_sc")
MIN_UMI, K, NPC, NHVG, SEED = 300, 7, 20, 2000, 1


def dispatched_genes():
    h = pd.read_csv(os.path.join(D, "hits.tsv"), sep="\t", header=None,
                    names=["q", "s", "pid", "len", "e", "bits"])
    h = h[~h.q.str.startswith("marker")]
    a = pd.read_csv(os.path.join(RES, "gene_tree_dispatched_assignments.csv"))
    grp = dict(zip(a.protein, a.subfamily))
    h["group"] = h.q.str.split("|").str[2].map(grp).fillna(h.q.str.split("|").str[1])
    best = h.sort_values("bits", ascending=False).drop_duplicates("s")
    best = best[(best.bits >= 150) & (best.group != "PTCH1")]
    return dict(zip(best.s, best.group))


def marker_genes():
    h = pd.read_csv(os.path.join(D, "hits_markers.tsv"), sep="\t", header=None,
                    names=["q", "s", "pid", "len", "e", "bits"])
    best = h.sort_values("bits", ascending=False).drop_duplicates("q")
    best["set"] = best.q.str.split("|").str[0]
    best["marker"] = best.q.str.split("|").str[1]
    return best[["set", "marker", "s", "bits"]]


def kmeans(X, k, restarts, rng):
    best = None
    for _ in range(restarts):
        c = X[rng.choice(len(X), k, replace=False)]
        for _ in range(100):
            lab = ((X[:, None, :] - c[None]) ** 2).sum(-1).argmin(1)
            new = np.array([X[lab == j].mean(0) if (lab == j).any() else c[j] for j in range(k)])
            if np.allclose(new, c):
                break
            c = new
        inertia = ((X - c[lab]) ** 2).sum()
        if best is None or inertia < best[0]:
            best = (inertia, lab)
    return best[1]


def load_cells(verbose=True):
    """Filtered log-normalised cells, gene index, PCs, k-means labels and marker-set scores."""
    rng = np.random.default_rng(SEED)
    m = scipy.io.mmread(gzip.open(os.path.join(D, "GSM5753523_matrix.mtx.gz"))).tocsc()
    genes = pd.read_csv(os.path.join(D, "GSM5753523_features.tsv.gz")).gene_short_name.values
    umi = np.asarray(m.sum(0)).ravel()
    keep = umi >= MIN_UMI
    X = np.asarray(m[:, keep].todense(), dtype=float)  # genes x cells
    X = np.log1p(X / X.sum(0) * 1e4)
    gi = {g: i for i, g in enumerate(genes)}
    if verbose:
        print(f"{keep.sum()} cells with >= {MIN_UMI} UMIs (paper: 632)")

    expressed = (X > 0).sum(1) >= 5
    mu, var = X.mean(1), X.var(1)
    disp = np.where(expressed, var / (mu + 1e-9), -1)
    hv = np.argsort(disp)[-NHVG:]
    Z = X[hv].T
    Z = (Z - Z.mean(0)) / (Z.std(0) + 1e-9)
    U, S, _ = np.linalg.svd(Z, full_matrices=False)
    pcs = U[:, :NPC] * S[:NPC]
    lab = kmeans(pcs, K, 20, rng)

    zs = (X - X.mean(1, keepdims=True)) / (X.std(1, keepdims=True) + 1e-9)
    mk = marker_genes()
    scores = {}
    for s, g in mk.groupby("set"):
        idx = [gi[x] for x in g.s.unique() if x in gi and X[gi[x]].any()]
        if not idx:  # e.g. the Caveolin hits are not detected in any cell
            if verbose:
                print(f"{s}: no marker detected in these cells; skipped")
            continue
        scores[s] = zs[idx].mean(0)
        if verbose:
            print(f"{s}: {len(idx)} genes -> {', '.join(g.marker)}")
    return X, gi, pcs, lab, scores


def main():
    X, gi, pcs, lab, scores = load_cells()
    dg = dispatched_genes()

    rows = []
    for j in range(K):
        c = lab == j
        r = dict(cluster=j, cells=int(c.sum()))
        for s in scores:
            r[f"{s}_score"] = scores[s][c].mean()
        for g, grp in dg.items():
            if g in gi:
                r[f"{g} ({grp.split(' (')[0]})"] = (X[gi[g], c] > 0).mean() * 100
        rows.append(r)
    cl = pd.DataFrame(rows)
    cl.to_csv(os.path.join(RES, "dreissena_singlecell_clusters.csv"), index=False, float_format="%.3f")
    pd.set_option("display.width", 250)
    print("\nClusters: gene-set scores and % of cells expressing each dispatched gene")
    print(cl.round(2).to_string(index=False))

    main_gene = max((g for g, grp in dg.items() if grp.startswith("mollusc") and g in gi),
                    key=lambda g: (X[gi[g]] > 0).sum())
    y = X[gi[main_gene]]
    print(f"\nPer cell, {main_gene} ({(y > 0).sum()} expressing cells) vs gene-set scores (Spearman):")
    for s in scores:
        rho, p = stats.spearmanr(y, scores[s])
        print(f"  {s:12s} rho = {rho:+.3f}  p = {p:.2g}")

    # figure: PCA coloured by cluster, cilia score, shell-field score, lineage gene
    fig, axes = plt.subplots(1, 4, figsize=(17, 4.2))
    panels = [("k-means cluster", lab, "tab10"), ("cilia score", scores["cilia"], "viridis"),
              ("shell-field score", scores["shellfield"], "viridis"),
              (f"{main_gene} (mollusc-specific dispatched-like)", y, "magma_r")]
    for ax, (t, v, cmap) in zip(axes, panels):
        order = np.argsort(v) if t != "k-means cluster" else np.arange(len(v))
        sc = ax.scatter(pcs[order, 0], pcs[order, 1], c=v[order], cmap=cmap, s=7, lw=0)
        ax.set_title(t, fontsize=9)
        ax.set_xticks([]); ax.set_yticks([])
        ax.set_xlabel("PC1"); ax.set_ylabel("PC2")
        if t == "k-means cluster":
            for j in range(K):
                ax.text(*np.median(pcs[lab == j, :2], 0), str(j), fontsize=10, weight="bold")
        else:
            fig.colorbar(sc, ax=ax, shrink=0.7)
    fig.suptitle("Dreissena rostriformis trochophore, single cells (GSE192624)", fontsize=10)
    fig.tight_layout()
    fig.savefig(os.path.join(RES, "dreissena_singlecell.png"), dpi=120)


if __name__ == "__main__":
    main()
