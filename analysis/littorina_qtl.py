"""
QTL scan: which parts of the Littorina saxatilis genome go with each Raup-type
shell parameter?

Data: Koch et al. 2021, Evolution Letters, Dryad doi:10.5061/dryad.zgmsbccb4.
F2 snails from Crab x Wave ecotype crosses, 13 families, ~19,800 SNP markers
on a linkage map, and shell shape measured with the growth model of
Larsson et al. 2020 (J. R. Soc. Interface 17:20190721).

Translation from the Larsson model to Raup's parameters
-------------------------------------------------------
  widthGrowth      = ln(g_w), g_w = width growth per radian
                     -> W = exp(2*pi*g_w)   (Raup's whorl expansion per turn)
  heightGrowth     = ln(g_h), g_h = height growth per radian
                     -> W_height = exp(2*pi*g_h); spire shape changes with
                        growth when g_h != g_w (T is not constant)
  aperturePosition = r0, distance of the opening's centre from the axis
  apertureSize     = a0, radius of the opening (both relative to shell length)
                     -> D_like = r0 / (r0 + a0)   (Raup's D measures the same
                        thing: how far the tube sits from the axis)
  apertureShape    = c, elongation of the opening -> S

Method (marker regression, as in Haley-Knott scans at the markers)
-----------------------------------------------------------------
For each trait and marker: y ~ covariates + additive + dominance,
LOD = n/2 * log10(RSS_null / RSS_marker). Covariates: family, batch, sex,
and (for shape traits, as in the original paper) shell length.
Genome-wide 5% thresholds from 1000 permutations of the covariate-adjusted
trait. Peaks are reported with 1.5-LOD support intervals.

Usage:  python analysis/littorina_qtl.py
Writes results to analysis/results/.
"""

import os

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data", "littorina_koch2021")
OUT = os.path.join(HERE, "results")
N_PERM = 1000
SEED = 1

# QTL-map linkage group -> Crab-map linkage group (from the authors' scripts)
CRAB_LG = {1: 12, 2: 2, 3: 1, 4: 3, 5: 17, 6: 6, 7: 4, 8: 7, 9: 9, 10: 14,
           11: 10, 12: 5, 13: 13, 14: 15, 15: 11, 16: 8, 17: 16}

TRAITS = {  # name: (description, adjust for shell length?)
    "W": ("Raup W: whorl expansion per turn (width)", True),
    "W_height": ("height expansion per turn", True),
    "D_like": ("Raup D: tube's distance from axis, r0/(r0+a0)", True),
    "S": ("Raup S: opening elongation (c)", True),
    "apertureSize": ("opening size a0 (relative)", True),
    "shellLength": ("shell length (size, for comparison)", False),
}


def load():
    p = pd.read_csv(os.path.join(DATA, "Pheno.QTL.csv"))

    ped = pd.read_csv(os.path.join(DATA, "Pedigree.txt"), sep="\t", header=None)
    fam = pd.Series(ped.iloc[0, 2:].values, index=ped.iloc[1, 2:].values)
    fam = fam[~fam.index.duplicated(keep=False)]  # parents appear in several families
    p["family"] = p["ind"].map(fam)

    p["W"] = np.exp(2 * np.pi * np.exp(p["widthGrowth"]))
    p["W_height"] = np.exp(2 * np.pi * np.exp(p["heightGrowth"]))
    p["D_like"] = p["aperturePosition"] / (p["aperturePosition"] + p["apertureSize"])
    p["S"] = p["apertureShape"]

    g = pd.read_csv(os.path.join(DATA, "Geno.QTL.csv"), header=None, low_memory=False)
    markers = pd.DataFrame({
        "marker": g.iloc[0, 1:].astype(int).values,
        "qtl_lg": g.iloc[1, 1:].astype(int).values,
        "cM": g.iloc[2, 1:].astype(float).values,
    })
    markers["LG"] = markers["qtl_lg"].map(CRAB_LG)
    geno = g.iloc[3:, 1:].values
    gid = g.iloc[3:, 0].astype(int).values
    assert (gid == p["id"].values).all(), "genotype rows must follow Pheno.QTL.csv order"
    return p, markers, geno


def design(p, use_length):
    cov = pd.get_dummies(p[["family"]].astype(str), drop_first=True).astype(float)
    cov = cov.join(pd.get_dummies(p["batch"].astype(int).astype(str), prefix="batch",
                                  drop_first=True).astype(float))
    cov["sex"] = p["sex"].astype(float)
    if use_length:
        cov["shellLength"] = p["shellLength"]
    cov.insert(0, "intercept", 1.0)
    return cov.values


def scan(y, C, Xa, Xd, n_perm=0, rng=None):
    """LOD for each marker; optionally genome-wide max-LOD permutation distribution."""
    n = len(y)
    Q, _ = np.linalg.qr(C)
    resid = lambda M: M - Q @ (Q.T @ M)
    ry, Ra, Rd = resid(y[:, None])[:, 0], resid(Xa), resid(Xd)
    saa, sdd, sad = (Ra * Ra).sum(0), (Rd * Rd).sum(0), (Ra * Rd).sum(0)
    det = saa * sdd - sad ** 2
    det[det <= 1e-12] = np.nan  # monomorphic marker in this subset

    def lod_for(Y):  # Y: n x k residualised traits
        say, sdy = Ra.T @ Y, Rd.T @ Y
        expl = (sdd[:, None] * say ** 2 - 2 * sad[:, None] * say * sdy
                + saa[:, None] * sdy ** 2) / det[:, None]
        rss0 = (Y ** 2).sum(0)
        return n / 2 * np.log10(rss0 / (rss0 - expl))

    lod = lod_for(ry[:, None])[:, 0]
    # additive effect at each marker (in trait units)
    say, sdy = Ra.T @ ry, Rd.T @ ry
    add = (sdd * say - sad * sdy) / det
    perm = None
    if n_perm:
        Yp = np.column_stack([rng.permutation(ry) for _ in range(n_perm)])
        perm = np.nanmax(lod_for(Yp), axis=0)
    return lod, add, perm


def peaks(lod, add, markers, thr, sd):
    rows = []
    for lg, idx in markers.groupby("LG").groups.items():
        idx = np.asarray(idx)
        l = lod[idx]
        if np.all(np.isnan(l)):
            continue
        j = idx[np.nanargmax(l)]
        if lod[j] < thr:
            continue
        within = idx[l >= lod[j] - 1.5]
        rows.append(dict(LG=lg, peak_cM=markers.cM[j], LOD=lod[j],
                         interval_cM=f"{markers.cM[within].min():.1f}-{markers.cM[within].max():.1f}",
                         additive_effect_SD=add[j] / sd))
    return rows


def main():
    os.makedirs(OUT, exist_ok=True)
    rng = np.random.default_rng(SEED)
    p, markers, geno = load()

    # sort markers along the genome in Crab-map order
    order = markers.sort_values(["LG", "cM"]).index.values
    markers = markers.loc[order].reset_index(drop=True)
    geno = geno[:, order]
    Xa_all = np.select([geno == "A", geno == "H", geno == "B"], [-1.0, 0.0, 1.0], np.nan)
    Xd_all = (geno == "H").astype(float)

    results, all_peaks, lods = {}, [], pd.DataFrame(markers[["marker", "LG", "cM"]])
    for t, (desc, use_len) in TRAITS.items():
        need = [t, "family", "batch", "sex"] + (["shellLength"] if use_len else [])
        keep = p[need].notna().all(axis=1).values
        sub = p[keep]
        y = sub[t].values.astype(float)
        C = design(sub, use_len)
        lod, add, perm = scan(y, C, Xa_all[keep], Xd_all[keep], N_PERM, rng)
        thr = float(np.quantile(perm, 0.95))
        lods[t] = lod
        results[t] = dict(n=int(keep.sum()), threshold=thr)
        for r in peaks(lod, add, markers, thr, y.std()):
            r.update(trait=t, n=int(keep.sum()), threshold=round(thr, 2),
                     pve_percent=100 * (1 - 10 ** (-2 * r["LOD"] / keep.sum())))
            all_peaks.append(r)
        print(f"{t:13s} n={keep.sum()}  5% threshold LOD={thr:.2f}  max LOD={np.nanmax(lod):.2f}")

    pk = pd.DataFrame(all_peaks)
    if len(pk):
        pk = pk[["trait", "LG", "peak_cM", "interval_cM", "LOD", "threshold",
                 "pve_percent", "additive_effect_SD", "n"]]
        pk.to_csv(os.path.join(OUT, "qtl_peaks.csv"), index=False, float_format="%.3f")
        print("\nSignificant QTL peaks (best marker per linkage group):")
        print(pk.to_string(index=False, float_format=lambda x: f"{x:.2f}"))
    lods.to_csv(os.path.join(OUT, "lod_scores.csv.gz"), index=False, float_format="%.3f")

    # trait correlations (do the parameters vary together?)
    corr = p[list(TRAITS)].corr()
    corr.to_csv(os.path.join(OUT, "trait_correlations.csv"), float_format="%.3f")
    print("\nTrait correlations:\n", corr.round(2).to_string())

    plot(lods, markers, results)


def plot(lods, markers, results):
    traits = list(TRAITS)
    lg_len = markers.groupby("LG").cM.max()
    offset = (lg_len + 10).cumsum().shift(fill_value=0)
    x = markers.cM + markers.LG.map(offset)
    fig, axes = plt.subplots(len(traits), 1, figsize=(14, 2.1 * len(traits)), sharex=True)
    for ax, t in zip(axes, traits):
        for lg in lg_len.index:
            m = (markers.LG == lg).values
            ax.plot(x[m], lods[t][m], lw=0.8, color="#8a3b12" if lg % 2 else "#d08a3c")
        ax.axhline(results[t]["threshold"], ls="--", color="grey", lw=0.8)
        ax.set_ylabel("LOD")
        ax.set_title(f"{t}: {TRAITS[t][0]}", fontsize=9, loc="left")
        ax.spines[["top", "right"]].set_visible(False)
    axes[-1].set_xticks(offset + lg_len / 2)
    axes[-1].set_xticklabels([f"LG{lg}" for lg in lg_len.index])
    fig.suptitle("Littorina saxatilis F2 (Koch et al. 2021): genome scans for Raup-type shell "
                 "parameters\n(dashed = genome-wide 5% threshold, 1000 permutations)", fontsize=11)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "qtl_scan.png"), dpi=130)
    print(f"\nsaved {os.path.join(OUT, 'qtl_scan.png')}")


if __name__ == "__main__":
    main()
