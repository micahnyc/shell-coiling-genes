"""
Step 1 (wild snails): do inversion karyotypes predict Raup-type shell shape in wild
Littorina saxatilis, beyond the Crab/Wave habitat each snail lives in?

Data: Koch et al. 2022, Evolution 76:2332, Dryad doi:10.5061/dryad.m905qfv4b
(Data_Transects.txt): ~1,100 snails from three Swedish Crab-Wave transects
(islands CZA, CZB, CZD) with Larsson growth-model shape parameters and karyotypes
(0/1/2) for 16 inversions, e.g. LGC6.1 (LG6) and LGC17.1 (LG17).

Traits (same Raup translation as littorina_qtl.py; here growth rates are per radian,
not log): W = exp(2*pi*Width-Growth), W_height = exp(2*pi*Height-Growth),
D_like = position/(position+size), S = Aperture-Shape, plus aperture size and shell length.

Model per trait and inversion (ordinary least squares):
  trait ~ island + sex + habitatPC1 + habitatPC1:island + [shell length] + karyotype (0/1/2)
habitatPC1 is the habitat axis (boulder/Crab vs bedrock/Wave). Inversions and shape both
track habitat, so without it any habitat-linked inversion would look like a shape locus.
Reported: karyotype effect per inversion copy (in SD of the trait), its t and p value,
and Bonferroni significance over 16 inversions x 6 traits.
A second model fits all 16 inversions at once (snails with complete karyotypes).

Usage:  python analysis/wild_inversions.py
Writes analysis/results/wild_inversions.csv and wild_inversions.png
"""

import os

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..")
RES = os.path.join(HERE, "results")
DATA = os.path.join(ROOT, "data", "koch2022", "Data_Transects.txt")

TRAITS = {  # name: (label, adjust for shell length?)
    "W": ("W, coil expansion", True), "W_height": ("height expansion", True),
    "D_like": ("D, distance from axis", True), "S": ("S, opening shape", True),
    "apertureSize": ("opening size", True), "shellLength": ("shell length", False),
}


def load():
    d = pd.read_csv(DATA, sep="\t")
    d.columns = [c.strip() for c in d.columns]
    d["W"] = np.exp(2 * np.pi * d["Width-Growth"])
    d["W_height"] = np.exp(2 * np.pi * d["Height-Growth"])
    d["D_like"] = d["Aperture-Position"] / (d["Aperture-Position"] + d["Aperture-Size"])
    d["S"] = d["Aperture-Shape"]
    d["apertureSize"] = d["Aperture-Size"]
    d["shellLength"] = d["Shell.Length"]
    invs = [c for c in d.columns if c.startswith("LGC")]
    return d, invs


def covariates(d, use_len):
    X = pd.get_dummies(d.location, prefix="isl", drop_first=True).astype(float)
    X["male"] = (d.sex == "m").astype(float)
    X["hab"] = d.habitatPC1
    for c in [c for c in X.columns if c.startswith("isl_")]:
        X[f"hab_{c}"] = X.hab * X[c]
    if use_len:
        X["len"] = d.shellLength
    X.insert(0, "const", 1.0)
    return X


def ols(y, X):
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    resid = y - X @ beta
    dof = len(y) - X.shape[1]
    s2 = resid @ resid / dof
    cov = s2 * np.linalg.inv(X.T @ X)
    se = np.sqrt(np.diag(cov))
    t = beta / se
    return beta, se, t, 2 * stats.t.sf(np.abs(t), dof), dof


def main():
    d, invs = load()
    print(f"{len(d)} snails; islands {d.location.value_counts().to_dict()}")
    rows = []
    for tr, (lab, use_len) in TRAITS.items():
        for inv in invs:
            need = [tr, inv, "location", "sex", "habitatPC1"] + (["shellLength"] if use_len else [])
            s = d[need].dropna()
            s = s[s.sex.isin(["m", "f"])]
            X = covariates(s, use_len)
            X["k"] = s[inv].astype(float)
            y = s[tr].values.astype(float)
            b, se, t, p, dof = ols(y, X.values.astype(float))
            j = list(X.columns).index("k")
            freq = s[inv].mean() / 2
            # same model without habitat, to show how much habitat confounds
            Xn = X.drop(columns=[c for c in X.columns if c.startswith("hab")])
            bn, _, tn, pn, _ = ols(y, Xn.values.astype(float))
            jn = list(Xn.columns).index("k")
            rows.append(dict(trait=tr, inversion=inv, n=len(s), arrangement_freq=round(freq, 3),
                             effect_SD=b[j] / y.std(), t=t[j], p=p[j],
                             effect_SD_no_habitat=bn[jn] / y.std(), p_no_habitat=pn[jn]))
    r = pd.DataFrame(rows)
    ntest = len(r)
    r["bonferroni_sig"] = r.p < 0.05 / ntest
    r.to_csv(os.path.join(RES, "wild_inversions.csv"), index=False, float_format="%.4g")
    pd.set_option("display.width", 200)
    print(f"\nSignificant after Bonferroni ({ntest} tests, p < {0.05 / ntest:.1e}), habitat-adjusted:")
    print(r[r.bonferroni_sig].sort_values("p")[["trait", "inversion", "n", "effect_SD", "t", "p",
                                                  "effect_SD_no_habitat", "p_no_habitat"]].to_string(index=False))
    print("\nLG6 and LG17 inversions, all traits:")
    print(r[r.inversion.isin(["LGC6.1", "LGC17.1"])][["trait", "inversion", "n", "arrangement_freq", "effect_SD",
                                                      "p", "effect_SD_no_habitat", "p_no_habitat"]].to_string(index=False))

    # joint model: all inversions together
    print("\nJoint model (all 16 inversions, complete cases):")
    jrows = []
    for tr, (lab, use_len) in TRAITS.items():
        s = d[list(dict.fromkeys([tr, "location", "sex", "habitatPC1", "shellLength"] + invs))].dropna()
        s = s[s.sex.isin(["m", "f"])]
        X = covariates(s, use_len)
        for inv in invs:
            X[inv] = s[inv].astype(float)
        y = s[tr].values.astype(float)
        b, se, t, p, dof = ols(y, X.values.astype(float))
        for inv in invs:
            j = list(X.columns).index(inv)
            jrows.append(dict(trait=tr, inversion=inv, n=len(s), effect_SD=b[j] / y.std(), p=p[j]))
    jr = pd.DataFrame(jrows)
    jr.to_csv(os.path.join(RES, "wild_inversions_joint.csv"), index=False, float_format="%.4g")
    print(jr[jr.p < 0.05 / ntest].sort_values("p").to_string(index=False))

    # figure: -log10 p heatmap, habitat-adjusted single-inversion models
    m = r.pivot(index="inversion", columns="trait", values="p")[list(TRAITS)]
    e = r.pivot(index="inversion", columns="trait", values="effect_SD")[list(TRAITS)]
    m = m.loc[sorted(m.index, key=lambda s: [int(x) for x in s[3:].split(".")])]
    e = e.loc[m.index]
    fig, ax = plt.subplots(figsize=(8, 7))
    im = ax.imshow(-np.log10(m.values), cmap="YlOrBr", aspect="auto", vmin=0, vmax=10)
    for i in range(m.shape[0]):
        for j in range(m.shape[1]):
            star = "*" if m.values[i, j] < 0.05 / ntest else ""
            ax.text(j, i, f"{e.values[i, j]:+.2f}{star}", ha="center", va="center", fontsize=7)
    ax.set_xticks(range(m.shape[1]))
    ax.set_xticklabels([TRAITS[c][0] for c in m.columns], rotation=25, ha="right", fontsize=8)
    ax.set_yticks(range(m.shape[0]))
    ax.set_yticklabels(m.index, fontsize=8)
    ax.set_title("Wild snails (Koch et al. 2022): inversion karyotype vs shell traits\n"
                 "colour = -log10 p, text = effect per copy (SD), * = Bonferroni; habitat-adjusted", fontsize=9)
    fig.colorbar(im, ax=ax, shrink=0.7, label="-log10 p (capped at 10)")
    fig.tight_layout()
    fig.savefig(os.path.join(RES, "wild_inversions.png"), dpi=120)


if __name__ == "__main__":
    main()
