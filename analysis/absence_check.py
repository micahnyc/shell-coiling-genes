"""
Are the absences of the mollusc-specific dispatched-like lineage real, or annotation gaps?

All dispatched-family proteins from the gene tree (labelled LIN = mollusc-specific lineage,
DISP12 = DISP1/2, DISP3, PTCH = patched outgroup) were aligned with miniprot directly to
the genomes of the five species whose annotation lacks the lineage, plus a positive
control (periwinkle LG17, which has 6 lineage genes and DISP1 annotated):
  Physella acuta GCF_028476545.2, Aplysia californica GCF_000002075.1,
  Berghia stephanieae GCF_034508935.2, Octopus bimaculoides GCF_046866305.1,
  Gigantopelta aegis GCF_016097555.1   (analysis/orthologs/run_absence_check.sh)
Nautilus (unannotated genome) is searched the same way by run_nautilus.sh.
Settings: miniprot -N 20 --outs=0.5. On the control, every lineage region and DISP1 is
found and correctly labelled; tandem copies a few kb apart merge into one locus, so this
counts presence reliably but undercounts copies.
Each genomic locus is labelled by the query group with the highest alignment score,
and reported with the length of its translated model.

Usage:  python analysis/absence_check.py
Writes analysis/results/absence_check.csv
"""

import glob
import os

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..")
RES = os.path.join(HERE, "results")
D = os.path.join(ROOT, "data", "absence_check")
NAMES = {"GCF_028476545.2": "Physella acuta", "GCF_000002075.1": "Aplysia californica",
         "GCF_034508935.2": "Berghia stephanieae", "GCF_046866305.1": "Octopus bimaculoides",
         "GCF_016097555.1": "Gigantopelta aegis", "control_LG17": "Littorina saxatilis LG17 (control)",
         "nautilus": "Nautilus pompilius GCA_018389105.1 (no annotation)"}


def parse(gff):
    rows, sta = [], None
    for line in open(gff):
        if line.startswith("##STA"):
            sta = line.split("\t")[1].strip()
        elif "\tmRNA\t" in line:
            c = line.rstrip("\n").split("\t")
            a = dict(x.split("=", 1) for x in c[8].split(";") if "=" in x)
            t = a["Target"].split()[0]
            rows.append(dict(contig=c[0], start=int(c[3]), end=int(c[4]), score=float(c[5]),
                             group=t.split("|")[0], target=t, identity=float(a.get("Identity", 0)),
                             frameshift=int(a.get("Frameshift", 0)), stop="StopCodon" in a,
                             aa=len(sta.replace("*", "")) if sta else 0))
    return pd.DataFrame(rows)


def loci(d, gap=5000):
    d = d.sort_values(["contig", "start"]).reset_index(drop=True)
    lid, cur, out = 0, None, []
    for r in d.itertuples():
        if cur is None or r.contig != cur[0] or r.start > cur[1] + gap:
            lid += 1
            cur = [r.contig, r.end]
        else:
            cur[1] = max(cur[1], r.end)
        out.append(lid)
    d["locus"] = out
    return d


def main():
    rows = []
    for gff in sorted(glob.glob(os.path.join(D, "*.gff"))):
        acc = os.path.basename(gff)[:-4]
        d = parse(gff)
        if d.empty:
            continue
        d = loci(d)
        for lid, g in d.groupby("locus"):
            best = g.sort_values("score", ascending=False).iloc[0]
            by = g.groupby("group").score.max()
            rows.append(dict(species=NAMES.get(acc, acc), locus=lid, contig=best.contig, start=best.start,
                             best_group=best.group, best_score=best.score,
                             LIN_score=by.get("LIN", 0), DISP12_score=by.get("DISP12", 0),
                             DISP3_score=by.get("DISP3", 0), PTCH_score=by.get("PTCH", 0),
                             model_aa=best.aa, identity=best.identity, frameshift=best.frameshift,
                             stop=best.stop))
    t = pd.DataFrame(rows)
    t.to_csv(os.path.join(RES, "absence_check.csv"), index=False)
    pd.set_option("display.width", 220)
    print(t.to_string(index=False))
    print("\nLoci whose best match is the mollusc-specific lineage (LIN):")
    print(t[t.best_group == "LIN"].groupby("species").size().to_string())


if __name__ == "__main__":
    main()
