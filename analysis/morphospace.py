"""
Step 3: where do the species with sequenced genomes sit in Raup's shell morphospace,
and do the candidate gene counts line up with shell geometry?

Shell parameters: Noshita, Asami & Ubukata 2012, Paleobiology 38:322, Dryad
doi:10.5061/dryad.6mk023c0, Supplementary Table 1 (PDF): Raup's W, T and D measured on
museum specimens of ~350 living gastropod species. (Noshita's D is signed: negative
means the tube overlaps the axis, as in shells with a closed umbilicus.)
Nautilus is added from Tajika & Klug 2020 (PeerJ 8:e8849): W ~2.9, T = 0, D ~0.

Genome species are matched to Noshita species directly where possible, otherwise to
congeners (or a same-family proxy), and flagged as such.
Gene counts come from the gene trees (gene_tree_counts.csv).

With so few species with both kinds of data, this is a description, not a test.

Usage:  python analysis/morphospace.py
Writes analysis/results/noshita2012_raup_parameters.csv, morphospace_genome_species.csv,
       morphospace.png
"""

import os
import re

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pypdf

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..")
RES = os.path.join(HERE, "results")
PDF = os.path.join(ROOT, "data", "noshita2012", "Table_without_sp+Supplementary+Table.+1.pdf")

ROW = re.compile(r"^(?P<species>[A-Z][a-z]+(?: [a-z][\w.-]*)+?(?: (?:form|var\.) [\w-]+)?)\s+"
                 r"(?P<repo>SPMN|SUM|SHIN|BPBM|UMUT)\s+(?P<habitat>[MFL])\s+"
                 r"(?P<W>-?[\d.]+)\s+(?P<T>-?[\d.]+)\s+(?P<D>-?[\d.]+)")

# genome species -> (Noshita species to use, match type)
MATCH = {
    "Pomacea canaliculata": (["Pomacea canaliculata"], "same species"),
    "Physella acuta": (["Physa acuta"], "same species (Physa = Physella)"),
    "Littorina saxatilis": (["Littorina littorea", "Littorina undulata", "Littorina zebra"], "congeners"),
    "Haliotis rufescens": (["Haliotis gigantea", "Haliotis diversicolor aquatilis"], "congeners"),
    "Biomphalaria glabrata": (["Planorbarius corneus"], "same family (Planorbidae)"),
}
REFERENCE = {"Lymnaea stagnalis": "Lsdia1 coiling gene; no RefSeq genome"}


def parse():
    text = "\n".join(p.extract_text() for p in pypdf.PdfReader(PDF).pages)
    rows, family = [], None
    for line in text.splitlines():
        line = line.strip()
        if re.fullmatch(r"[A-Z][a-z]+idae", line):
            family = line
        m = ROW.match(line)
        if m:
            r = m.groupdict()
            rows.append(dict(species=r["species"], family=family, repository=r["repo"],
                             habitat=r["habitat"], W=float(r["W"]), T=float(r["T"]), D=float(r["D"])))
    return pd.DataFrame(rows)


def main():
    t = parse()
    t.to_csv(os.path.join(RES, "noshita2012_raup_parameters.csv"), index=False)
    sp = t.groupby("species")[["W", "T", "D"]].mean()
    print(f"parsed {len(t)} specimens, {len(sp)} species, {t.family.nunique()} families")

    counts = pd.read_csv(os.path.join(RES, "gene_tree_counts.csv"))
    disp = counts[counts.family == "dispatched"].set_index("species")
    form = counts[counts.family == "formin"].set_index("species")
    rows = []
    for g, (names, how) in MATCH.items():
        v = sp.loc[names].mean()
        rows.append(dict(genome_species=g, noshita_species="; ".join(names), match=how,
                         W=v["W"], T=v["T"], D=v["D"],
                         mollusc_specific_dispatched=disp.loc[g, "mollusc-specific dispatched-like lineage"],
                         diaphanous=form.loc[g, "diaphanous"], shell=disp.loc[g, "shell"]))
    rows.append(dict(genome_species="Nautilus pompilius", noshita_species="Tajika & Klug 2020", match="literature",
                     W=2.9, T=0.0, D=0.0, mollusc_specific_dispatched=np.nan, diaphanous=np.nan,
                     shell="planispiral, chambered (no annotated genome)"))
    g = pd.DataFrame(rows)
    ly = sp.loc["Lymnaea stagnalis"]
    g.to_csv(os.path.join(RES, "morphospace_genome_species.csv"), index=False, float_format="%.3f")
    pd.set_option("display.width", 220)
    print(g.round(3).to_string(index=False))
    print(f"\nLymnaea stagnalis (reference, Lsdia1): W={ly['W']:.3f} T={ly['T']:.3f} D={ly['D']:.3f}")

    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))
    for ax, (x, y, xl, yl) in zip(axes, [("W", "T", "W, whorl expansion per turn (log scale)", "T, translation (spire height)"),
                                         ("W", "D", "W, whorl expansion per turn (log scale)", "D (Noshita, signed)")]):
        ax.scatter(sp[x], sp[y], s=10, c="#bbbbbb", lw=0, label=f"{len(sp)} gastropod species (Noshita 2012)")
        for r in g.itertuples():
            ax.scatter(getattr(r, x), getattr(r, y), s=60, c="#b8741a", edgecolor="black", lw=0.6, zorder=3)
            ax.annotate(r.genome_species + ("" if r.match.startswith("same species") or r.match == "literature" else "*"),
                        (getattr(r, x), getattr(r, y)), fontsize=8, xytext=(5, 4), textcoords="offset points")
        ax.scatter(ly[x], ly[y], s=60, c="#3e6f8e", edgecolor="black", lw=0.6, zorder=3)
        ax.annotate("Lymnaea stagnalis", (ly[x], ly[y]), fontsize=8, xytext=(-8, 0), textcoords="offset points",
                    ha="right", va="center", color="#3e6f8e")
        ax.set_xscale("log")
        ax.set_xlabel(xl)
        ax.set_ylabel(yl)
    axes[0].legend(fontsize=8, loc="upper right")
    fig.suptitle("Raup morphospace with the genome species (orange; * = congener or family proxy) "
                 "and Lymnaea (blue)", fontsize=10)
    fig.tight_layout()
    fig.savefig(os.path.join(RES, "morphospace.png"), dpi=120)


if __name__ == "__main__":
    main()
