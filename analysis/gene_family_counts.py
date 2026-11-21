"""
Step 2: do the candidate gene families from the Littorina shell-shape QTL
regions vary across molluscs with different shells?

For each mollusc genome with an NCBI RefSeq annotation, count protein-coding
genes in each family by their RefSeq gene names (NCBI Datasets
annotation_report API). Counts by name are approximate: RefSeq names come from
the same automated pipeline for every species, which keeps them comparable,
but a name like "inverted formin-2-like" is a similarity label, not proven
orthology.

Usage:  python analysis/gene_family_counts.py
Writes analysis/results/gene_family_counts.csv and gene_family_counts.png
"""

import json
import os
import re
import time
import urllib.parse
import urllib.request

import pandas as pd
import matplotlib.pyplot as plt

RES = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
API = "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/{acc}/annotation_report"

# (accession, species, group, shell)
SPECIES = [
    ("GCF_037325665.1", "Littorina saxatilis", "snail", "coiled, right-handed"),
    ("GCF_003073045.1", "Pomacea canaliculata", "snail", "coiled, right-handed"),
    ("GCF_041734735.1", "Babylonia areolata", "snail", "coiled, right-handed"),
    ("GCF_016097555.1", "Gigantopelta aegis", "snail", "coiled, right-handed"),
    ("GCF_947242115.1", "Biomphalaria glabrata", "snail", "coiled flat, left-handed"),
    ("GCF_028476545.2", "Physella acuta", "snail", "coiled, left-handed"),
    ("GCF_023055435.1", "Haliotis rufescens", "snail", "ear-shaped (abalone)"),
    ("GCF_000327385.1", "Lottia gigantea", "snail", "cap (limpet)"),
    ("GCF_932274485.2", "Patella vulgata", "snail", "cap (limpet)"),
    ("GCF_000002075.1", "Aplysia californica", "snail", "reduced, internal"),
    ("GCF_034508935.2", "Berghia stephanieae", "snail", "none (sea slug)"),
    ("GCF_021730395.1", "Mercenaria mercenaria", "clam", "two valves"),
    ("GCF_026914265.1", "Mya arenaria", "clam", "two valves"),
    ("GCF_026571515.1", "Ruditapes philippinarum", "clam", "two valves"),
    ("GCF_053477285.1", "Crassostrea virginica", "clam", "two valves (oyster)"),
    ("GCF_902652985.1", "Pecten maximus", "clam", "two valves (scallop)"),
    ("GCF_032854445.1", "Liolophura sinensis", "chiton", "eight plates"),
    ("GCF_053919585.1", "Euprymna scolopes", "cephalopod", "reduced, internal"),
    ("GCF_046866305.1", "Octopus bimaculoides", "cephalopod", "none"),
    ("GCF_006345805.1", "Octopus sinensis", "cephalopod", "none"),
]

# family: (search text for the API, regex the gene name must match)
FAMILIES = {
    "diaphanous formin": ("diaphanous", r"diaphanous"),
    "inverted formin (INF2)": ("inverted formin", r"inverted formin"),
    "formin-2 (FMN2)": ("formin-2", r"(?<!inverted )formin-2"),
    "dispatched": ("dispatched", r"dispatched"),
    "perlucin": ("perlucin", r"perlucin"),
}


def get(url):
    for attempt in range(5):
        try:
            with urllib.request.urlopen(url, timeout=60) as r:
                return json.load(r)
        except Exception:
            time.sleep(2 + 3 * attempt)
    raise RuntimeError(url)


def count(acc, search, pattern):
    names, token = set(), None
    while True:
        q = {"search_text": search, "page_size": 1000}
        if token:
            q["page_token"] = token
        d = get(API.format(acc=acc) + "?" + urllib.parse.urlencode(q))
        for rep in d.get("reports", []):
            a = rep["annotation"]
            if a.get("gene_type") == "protein-coding" and re.search(pattern, a.get("name", ""), re.I):
                names.add(a["gene_id"])
        token = d.get("next_page_token")
        if not token:
            return len(names)
        time.sleep(0.35)


def main():
    rows = []
    for acc, sp, grp, shell in SPECIES:
        row = dict(species=sp, group=grp, shell=shell, accession=acc)
        for fam, (search, pat) in FAMILIES.items():
            row[fam] = count(acc, search, pat)
            time.sleep(0.35)  # stay under NCBI's 3 requests/second
        rows.append(row)
        print(sp, {f: row[f] for f in FAMILIES})
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(RES, "gene_family_counts.csv"), index=False)

    fam = list(FAMILIES)
    fig, ax = plt.subplots(figsize=(10, 8))
    im = ax.imshow(df[fam].values, cmap="YlOrBr", aspect="auto")
    for i in range(len(df)):
        for j, f in enumerate(fam):
            ax.text(j, i, df[f].iloc[i], ha="center", va="center", fontsize=9)
    ax.set_xticks(range(len(fam)))
    ax.set_xticklabels(fam, rotation=20, ha="right")
    ax.set_yticks(range(len(df)))
    ax.set_yticklabels([f"{r.species} ({r.group}; {r.shell})" for r in df.itertuples()], fontsize=8)
    ax.set_title("Genes per family (RefSeq annotation names)")
    fig.colorbar(im, ax=ax, shrink=0.6)
    fig.tight_layout()
    fig.savefig(os.path.join(RES, "gene_family_counts.png"), dpi=120)
    print("saved gene_family_counts.png")


if __name__ == "__main__":
    main()
