"""
Step 2b: gene trees for the two candidate families (formins incl. diaphanous,
and dispatched), to replace name-based counts with tree-based orthology.

1. fetch  For every species in gene_family_counts.SPECIES plus human and fruit fly
          (well-studied references), collect protein-coding genes whose RefSeq
          names belong to the family (NCBI Datasets annotation_report), keep the
          longest protein isoform per gene, download sequences (NCBI efetch).
          Lymnaea Lsdia1/Lsdia2 (the pond-snail coiling genes) are added.
2. tree   (in WSL) MAFFT alignment -> trimAl -automated1 -> IQ-TREE
          (ModelFinder over LG/WAG/JTT, 1000 ultrafast bootstraps).
3. read   Root on the outgroup, assign each mollusc gene to the nearest
          reference subfamily (human/fly gene names) by tree distance, count per
          species, and draw the trees.

Usage:
  python analysis/gene_trees.py fetch
  wsl -d Ubuntu -- bash /mnt/c/micah/shells/analysis/run_trees.sh
  python analysis/gene_trees.py read
"""

import json
import os
import re
import sys
import time
import urllib.parse
import urllib.request

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gene_family_counts import SPECIES  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..")
RES = os.path.join(HERE, "results")
OUT = os.path.join(ROOT, "data", "gene_trees")
API = "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/{acc}/annotation_report"
EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
MIN_LEN = 250  # drop fragments

REFS = [("GCF_000001405.40", "Homo sapiens", "human", "reference"),
        ("GCF_000001215.4", "Drosophila melanogaster", "fly", "reference")]

FAMILIES = {
    "formin": dict(
        searches=["diaphanous", "formin", "dishevelled", "activator of morphogenesis", "FH2",
                  "delphilin", "cappuccino", "Fhos", "inverted formin"],
        name_re=r"diaphanous|formin(?!g)|activator of morphogenesis|FH2 domain|delphilin|"
                r"cappuccino|\bFhos\b|grid2-interacting|formin homology",
        exclude_re=r"forming|binding protein|interacting protein(?! .*grid2)|fhod-?interacting"),
    "dispatched": dict(
        searches=["dispatched"],
        name_re=r"dispatched",
        exclude_re=r"$^"),
}
# outgroup for the dispatched tree: human patched 1 (same sterol-sensing protein family)
EXTRA = {"dispatched": [("GCF_000001405.40", "Homo sapiens", "patched 1", r"^patched 1$")]}
LYMNAEA = {"AOV18871.1": "Lymnaea_stagnalis|Lsdia1", "AOV18874.1": "Lymnaea_stagnalis|Lsdia2"}


def get(url):
    for attempt in range(6):
        try:
            with urllib.request.urlopen(url, timeout=90) as r:
                return r.read()
        except Exception:
            time.sleep(2 + 3 * attempt)
    raise RuntimeError(url)


def reports(acc, search):
    token = None
    while True:
        q = {"search_text": search, "page_size": 1000}
        if token:
            q["page_token"] = token
        d = json.loads(get(API.format(acc=acc) + "?" + urllib.parse.urlencode(q)))
        yield from d.get("reports", [])
        token = d.get("next_page_token")
        time.sleep(0.35)
        if not token:
            return


def longest_protein(a):
    best = None
    for t in a.get("transcripts", []):
        p = t.get("protein")
        if p and (best is None or p.get("length", 0) > best.get("length", 0)):
            best = p
    return best


def efetch(accs):
    seqs = {}
    for i in range(0, len(accs), 150):
        data = urllib.parse.urlencode({"db": "protein", "id": ",".join(accs[i:i + 150]),
                                       "rettype": "fasta", "retmode": "text"}).encode()
        for attempt in range(6):
            try:
                with urllib.request.urlopen(EFETCH, data=data, timeout=120) as r:
                    txt = r.read().decode()
                break
            except Exception:
                time.sleep(3 + 3 * attempt)
        name = None
        for line in txt.splitlines():
            if line.startswith(">"):
                name = line[1:].split()[0]
                seqs[name] = []
            elif name:
                seqs[name].append(line.strip())
        time.sleep(0.4)
    return {k: "".join(v) for k, v in seqs.items()}


def fetch():
    os.makedirs(OUT, exist_ok=True)
    for fam, spec in FAMILIES.items():
        rows = []
        for acc, sp, grp, shell in SPECIES + REFS:
            seen = set()
            for s in spec["searches"]:
                for rep in reports(acc, s):
                    a = rep["annotation"]
                    name = a.get("name", "")
                    if a.get("gene_type") != "protein-coding" or a["gene_id"] in seen:
                        continue
                    if not re.search(spec["name_re"], name, re.I) or re.search(spec["exclude_re"], name, re.I):
                        continue
                    p = longest_protein(a)
                    if not p or p.get("length", 0) < MIN_LEN:
                        continue
                    seen.add(a["gene_id"])
                    chrom = (a.get("chromosomes") or [""])[0]
                    rows.append(dict(family=fam, species=sp, group=grp, shell=shell, gene_id=a["gene_id"],
                                     symbol=a.get("symbol"), name=name, chrom=chrom,
                                     protein=p["accession_version"], length=p["length"]))
            print(f"{fam:10s} {sp:28s} {len(seen)}")
        for acc, sp, search, pat in EXTRA.get(fam, []):
            for rep in reports(acc, search):
                a = rep["annotation"]
                if re.search(pat, a.get("name", ""), re.I):
                    p = longest_protein(a)
                    rows.append(dict(family=fam, species=sp, group="human", shell="outgroup", gene_id=a["gene_id"],
                                     symbol=a.get("symbol"), name=a["name"], chrom="",
                                     protein=p["accession_version"], length=p["length"]))
                    break
        df = pd.DataFrame(rows)
        seqs = efetch(df.protein.tolist() + (list(LYMNAEA) if fam == "formin" else []))
        df["label"] = [f"{r.species.replace(' ', '_')}|{r.symbol}|{r.protein}" for r in df.itertuples()]
        with open(os.path.join(OUT, f"{fam}.faa"), "w") as f:
            for r in df.itertuples():
                f.write(f">{r.label}\n{seqs[r.protein]}\n")
            if fam == "formin":
                for acc, lab in LYMNAEA.items():
                    f.write(f">{lab}|{acc}\n{seqs[acc]}\n")
        df.to_csv(os.path.join(OUT, f"{fam}_genes.csv"), index=False)
        print(f"{fam}: {len(df)} genes written ({len(seqs)} sequences fetched)")


# reference gene symbol -> subfamily
SUBFAMILY = [
    (r"^(DIAPH\d|dia|Lsdia\d)$", "diaphanous"), (r"^DAAM\d?$", "DAAM"), (r"^(FMNL\d|Frl)$", "FMNL"),
    (r"^(FHOD\d|Fhos)$", "FHOD"), (r"^(INF2|form3)$", "INF2"), (r"^(FMN\d|capu)$", "FMN"),
    (r"^FHDC1$", "FHDC"), (r"^DISP1$", "DISP1"), (r"^DISP2$", "DISP2"), (r"^DISP3$", "DISP3"),
    (r"^disp$", "fly disp"), (r"^PTCH1$", "patched (outgroup)"),
]
REF_SPECIES = ("Homo_sapiens", "Drosophila_melanogaster", "Lymnaea_stagnalis")
CAENO = ("Littorina_saxatilis", "Pomacea_canaliculata", "Babylonia_areolata", "Gigantopelta_aegis")


def subfamily(symbol):
    for pat, lab in SUBFAMILY:
        if re.match(pat, symbol):
            return lab
    return None


def read():
    from Bio import Phylo
    import matplotlib.pyplot as plt

    summary = []
    for fam in FAMILIES:
        tree = Phylo.read(os.path.join(OUT, f"{fam}.contree"), "newick")
        leaves = tree.get_terminals()
        if fam == "dispatched":
            tree.root_with_outgroup([l for l in leaves if "|PTCH1|" in l.name][0])
        else:
            tree.root_at_midpoint()
        refs = {l.name: subfamily(l.name.split("|")[1]) for l in leaves
                if l.name.startswith(REF_SPECIES) and subfamily(l.name.split("|")[1])}

        # monophyly of each reference subfamily (with bootstrap of its smallest clade)
        print(f"\n=== {fam} tree: {len(leaves)} sequences ===")
        for sf in sorted(set(refs.values())):
            members = [n for n, s in refs.items() if s == sf]
            if len(members) < 2:
                continue
            ca = tree.common_ancestor(members)
            others = {refs[l.name] for l in ca.get_terminals() if l.name in refs} - {sf}
            print(f"  {sf:12s} refs={len(members)} clade size={len(ca.get_terminals()):3d} "
                  f"bootstrap={ca.confidence} other ref subfamilies inside={sorted(others) or 'none'}")

        rows = []
        ref_leaves = [l for l in leaves if l.name in refs]
        for l in leaves:
            if l.name.startswith(REF_SPECIES) or "|PTCH1|" in l.name:
                continue
            d = {r.name: tree.distance(l, r) for r in ref_leaves}
            nearest = min(d, key=d.get)
            sp, sym, prot = l.name.split("|")
            rows.append(dict(family=fam, species=sp.replace("_", " "), symbol=sym, protein=prot,
                             nearest_ref=nearest.split("|")[1], subfamily=refs[nearest],
                             distance=round(d[nearest], 3)))
        a = pd.DataFrame(rows)
        genes = pd.read_csv(os.path.join(OUT, f"{fam}_genes.csv"))
        a = a.merge(genes[["protein", "name", "chrom", "group", "shell"]], on="protein", how="left")
        a.to_csv(os.path.join(RES, f"gene_tree_{fam}_assignments.csv"), index=False)
        counts = a.pivot_table(index=["species", "group", "shell"], columns="subfamily",
                               values="protein", aggfunc="count", fill_value=0)
        print(counts.to_string())
        counts.insert(0, "family", fam)
        summary.append(counts.reset_index())

        # figure
        n = len(leaves)
        fig, ax = plt.subplots(figsize=(11, max(8, n * 0.16)))
        colour = lambda name: ("#c0392b" if name.startswith("Littorina") else
                               "#1f6f8b" if name.startswith(REF_SPECIES) else "black")
        Phylo.draw(tree, axes=ax, do_show=False, label_colors=colour,
                   label_func=lambda c: c.name.replace("_", " ") if c.is_terminal() else "",
                   branch_labels=lambda c: (f"{c.confidence:.0f}" if c.confidence and not c.is_terminal()
                                            and c.confidence >= 95 and len(c.get_terminals()) > 3 else ""))
        for t in ax.texts:
            t.set_fontsize(6)
        ax.set_title(f"{fam}: maximum-likelihood tree (IQ-TREE, numbers = ultrafast bootstrap >= 95). "
                     f"Red = Littorina, blue = human/fly/Lymnaea references", fontsize=9)
        ax.set_ylabel("")
        fig.tight_layout()
        fig.savefig(os.path.join(RES, f"gene_tree_{fam}.png"), dpi=110)
        plt.close(fig)
    pd.concat(summary).to_csv(os.path.join(RES, "gene_tree_counts.csv"), index=False)


if __name__ == "__main__":
    {"fetch": fetch, "read": read}[sys.argv[1]]()
