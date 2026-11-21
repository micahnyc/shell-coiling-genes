"""
Place the QTL markers on the 2024 chromosome-level genome and list the genes
inside each QTL region.

The Koch et al. 2021 markers are named after contigs of the 2018 genome
(v1, Westram et al. 2018). Following De Jode et al. 2024, we take 200 bp of v1
sequence centred on each SNP and find where it lies on the 2024 genome
(v2, GCF_037325665.1, chromosomes named LG1-LG17). Instead of an aligner we
use exact 27-mer matching (no extra software needed): a fragment is placed
where most of its 27-mers hit, if that spot clearly beats any other.

Genes come from NCBI's RefSeq annotation of v2 (GCF_037325665.1-RS_2024_12).

Inputs
  data/littorina_koch2021/QTL.map.txt        marker -> v1 contig + position
  data/littorina_v1/<assembly>.fa            v1 genome (Zenodo 5003441)
  data/littorina_v2/<chrom>.fa, genomic.gff  v2 chromosomes and annotation
  analysis/results/lod_scores.csv.gz         from littorina_qtl.py

Outputs (analysis/results/)
  marker_v2_positions.csv   each marker's position on v2
  qtl_regions.csv           each QTL region in cM and in base pairs
  qtl_region_genes.csv      every gene in each region, with candidate flags
  cM_vs_Mb.png              genetic vs physical map (flat runs = low recombination,
                            typical of inversions)
"""

import glob
import os
import re
import zipfile

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..")
RES = os.path.join(HERE, "results")
K = 27
FLANK = 100

CHROM = {17: "NC_090261.1", 6: "NC_090250.1"}  # Crab-map LG -> v2 RefSeq chromosome
REGIONS = [  # (trait, LG): 1.5-LOD support interval around the LG's highest LOD
    ("W", 17), ("D_like", 17), ("W_height", 17), ("S", 6), ("shellLength", 6),
]

# Genes known or suspected to act in shell formation or in setting growth
# direction in molluscs, and general growth/patterning signals.
CANDIDATE_PATTERNS = [
    r"decapentaplegic|bone morphogenetic|\bBMP|TGF|transforming growth factor|nodal|lefty|activin",
    r"\bSMAD|mothers against", r"\bwnt|frizzled|dishevelled", r"hedgehog|patched|dispatched",
    r"notch|delta-like protein|jagged", r"homeobox|\bhox|engrailed|pitx|goosecoid|brachyury|dachshund",
    r"chitin|chitinase", r"tyrosinase", r"carbonic anhydrase", r"perlucin|nacrein|shematrin|pif|lustrin",
    r"collagen|matrilin|von willebrand|fibronectin", r"\bmucin",
    r"calmodulin|calcium-binding|calcium.*(channel|transporter|ATPase)|sodium/calcium",
    r"formin|inverted formin|formin-", r"insulin|growth factor", r"\bEGF\b|epidermal growth factor", r"peroxidase",
]
CAND_RE = re.compile("|".join(CANDIDATE_PATTERNS), re.I)

LUT = np.full(256, 4, np.uint8)
for i, c in enumerate(b"ACGT"):
    LUT[c] = i
    LUT[c + 32] = i  # lower case


def read_fasta_stream(lines, wanted=None):
    name, buf = None, []
    for line in lines:
        if line.startswith(">"):
            if name is not None and (wanted is None or name in wanted):
                yield name, "".join(buf)
            name, buf = line[1:].split()[0], []
        elif wanted is None or name in wanted:
            buf.append(line.strip())
    if name is not None and (wanted is None or name in wanted):
        yield name, "".join(buf)


def v1_fragments(map_df):
    """200 bp around each SNP from the v1 assembly."""
    wanted = set(map_df.contig)
    files = glob.glob(os.path.join(ROOT, "data", "littorina_v1", "*.zip"))
    z = zipfile.ZipFile(files[0])
    fa = [n for n in z.namelist() if re.search(r"\.(fa|fasta|fna)$", n)][0]
    with z.open(fa) as fh:
        lines = (l.decode("ascii") for l in fh)
        seqs = dict(read_fasta_stream(lines, wanted))
    out = {}
    for r in map_df.itertuples():
        s = seqs.get(r.contig)
        if s is None:
            continue
        p = int(r.pos) - 1  # 1-based SNP position
        out[r.marker] = s[max(0, p - FLANK): p + FLANK]
    return out


def kmer_codes(seq_bytes):
    """All K-mers of a sequence as uint64 (2 bits per base); invalid -> max value."""
    a = LUT[np.frombuffer(seq_bytes, np.uint8)]
    n = len(a) - K + 1
    if n <= 0:
        return np.zeros(0, np.uint64)
    bad = np.concatenate([[0], np.cumsum(a == 4)])
    invalid = (bad[K:] - bad[:-K]) > 0
    v = np.zeros(n, np.uint64)
    a64 = a.astype(np.uint64) & np.uint64(3)
    for j in range(K):
        v = (v << np.uint64(2)) | a64[j:j + n]
    v[invalid] = np.iinfo(np.uint64).max
    return v


def revcomp(s):
    return s.translate(str.maketrans("ACGTacgtN", "TGCAtgcaN"))[::-1]


def place(fragments, chrom_fa):
    # index query k-mers (both strands); drop k-mers shared by several fragments
    owner = {}
    for m, s in fragments.items():
        for strand in (s, revcomp(s)):
            for k in set(kmer_codes(strand.encode()).tolist()):
                if k == np.iinfo(np.uint64).max:
                    continue
                owner[k] = m if owner.get(k, m) == m else None
    owner = {k: m for k, m in owner.items() if m is not None}
    keys = np.array(sorted(owner), np.uint64)
    owners = np.array([owner[k] for k in keys.tolist()])

    seq = "".join(l.strip() for l in open(chrom_fa) if not l.startswith(">")).encode()
    hits_m, hits_p = [], []
    step = 10_000_000
    for start in range(0, len(seq), step):
        chunk = seq[start: start + step + K - 1]
        v = kmer_codes(chunk)
        idx = np.searchsorted(keys, v)
        idx[idx >= len(keys)] = 0
        ok = keys[idx] == v
        hits_m.append(owners[idx[ok]])
        hits_p.append(np.nonzero(ok)[0] + start)
    hm, hp = np.concatenate(hits_m), np.concatenate(hits_p)

    rows = []
    df = pd.DataFrame(dict(marker=hm, pos=hp, bin=hp // 1000))
    for m, g in df.groupby("marker"):
        counts = g.bin.value_counts()
        best = counts.iloc[0]
        second = counts.iloc[1] if len(counts) > 1 else 0
        if best >= 20 and best >= 3 * second:
            b = counts.index[0]
            near = g[(g.bin >= b - 1) & (g.bin <= b + 1)].pos
            rows.append(dict(marker=m, bp=int(near.median()), kmer_hits=int(best)))
    return pd.DataFrame(rows)


def genes_in(gff, chrom, start, end):
    rows = []
    with open(gff, encoding="utf-8") as f:
        for line in f:
            if not line.startswith(chrom + "\t"):
                continue
            c = line.rstrip("\n").split("\t")
            if c[2] != "gene":
                continue
            s, e = int(c[3]), int(c[4])
            if e < start or s > end:
                continue
            a = dict(kv.split("=", 1) for kv in c[8].split(";") if "=" in kv)
            rows.append(dict(gene=a.get("Name"), gene_id=a.get("Dbxref", "").replace("GeneID:", ""),
                             start=s, end=e, strand=c[6], biotype=a.get("gene_biotype"),
                             description=a.get("description", "").replace("%2C", ",")))
    return pd.DataFrame(rows)


def main():
    qmap = pd.read_csv(os.path.join(ROOT, "data", "littorina_koch2021", "QTL.map.txt"), sep="\t")
    qmap["marker"] = qmap.index
    lods = pd.read_csv(os.path.join(RES, "lod_scores.csv.gz"))
    lods = lods.merge(qmap[["marker", "contig", "pos"]], on="marker")

    placed = []
    for lg, chrom in CHROM.items():
        sub = lods[lods.LG == lg]
        frags = v1_fragments(sub)
        print(f"LG{lg}: {len(sub)} markers, {len(frags)} v1 fragments extracted")
        pl = place(frags, os.path.join(ROOT, "data", "littorina_v2", chrom + ".fa"))
        pl["LG"], pl["chrom"] = lg, chrom
        print(f"LG{lg}: {len(pl)} markers placed uniquely on {chrom}")
        placed.append(sub.merge(pl, on=["marker", "LG"]))
    placed = pd.concat(placed)
    placed[["marker", "LG", "cM", "contig", "pos", "chrom", "bp", "kmer_hits"]].to_csv(
        os.path.join(RES, "marker_v2_positions.csv"), index=False)

    gff = os.path.join(ROOT, "data", "littorina_v2", "genomic.gff")
    regions, gene_tables = [], []
    for trait, lg in REGIONS:
        L = lods[lods.LG == lg]
        peak = L[trait].max()
        lo, hi = L.cM[L[trait] >= peak - 1.5].agg(["min", "max"])
        p = placed[(placed.LG == lg) & placed.cM.between(lo, hi)]
        if len(p) < 2:
            print(f"{trait} LG{lg}: too few placed markers in {lo}-{hi} cM")
            continue
        start, end = int(p.bp.quantile(0.05)), int(p.bp.quantile(0.95))
        g = genes_in(gff, CHROM[lg], start, end)
        g["candidate"] = g.description.apply(lambda d: bool(CAND_RE.search(d)))
        g.insert(0, "region", f"{trait}_LG{lg}")
        gene_tables.append(g)
        pc = (g.biotype == "protein_coding").sum()
        regions.append(dict(trait=trait, LG=lg, chrom=CHROM[lg], peak_LOD=round(peak, 2),
                            cM_from=lo, cM_to=hi, markers_placed=len(p),
                            Mb_from=round(start / 1e6, 2), Mb_to=round(end / 1e6, 2),
                            protein_coding_genes=int(pc), candidate_genes=int(g.candidate.sum())))
    reg = pd.DataFrame(regions)
    reg.to_csv(os.path.join(RES, "qtl_regions.csv"), index=False)
    genes = pd.concat(gene_tables)
    genes.to_csv(os.path.join(RES, "qtl_region_genes.csv"), index=False)
    print("\n", reg.to_string(index=False))
    print("\nCandidate genes:")
    c = genes[genes.candidate & (genes.biotype == "protein_coding")]
    print(c.drop_duplicates("gene")[["region", "gene", "start", "description"]].to_string(index=False))

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    for ax, lg in zip(axes, CHROM):
        p = placed[placed.LG == lg]
        ax.plot(p.bp / 1e6, p.cM, ".", ms=3, color="#8a3b12")
        for r in regions:
            if r["LG"] == lg:
                ax.axvspan(r["Mb_from"], r["Mb_to"], alpha=0.15, label=f"{r['trait']} QTL")
        ax.set_xlabel(f"position on 2024 genome {CHROM[lg]} (Mb)")
        ax.set_ylabel("genetic map position (cM)")
        ax.set_title(f"LG{lg}")
        ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(RES, "cM_vs_Mb.png"), dpi=120)


if __name__ == "__main__":
    main()
