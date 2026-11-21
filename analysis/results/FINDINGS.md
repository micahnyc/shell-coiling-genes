# Littorina saxatilis: genome regions linked to Raup-type shell parameters

Data: Koch et al. 2021 (Evolution Letters 5:196-213), Dryad doi:10.5061/dryad.zgmsbccb4.
F2 snails from Crab x Wave ecotype crosses: 350 with shape data, 13 families,
19,849 SNP markers on 17 linkage groups (LG, Crab-map numbering).
Shape measured with the Larsson et al. 2020 growth model, translated to Raup terms.
Script: `analysis/littorina_qtl.py`.

## Reproduction check
Using the paper's covariates (batch, sex, shell length), our marker regression gives
the same peaks as the paper's Haley-Knott scans in R/qtl2:

| Trait | Ours (LOD, LG) | Paper (LOD, LG) |
|---|---|---|
| heightGrowth | 4.38, LG17 | 4.16, LG17 |
| aperturePosition | 6.36, LG17 | 6.16, LG17 |
| apertureShape | 4.07, LG6 | 4.36, LG6 |

## Raup parameters (stricter scan: family added as a covariate, 1000 permutations)

| Raup parameter | Strongest region | LOD (5% threshold) | Status |
|---|---|---|---|
| **W** (expansion per turn) | LG17, 64-67 cM | 4.43 (3.88) | significant, 5.7% of variance |
| **D** (distance from axis) | LG17, 48-64 cM | 3.95 (3.93) | at threshold (raw aperture position: 5.56, clearly significant) |
| **S** (opening elongation) | LG6 | 3.82 (3.99) | just below threshold |
| Height expansion (T-related) | LG17 start / LG2 | 2.86 (3.94) | not significant with family included |
| Opening size | LG12 | 3.44 (3.89) | suggestive (paper also saw LG12) |
| Shell length (size) | LG6 | 3.23 (3.92) | suggestive |

## What it means
- **Coiling geometry (W and D) maps to one region, LG17.** The paper places the
  W and aperture-position QTLs inside a putative chromosomal inversion there.
- **Opening shape (S) maps to a different region, LG6,** together with size and weight,
  also an inversion region.
- So Raup's parameters are not controlled as one block: the coil and the opening
  have partly separate genetic bases, each tied to an inversion. Inversions keep
  many genes inherited together, so the next question is which genes inside them.
- Effects are modest (about 5-7% of variance each); shape is polygenic.

## Step 1: genes inside the QTL regions
Script: `analysis/map_qtl_genes.py`. Markers were placed on the 2024 chromosome-level
genome (GCF_037325665.1) by matching 200 bp of 2018-assembly sequence around each SNP
(358 of 1,370 LG17 markers and 266 of 1,174 LG6 markers placed uniquely). Genes are from
NCBI's RefSeq annotation. Outputs: `qtl_regions.csv`, `qtl_region_genes.csv`, `cM_vs_Mb.png`.

| Region | Location on 2024 genome | Protein-coding genes |
|---|---|---|
| W (coil expansion), LG17 | 63.9-70.3 Mb | 174 |
| D (distance from axis), LG17 | 54.0-67.3 Mb | 347 |
| S (opening shape), LG6 | 1.4-25.7 Mb | 539 |
| Shell length, LG6 | 1.4-17.7 Mb | 384 |

- **The LG6 region is a ~17 Mb block with no recombination** (genetic map flat at ~53 cM
  from 0 to 17 Mb): the signature of the LG6 inversion. S and size QTLs sit in it.
- **Shell/signalling genes in general are not enriched** in any region (2.0-2.9% vs 2.9%
  genome-wide), so a long list of "shell-like" names in a region means little by itself.
- **Two gene families are concentrated** (found after looking, so p-values are optimistic,
  and copies sit side by side, so they are not independent):
  - **Formins in the LG6 shape/size region:** formin-2-like (1.8 Mb), four inverted
    formin-2-like genes (3.3, 9.9, 10.0, 10.3 Mb) and **protein diaphanous homolog 2-like
    (6.4 Mb)**. 5 of the genome's 12 formin-named genes, p ~ 2e-6. A diaphanous-type
    formin (Lsdia1) sets coiling direction in the pond snail Lymnaea (Davison et al. 2016).
  - **Dispatched in the LG17 W region:** 4 of the genome's 8 dispatched genes (p ~ 1e-7);
    7 of the 8 are on LG17 overall. Dispatched releases Hedgehog, a growth signal.
  - Perlucin, a nacre shell-matrix protein, also lies in the LG6 region (5.9 Mb).

### The LG6 diaphanous gene is Littorina's version of the pond-snail coiling gene
Global protein alignment (BLOSUM62), sequences in `data/formin_proteins/diaphanous.faa`:

| | Lymnaea Lsdia1 (AOV18871) | Lymnaea Lsdia2 (AOV18874) |
|---|---|---|
| Littorina dia-2-like, LG6 QTL (XP_070197016.1) | **65.7%** | **69.2%** |
| Littorina dia-1-like, LG3 (XP_070189802.1) | 34.8% | 33.7% |

The gene tree (below) confirms this: diaphanous is single-copy in most molluscs, and the
LG6 gene is Littorina's only member of the diaphanous clade. The LG3 "diaphanous-1-like"
gene is a misnamed 436-aa fragment that falls outside the clade.
Lsdia1 and Lsdia2 are 87.9% identical to each other (a duplication in the Lymnaea
lineage), and Littorina has one gene matching the pair, inside the LG6 inversion where
the opening-shape (S) and size QTLs map. In Lymnaea, losing Lsdia1 function flips the
coil from right- to left-handed (Davison et al. 2016). Whether this gene affects shell
shape in Littorina is untested: the LG6 region holds 539 genes inherited as one block.

### Which QTL genes are switched on in the mantle (shell-building tissue)
Script: `analysis/mantle_expression.py`. Gene-level TPMs from NCBI's RefSeq annotation
for a 10-sample tissue panel (PRJNA550990: mantle x2, foot x2, hepatopancreas x2,
reproductive x2, penis, head). Each gene is ranked by mantle bias, log2((mantle+1)/(others+1)),
against all 23,211 expressed genes. One or two samples per tissue: a ranking, not a test.
Outputs: `mantle_expression_genome.csv.gz`, `mantle_expression_qtl_genes.csv`, `mantle_expression.png`.

- **The data behave:** known shell-forming genes are at the top of the genome-wide ranking
  (carbonic anhydrases, tyrosinase, chitin synthase, nacrein, perlucin: 98th-99th percentile).
- **The QTL regions as a whole are not enriched** for mantle-biased genes (3-4% of each
  region's genes in the genome's top 5%).
- **Diaphanous-2-like (LG6)** is expressed in the mantle (21 TPM), but about equally in foot
  and head: a general cytoskeleton gene, not mantle-specific. This does not argue against it:
  in Lymnaea, Lsdia1 sets coiling direction from the egg (maternal expression, 8-cell
  stage), which adult tissues cannot show.
- **The LG17 dispatched cluster is nearly silent in the mantle** (0-1.4 TPM; what
  little expression there is is in hepatopancreas and reproductive organs). In adults at least,
  this weakens dispatched as a shell-growth candidate. The broadly expressed dispatched
  gene (dispatched-3-like, 20 TPM in mantle) is on another chromosome.

### The LG6 shape region in detail
The S (opening shape) LOD peaks at 14.8 Mb, inside the 0-17.7 Mb no-recombination block
(the inversion), but stays high to ~22 Mb (LOD 3.65 at 20-22 Mb, within noise of the
peak). Shell length is confined to the inversion (LOD drops from 3.2 to 2.1 past 17.7 Mb).
14 genes in the region are among the genome's 5% most mantle-biased:

| Mb | In inversion? | Gene (RefSeq name) | Mantle TPM | Other tissues |
|---|---|---|---|---|
| 5.94 | yes | perlucin-like protein (nacre shell-matrix protein family) | 16.8 | 6.7 |
| 9.42 | yes | UPAR/Ly6 domain-containing protein crok-like | **87.9** | 2.0 |
| 9.48 | yes | protocadherin-9-like | 9.9 | 2.6 |
| 10.20 | yes | S-acyl fatty acid synthase thioesterase-like | 4.0 | 0.1 |
| 13.16 | yes | beta-1,3-galactosyltransferase 1-like | 5.0 | 1.5 |
| 14.81 | yes | cytochrome P450 2U1-like (next to S peak) | 6.7 | 2.4 |
| 14.84 | yes | CD109 antigen-like (next to S peak) | 8.6 | 0.2 |
| 17.59 | yes | rho GTPase-activating protein 6-like | 11.4 | 1.5 |
| 18.87 | no | murinoglobulin-2-like | 8.1 | 0.4 |
| 19.25 | no | trichohyalin-like | 2.9 | 0.1 |
| 21.36 | no | CD109 antigen-like | 5.0 | 0.2 |
| 21.49 | no | uncharacterized LOC138968731 | 10.7 | 0.9 |
| 22.81 | no | solute carrier family 4 member 11-like (bicarbonate-transporter family) | 11.1 | 1.8 |
| 23.68 | no | arachidonate 12-lipoxygenase-like | 9.0 | 0.9 |

Reading this list:
- Mantle-specific is a better filter than gene names, but the mantle does many things
  besides building shell (mucus, immunity, sensing), so CD109/murinoglobulin
  (alpha-macroglobulin family) and lipoxygenase may be immune or mucus genes.
- An F2 cross cannot separate genes within the inversion: the 14.8 Mb peak is not
  evidence for the genes beside it over others in the block.
- Testing them needs either recombinants (wild snails where the inversion is
  broken up, or arrangement-matched comparisons) or expression at the growing shell edge
  in Crab vs Wave snails.

## Step 2: the same families across molluscs
Script: `analysis/gene_family_counts.py`; counts of RefSeq gene names in 20 molluscs
(`gene_family_counts.csv`, `gene_family_counts.png`). Nautilus has no NCBI-annotated
genome; Lottia's annotation uses different names and returned 0 for everything (ignore).

- **Dispatched** is highest in caenogastropod snails (Littorina 8, Babylonia 8, Pomacea 6),
  3-5 in clams, other snails and the chiton, 1-2 in shell-less or reduced-shell species
  (octopuses 1, Euprymna 1, Aplysia 2, Berghia 2). Expansion with shelled lineages fits
  a role in shell growth, but needs a gene tree to confirm.
- **Diaphanous formins** are 1-3 copies by name. (The name-based hint that left-coiling
  snails have fewer copies does not survive the gene tree below.)
- **Inverted formin (INF2)** is most numerous in clams (Mercenaria 7, oyster 6, Ruditapes 5).
- **Perlucin-like** genes are hugely expanded in clams (69-153) and abalone (98) but this
  name is given to many C-type lectins, which molluscs also expand for immunity, so this
  count is not a clean shell signal.

### Gene trees (replacing name-based counts)
Script: `analysis/gene_trees.py` + `analysis/run_trees.sh`. For the 19 molluscs with usable
annotation plus human and fruit fly, all formin-family genes (295 proteins) and all
dispatched genes (75, with human patched-1 as outgroup), longest isoform per gene;
MAFFT alignment, trimAl, IQ-TREE (best models LG+R6 and LG+R5, 1000 ultrafast bootstraps).
Genes were assigned by membership in clades anchored on human/fly genes.
Outputs: `gene_tree_formin.png`, `gene_tree_dispatched.png`, `gene_tree_*_assignments.csv`,
`gene_tree_counts.csv`.

**Diaphanous** (clade with human DIAPH1-3, fly dia and Lymnaea Lsdia1/2, bootstrap 92):
single-copy in 13 of 19 molluscs, including Littorina (the LG6 gene) and both left-coiling
snails. Duplicates in Pomacea, Babylonia, Gigantopelta, Patella, oyster and O. sinensis.
Lsdia1/2 group with Physella and Biomphalaria (bootstrap 100), as expected for pulmonates.
So copy number does not track coiling direction; the name-based hint was an artefact.

**Dispatched** falls into three groups:

| Group | Human/fly member | Littorina copies | Pattern across molluscs |
|---|---|---|---|
| DISP1/2 (Hedgehog-releasing) | DISP1, DISP2, fly disp | 1 (LG17, 70.0 Mb, in W region) | exactly 1 in every species |
| DISP3 | DISP3 | 1 (LG10) | 0-2 per species |
| Mollusc-specific dispatched-like lineage, sister to DISP1/2 | none | **6, all LG17** (3 at 64.1-64.3 Mb in the W region) | see below |

The mollusc-specific lineage:
- **Present in 12 of 14 species with an external shell** (snails, limpets, abalone, clams,
  chiton; missing in Physella and Gigantopelta).
- **Absent from all 5 shell-less or internal-shell species** (both octopuses, Euprymna,
  Aplysia, Berghia). Taken at face value p ~ 0.002 (Fisher), but the absences reflect
  only 2-3 independent losses (cephalopods; within heterobranch snails), and a missing
  gene can be an annotation gap.
- **Expanded in caenogastropods** by duplications shared by Littorina, Pomacea and Babylonia
  (bootstrap 100 clades), followed by tandem copies within each species: Littorina's 6 sit
  as two groups of 3 on LG17.

So the caenogastropod "dispatched expansion" is real, and it is this mollusc-specific
lineage, not the Hedgehog-releasing DISP1. Against a shell role in adults: Littorina's
copies are barely expressed in the mantle (0-1.4 TPM). The pattern makes the lineage worth
checking in developing shell gland or juvenile mantle data.

## Caveats
- Adding family as a covariate is conservative: it removes between-family genetic
  differences, which lowers LODs. Both versions are reported.
- W is computed from the published widthGrowth by a one-to-one transform, so
  its QTL is the paper's widthGrowth QTL (suggestive there, P = 0.053).
- The genetic map is from F2 crosses, so QTL intervals span many cM and many genes.
