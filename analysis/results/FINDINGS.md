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

**Nautilus, the shelled cephalopod** (`analysis/orthologs/run_nautilus.sh`). Nautilus
pompilius has a genome (GCA_018389105.1, Zhang et al. 2021) but no gene annotation, so the
dispatched and diaphanous proteins were aligned to it with miniprot and the best model per
locus was added to the dispatched tree (`results/trees/dispatched_nautilus.*`).
- Nautilus has one gene in the **mollusc-specific lineage** (982 aa, no stop codons or
  frameshifts; inside the lineage clade, bootstrap 96), plus one DISP1/2 that groups with
  octopus and squid DISP1/2 (bootstrap 100), one diaphanous and one patched.
- So the lineage is present in the one cephalopod with an external shell and absent in the
  three shell-less ones (two octopuses, Euprymna). Within cephalopods its loss coincides
  with loss of the external shell, a cleaner contrast than across classes.
- **Absences checked in the genomes themselves** (`analysis/absence_check.py`,
  `orthologs/run_absence_check.sh` + `run_absence_align.sh`, table `absence_check.csv`).
  All dispatched-family proteins were aligned with miniprot (-N 20 --outs=0.5) to the five
  annotated genomes lacking the lineage, and to periwinkle LG17 as a control (every lineage
  region and DISP1 recovered and correctly labelled; tandem copies merge, so this counts
  presence, not copy number):
  - Octopus bimaculoides: only ~140-aa fragments at ~32% identity (lineage-like and
    DISP3-like alike), not genes; intact DISP1/2 and patched found. **Absent.**
  - Berghia stephanieae (no shell): **absent**; DISP1/2, DISP3, patched found.
  - Gigantopelta aegis (external shell): **absent**; DISP1/2, DISP3, patched found.
  - Physella acuta (external shell): **present but split**, two pieces (480 and 363 aa, 72%
    identity, frameshifts) on small unplaced scaffolds. Missed by the annotation.
  - Aplysia californica (small internal shell): **one partial copy** (668 aa, 45% identity,
    2 frameshifts). Could be a decaying gene or errors in the older AplCal3.0 assembly.
  - Euprymna and Octopus sinensis were not re-checked (annotation only).
  Revised tally: present in 14 of 15 species with an external shell (absent in
  Gigantopelta); in the 5 without one, absent (octopuses, Euprymna, Berghia) or a broken
  partial copy (Aplysia).
- This sits uneasily with the expression evidence (gill- and cilia-biased in bivalves and
  Dreissena larvae). Both could be true, e.g. a role in a ciliated epithelium that shell-less
  lineages no longer need, so the pattern is a lead, not a shell function.

So the caenogastropod "dispatched expansion" is real, and it is this mollusc-specific
lineage, not the Hedgehog-releasing DISP1. Against a shell role in adults: Littorina's
copies are barely expressed in the mantle (0-1.4 TPM). The pattern makes the lineage worth
checking in developing shell gland or juvenile mantle data.

### Expression during development (bivalves)
Script: `analysis/dev_expression.py`. NCBI RefSeq TPM tables for the two bivalves whose
annotation includes larval RNA-seq: Crassostrea virginica (1-day larvae, when the first
shell forms, 23 runs; 22-day juveniles; adult mantle, gill, gonad, hemocytes) and Ruditapes
philippinarum (larvae, umbo larvae, mantle and 13 other tissues). No gastropod in our set
has embryo or larval samples in its RefSeq annotation. Mercenaria ships raw counts only.
Output: `dev_expression_long.csv`, `dev_expression.png`.

Highest TPM of any gene in each group:

| Group | Oyster: larvae (1 d) | Oyster: mantle | Oyster: top tissue | Ruditapes: larva | Ruditapes: mantle | Ruditapes: top tissue |
|---|---|---|---|---|---|---|
| Mollusc-specific dispatched-like | 10.9 | 9.8 | gill 49.6 | 6.1 | 6.3 | gill 41.1 |
| DISP1/2 | 3.9 | 6.2 | hemocytes 22.5 | 6.1 | 1.8 | foot/whole 7-8 |
| Diaphanous | 15.5 | 16.8 | whole oyster 30.3 | 6.0 | 5.2 | whole-organism 10.4 |

- **The mollusc-specific dispatched-like lineage is gill-biased in both bivalves** (and
  hepatopancreas/reproductive in Littorina), with only moderate expression in larvae and
  mantle. That does not look like a shell-forming gene. Its absence in shell-less species
  may track something else they share; the lead is weakened.
- **Diaphanous is expressed at moderate levels in every tissue and stage,** larvae included,
  as expected for a general actin regulator. Its known role (chirality in Lymnaea) is set in
  the first cleavages, earlier than any stage sampled here.

### Developmental time courses (mussel and apple snail)
Script: `analysis/dev_timecourse.py`, figure `dev_timecourse.png`. Expression in CPM
(compare one gene across stages, not genes with each other). Reference curve: mean scaled
profile of shell-forming gene names reaching >= 5 CPM.
- *Mytilus galloprovincialis* (Miglioli et al. 2024, PRJNA996031; NCBI RefSeq counts),
  0-72 h, 15 stages x 2. Added to the dispatched tree: 4 genes in the mollusc-specific
  lineage, 1 DISP1/2, 1 DISP3.
- *Pomacea canaliculata* (Accorsi et al. 2025, GEO GSE240084), stages 2-10 and hatchling, 3-4
  replicates. The stage-to-morphology mapping is not in the deposit.

| Species | Gene (group) | Pattern |
|---|---|---|
| Mussel | LOC143069148, LOC143069150 (mollusc-specific) | ~0.5 CPM to 10-16 CPM between 12 and 16 h, as the shell field forms; stay at 5-10 |
| Mussel | LOC143069149 (mollusc-specific) | switches on earlier, peak 120 CPM at 8 h (gastrula), then 40-80 |
| Mussel | LOC143051730 (mollusc-specific) | maternal (34 CPM in egg), dips, then 15-28 |
| Mussel | DISP1/2, DISP3 | present throughout, no shell-field change |
| Mussel | shell-gene reference | rises steadily from ~16 h to a maximum at 72 h (D-veliger) |
| Apple snail | LOC112571125 (mollusc-specific) | 0.1 to ~2-3 CPM from stage 7 on, rising with the shell-gene reference; other 4 copies near zero |
| Apple snail | DISP1/2 | flat, 35-54 CPM |
| Apple snail | diaphanous LOC112560687 | high throughout (60-115 CPM) |

Reading this:
- In the mussel, two copies of the mollusc-specific lineage switch on exactly as the shell
  field forms, and before the bulk of shell genes. That fits a shell-field role. In the
  snail the lineage is barely expressed, though one copy tracks shell formation.
- Whole-embryo data cannot say which cells express a gene. The 12-16 h window also brings
  the ciliated prototroch, and adult expression is gill-biased in both bivalves (ciliated
  tissue). So "ciliated epithelium gene" fits as well as "shell gene".
- Telling these apart needs spatial data: a single-cell atlas of a trochophore with a
  shell-field cluster (for example Dreissena, GSE192624) or in situ hybridisation.

### Single cells: shell field or ciliated cells?
Script: `analysis/dreissena_singlecell.py`, figure `dreissena_singlecell.png`. Data: Dreissena
rostriformis trochophore (13-15 h), 10x single-cell RNA-seq, GEO GSE192624 (Salamanca-Diaz et
al. 2022), whose authors describe shell-field, prototroch and other ciliated populations.
- The deposited gene models were translated and searched with DIAMOND. By best hit against
  our dispatched set, Dreissena has 4 genes in the mollusc-specific lineage, 1 DISP3 and
  fragments of DISP1/2 (best-hit assignment, not a tree).
- Marker sets were taken as best hits of human proteins: cilia (FOXJ1, 3 tektins, 3 radial
  spoke proteins, axonemal dyneins DNAI1/2/DNALI1, SPAG6, CFAP52) and shell field (engrailed,
  Hox1, BMP2/4, GATA2/3, SMAD1, CHRDL1). Caveolin hits were not detected in any cell.
- 575 cells (>= 300 UMIs), PCA, k-means (k = 7).

Result: the most widely expressed lineage gene (Gene.32991, 113 cells) is concentrated in the
cluster with the strongest cilia signature (40% of its cells express it) and nearly absent
from the two clusters leaning to shell field (0% and 3%). Per cell it correlates with the
cilia score (Spearman rho = 0.36, p = 7e-19) much more than with the shell-field score
(rho = 0.15). Together with gill-biased expression in adult oyster and Ruditapes, this points
to a role in ciliated epithelia rather than shell secretion.

Limits: the shell-field marker set is inferred from human genes and weak in this data; the
clustering is coarse (one large mixed cluster); 3 of the 4 lineage genes are detected in too
few cells to place; one species and one stage. This weakens, but does not rule out, a shell
role in gastropods.

### Do the LG6 candidates behave like shell genes in other molluscs?
Script: `analysis/lg6_orthologs.py` (orthologs by reciprocal best DIAMOND hits,
`analysis/orthologs/run_rbh.sh`), figure `lg6_orthologs.png`, table `lg6_orthologs.csv`.
The 14 mantle-specific genes of the LG6 shape region and the LG6 diaphanous gene were
matched to oyster (adult tissues + 1-day larvae), mussel (0-72 h development) and Dreissena
(larval single cells). Nine known shell genes that are mantle-biased in periwinkle
(carbonic anhydrases, tyrosinase, chitin synthase, nacrein, perlucin) went through the same
tests as controls.

- **Calibration:** the oyster mantle-bias test works: the controls' oyster matches have median
  mantle-bias percentile 94 (8 genes). The mussel time course and the Dreissena cells are not
  consistent for the controls (several carbonic anhydrases are flat or fall in the mussel;
  most controls are too sparse in the single-cell data), so they carry little weight.
- **Most LG6 candidates do not behave like shell genes:** median oyster mantle-bias percentile
  57 (13 genes) vs 94 for controls (Mann-Whitney p = 3e-4).
- **Two exceptions inside the inversion:**
  - **UPAR/Ly6 domain-containing protein crok-like** (9.42 Mb): oyster match at the 92nd percentile
    for mantle bias and the 99th for level in 1-day larvae. Periwinkle has 4 crok-like genes;
    the LG6 one is the most mantle-specific (88 vs 2 TPM), and the oyster gene matches the
    family rather than this copy alone (best reciprocal hit is another periwinkle crok-like
    gene, 139 vs 129 bits). So: a mantle-expressed gene family in both species, with the
    most mantle-specific periwinkle copy inside the LG6 inversion.
  - **rho GTPase-activating protein 6-like** (17.59 Mb, near the inversion edge): reciprocal
    best hits in all three species; oyster mantle bias 89th percentile, larvae 90th; rises
    into the mussel veliger (+1.65 log2).
- Protocadherin-9-like rises ~85-fold into the mussel shell-field window, but its mussel
  match is not reciprocal and its oyster mantle bias is average (57th percentile).
- Diaphanous: average in every test, as before.
- CD109 and murinoglobulin (alpha-macroglobulin family) share one non-reciprocal mussel match,
  so their mussel values are the same gene; not informative.

**Ly6/uPAR family tree** (`gene_trees.py fetch ly6` + `run_trees.sh ly6`; 104 proteins of
80-400 aa named crok/UPAR-Ly6/Ly6-PLAUR/lymphocyte antigen 6, incl. 22 human; WAG+R7;
files in `results/trees/ly6.*`). Names are inconsistent across annotations (Mya, Patella,
Biomphalaria show 0), so only the tree topology is used, not counts.
- The LG6 crok-like gene's closest relative is a Babylonia gene (bootstrap 100): a
  one-to-one caenogastropod ortholog pair.
- Its sister pair is periwinkle's LG2 crok-like gene (LOC138958607) plus another Babylonia
  gene (bootstrap 81), so the LG6 copy arose by a duplication within caenogastropods.
- The wider clade (bootstrap 84-91) also holds octopus, squid, sea slug, chiton and the
  oyster's mantle-biased crok-like gene: the family is old and not limited to shelled species.
- After the duplication the copies diverged in expression: the LG6 copy is mantle-specific
  (88 vs 2 TPM elsewhere, 99.5th percentile) and the LG2 copy is broad (90 vs 70 TPM,
  highest in head). That fits a mantle role taken on by the LG6 copy (an interpretation, not
  a test).

**Shell proteomes** (`analysis/shell_proteome_homologs.py`, `orthologs/run_shellprot.sh`).
Published shell-matrix proteomes (proteins found by mass spectrometry in shell): Lottia
gigantea (Mann et al. 2012, 600 consensus proteins), Crepidula fornicata (Batzel et al. 2022,
185), C. atrasolea (2024, 18), Haliotis laevigata (Mann et al. 2018, 80), Pinctada fucata
(Liu et al. 2015, 43), Cepaea nemoralis (Mann & Jackson 2014, 59 contigs). Searched with DIAMOND
against the periwinkle proteome (e <= 1e-10, best hit per shell protein): 528 periwinkle genes.
- None of the shell proteins hits the three shortlisted genes.
- Lottia's shell proteome has two Ly6-family proteins (ESO89087 "Sleepless/QVR" type, in the
  acid-insoluble consensus; ESP04019 Boudin-like, acid-soluble). Their periwinkle
  counterparts are an LG2 crok-like gene (LOC138958625, 71% identity, 26.42 Mb) and an LG2
  bou-like gene (62%), not the LG6 copy (29%, family-level). LOC138958625 sits beside
  LOC138958607, the LG6 copy's sister in the Ly6 tree. So the LG2 crok-like cluster holds
  the gene whose limpet version is a shell protein, and the LG6 copy is a caenogastropod
  duplicate of that cluster that became mantle-specific. This makes the LG6 copy more
  plausible as a shell-related gene, without showing it is in shell.
- The inversions are not enriched for shell-proteome homologs: LG6 inversion 5 of 384 genes
  (expected 7.3, mostly histones and a ribosomal protein, i.e. cellular carry-over); LG17
  inversion 1 of 194 (expected 3.7).
- Check on the mantle ranking: shell-proteome homologs are 3x over-represented among the
  periwinkle genome's 5% most mantle-biased genes (88 of 514, 17%).
- Formins (diaphanous) are absent from all shell proteomes, as expected for an
  intracellular protein; any effect on coiling would act through development.

**Crab vs Wave expression data:** none usable. A search found no RNA-seq with per-gene
counts for Crab vs Wave or for known LG6 karyotypes. A NimbleGen foot-muscle array (GEO
GSE120697) has no probes on the three shortlisted genes (probe-sequence check). The only
option is Westram et al. 2014 pooled whole-body + embryo RNA-seq (Dryad doi:10.5061/dryad.21pf0,
~27 GB, 2 pools per ecotype per country), not pursued.

Updated candidate shortlist for the LG6 inversion: the crok-like Ly6 gene and rhoGAP6 as
mantle/shell-type genes, plus diaphanous as the coiling-direction candidate.

## Wild snails: inversion karyotypes vs shell shape
Script: `analysis/wild_inversions.py`. Data: Koch et al. 2022 (Evolution 76:2332, Dryad
doi:10.5061/dryad.m905qfv4b): 1,130 snails from three Swedish Crab-Wave transects, with the
same Larsson shape parameters as the F2 cross and karyotypes (0/1/2) for 16 inversions.
Model: trait ~ island + sex + habitat axis (habitatPC1, and its interaction with island)
+ shell length (for shape traits) + karyotype. Without the habitat term, any inversion that
tracks the Crab/Wave habitat would look like a shape locus; the habitat-free effects are
reported alongside and are 1.5-3x larger. Bonferroni threshold over 96 tests: p < 5.2e-4.
A joint model fits all 16 inversions at once (1,053 snails with complete karyotypes).
Outputs: `wild_inversions.csv`, `wild_inversions_joint.csv`, `wild_inversions.png`.

| Inversion | F2 cross (Koch 2021) | Wild snails, one inversion at a time | Wild, joint model |
|---|---|---|---|
| **LGC6.1** (LG6) | S and size (S just below threshold) | size -0.54 SD/copy (p 1e-69); W +0.23 (p 3e-11); height expansion +0.25; opening size +0.16; S -0.15 (p 7e-4) | size -0.25 (p 4e-17), height expansion +0.18 (p 1e-7), W +0.16 (p 5e-6) |
| **LGC17.1** (LG17) | W and D | S +0.22 (p 1e-5); size +0.42; W -0.10 (p 0.015); D +0.11 (p 0.02) | size +0.12 (p 3e-4) only |
| LGC12.2 (LG12) | opening size (suggestive) | height expansion -0.30 (p 4e-23), W -0.22, opening size -0.19 | height expansion, W, opening size all significant |
| LGC14.1 (LG14) | - | size -0.58, W +0.28, S -0.24 | size, W |

What this adds:
- **LG6 is confirmed in the wild as a shell-size inversion, and it also goes with coil
  expansion (W).** In the F2 cross W mapped to LG17; in the wild the LG6 effect on W survives
  adjusting for habitat and for all other inversions.
- **LG17's F2 effect on W and D does not replicate clearly** in the wild (W p = 0.015, not
  significant after correction). LG17 instead goes with opening shape (S) and size, and
  those mostly vanish once the other inversions are in the model.
- **The F2 LG17 QTL lies inside the LGC17.1 inversion**, so that is not the explanation.
  Koch et al. 2021 Table S1 places inv17.1 at 58.3-66.1 cM on the F2 map; through our marker
  placement that is 63.0-70.3 Mb on the 2024 genome, which contains the whole W region
  (63.9-70.3 Mb, peak 65.5 cM). Koch 2021 also found inv17.1 genotype affects width growth in
  the F2 (FDR p = 0.009; their Table S7). Check of the method: their inv6.1/2 bounds
  (0.75-19.65 cM) convert to 1.4-15.0 Mb, matching the flat block. Note that Koch's
  supplementary tables count LG6 from the opposite end to Geno.QTL.csv (their aperture-shape
  peak 18.6 cM = 69.1 - our 50.5 cM); LG17 runs the same way in both.
- The likely reason for the weak wild result: inv17.1 frequency tracks the Crab/Wave
  habitat, so the habitat term absorbs much of its effect. Without the habitat term the wild
  W effect is p = 0.0015. The F2 cross breaks that link (arrangements are shuffled across
  families), which is why it detects the effect. So the LG17 result is "consistent but not
  independently confirmed", rather than "contradicted". The dispatched cluster (64.1-64.3 Mb)
  and DISP1 (70.0 Mb) lie inside the inversion.
- **LG12 replicates**: the F2 cross had a suggestive opening-size QTL there; in the wild
  LGC12.2 is one of the strongest shape inversions.
- Size, W and S are correlated, and inversions are correlated with each other and with the
  habitat cline, so effects on single traits are not independent. Residual habitat
  confounding is possible: habitatPC1 is a coarse summary of a continuous cline.
- Karyotypes test whole inversions, not genes. Recombination is suppressed inside an
  inversion, so wild data cannot narrow LG6 below the inversion itself.

## Step 3: genome species in Raup's morphospace
Script: `analysis/morphospace.py`. Raup W, T and D for 389 museum specimens of 358 living
gastropod species from Noshita, Asami & Ubukata 2012 (Paleobiology 38:322; Dryad
doi:10.5061/dryad.6mk023c0, Supplementary Table 1, parsed from PDF to
`noshita2012_raup_parameters.csv`). Nautilus from Tajika & Klug 2020.
Outputs: `morphospace_genome_species.csv`, `morphospace.png`.

| Genome species | Shape source | W | T | D (Noshita, signed) | Mollusc-specific dispatched | Diaphanous |
|---|---|---|---|---|---|---|
| Pomacea canaliculata | same species | 2.18 | 0.56 | 0.07 | 5 | 2 |
| Physella acuta | same species (Physa acuta) | 1.93 | 2.29 | -0.31 | 0 | 1 |
| Littorina saxatilis | 3 congeners | 2.07 | 1.86 | -0.05 | 6 | 1 |
| Haliotis rufescens | 2 congeners | 14.4 | 0.57 | -0.21 | 2 | 1 |
| Biomphalaria glabrata | Planorbarius (same family) | 1.24 | 0.39 | 0.32 | 1 | 1 |
| Nautilus pompilius | Tajika & Klug 2020 | ~2.9 | 0 | ~0 | no genome | no genome |
| *Lymnaea stagnalis* (reference) | same species | 2.03 | 2.33 | -0.16 | no RefSeq genome | Lsdia1/2 |

- Cross-check: wild L. saxatilis W from the Larsson model (Koch 2022) has median 2.04
  (IQR 1.87-2.22), close to its congeners' Raup-measured 2.07. The two methods agree.
- Only 5 genome species can be placed, and two of them via relatives. That is far too few
  to relate gene counts to geometry; no pattern is claimed. More power needs many more
  species with both a shell measurement and an annotated genome.

## Caveats
- Adding family as a covariate is conservative: it removes between-family genetic
  differences, which lowers LODs. Both versions are reported.
- W is computed from the published widthGrowth by a one-to-one transform, so
  its QTL is the paper's widthGrowth QTL (suggestive there, P = 0.053).
- The genetic map is from F2 crosses, so QTL intervals span many cM and many genes.
