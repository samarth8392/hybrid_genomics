## Introduction

Hybridization and introgression are pervasive evolutionary processes that can reshape genetic variation within and between species, influence adaptation, and alter the trajectory of speciation. Genome-wide datasets increasingly show that hybrid genomes are not uniform mixtures of parental lineages; instead, they are **mosaics** of ancestry blocks whose genomic distribution is determined by recombination, demographic history, and selection (Abbott et al. 2013; Harrison & Larson 2016). Selection can remove introgressed ancestry near loci involved in reproductive isolation (producing ancestry “deserts”) while permitting or favoring introgressed haplotypes that are neutral or adaptive (producing ancestry “islands”). Consequently, *local ancestry* profiles across hybrid genomes provide a direct route to identifying candidate barrier regions, adaptive introgression, and genomic heterogeneity in gene flow.

### Introgression deserts as evolutionary signals

In an admixed population, **introgression deserts** are genomic intervals where ancestry from a donor lineage is strongly depleted relative to the genome-wide baseline. Operationally, deserts can be defined as contiguous segments where local ancestry is confidently near one parental extreme (e.g., near-pure *parent1* ancestry) in individuals or populations that are otherwise admixed. Evolutionarily, deserts are typically interpreted as **regions resisting gene flow**, and can arise through several (non-exclusive) mechanisms:

1. **Selection against introgressed alleles** due to intrinsic incompatibilities (e.g., Dobzhansky–Muller incompatibilities) or disruption of co-adapted gene complexes, yielding low introgression at barrier loci (Harrison & Larson 2016; Ravinet et al. 2017).  
2. **Ecological selection** against alleles that are maladaptive in the recipient environment, which can also generate localized resistance to introgression.  
3. **Linkage and recombination effects**: selection against introgressed alleles is expected to purge larger linked ancestry blocks in low-recombination regions, producing broad deserts even if only a subset of loci are directly selected (Harrison & Larson 2016; Ravinet et al. 2017).  
4. **Structural variation** (e.g., inversions) that reduces effective recombination between ancestry backgrounds, reinforcing restricted introgression in affected regions (Ravinet et al. 2017).  

Because technical factors (marker density, mapping bias, local filtering stringency) and demographic history can also generate heterogeneity, deserts should be treated as **candidate regions** consistent with selection and/or barrier loci, and prioritized for follow-up rather than interpreted as definitive evidence by themselves.

### Empirical examples linking deserts to reproductive barriers and phenotype

Reduced introgression in specific genomic regions—often enriched on sex chromosomes—has been documented across vertebrates and is frequently linked to reproductive isolation and trait maintenance:

- **Mammals (house mice): reproductive isolation and hybrid sterility**. In the European *Mus musculus domesticus* × *M. m. musculus* hybrid zone, the X chromosome contains loci with markedly reduced introgression, consistent with the “large X-effect” and implicating X-linked regions in reproductive isolation (Payseur et al. 2004).  
- **Mammals (humans and primates): deserts of archaic ancestry**. Genome-wide maps of Neanderthal introgression in present-day humans show long intervals (“deserts”) depleted of Neanderthal ancestry, including strong depletion on the X chromosome—patterns consistent with selection against introgressed variation (Sankararaman et al. 2014). In a natural howler monkey hybrid zone (*Alouatta palliata* × *A. pigra*), an X-linked region with reduced introgression overlaps a previously described archaic-ancestry desert on the human X, suggesting recurring genomic features of reproductive isolation across primates (Baiz et al. 2020).  
- **Birds (phenotype-linked restricted introgression): plumage integrity under gene flow**. Hooded and carrion crows exhibit extensive gene flow across most of the genome, but a narrow set of genomic regions shows strong differentiation and limited introgression, associated with the maintenance of plumage phenotypes and assortative mating—an example where restricted introgression aligns with a key phenotype mediating reproductive isolation (Poelstra et al. 2014). In tidal-marsh sparrows, steeper clines and reduced introgression at markers linked to tidal-marsh adaptations (including melanin-based pigmentation) support selection maintaining locally adaptive phenotypes despite extensive admixture (Walsh et al. 2016).  
- **Birds (sex-chromosome barriers)**. In a spotted eagle hybrid zone (*Aquila clanga* × *A. pomarina*), introgression was reduced on Z-linked loci relative to autosomes, consistent with a disproportionate role of sex-linked loci in barrier formation (Väli et al. 2011).  
- **Reptiles (sexual traits and barrier loci)**. In wall lizards (*Podarcis muralis*), exaggerated male morphology and coloration are associated with reduced and asymmetric gene flow in secondary contact, and genomic analyses identify regions of consistently low exchange (e.g., around *ATXN1*), implicating behavioral and sexual traits in reproductive isolation (Yang et al. 2020).  
- **Reptiles (functional trait loci resisting introgression)**. In a rattlesnake hybrid zone, toxin-gene presence is tightly associated with venom phenotype, yet toxin genes do not introgress broadly beyond the hybrid zone, indicating that even with extensive hybridization, strong selection can restrict introgression of functionally important loci (Zancolli et al. 2016). In musk turtles, selection can inhibit introgression in portions of the genome while promoting it elsewhere, highlighting how ecological and intrinsic processes jointly shape genomic permeability (Scott et al. 2019).  

Together, these systems motivate using introgression deserts to prioritize candidate genomic regions likely shaped by selection, reproductive barriers, or ecologically important trait variation.

### Hybrid genomes as tools for functional genomics

Hybrid populations provide a powerful substrate for **functional genomics** because introgressed haplotypes act as naturally occurring genetic perturbations segregating across many recombinant individuals. When local ancestry is quantified along the genome, hybrid systems can be leveraged to: (i) associate ancestry tracts with phenotypes (admixture mapping), (ii) prioritize candidate genes within islands/deserts, and (iii) integrate ancestry with molecular traits (e.g., gene expression, splicing, chromatin accessibility) to identify lineage-specific regulatory and coding effects. In this framework, ancestry mosaics can guide targeted functional assays by selecting individuals carrying contrasting ancestry at specific regions and nominating candidate loci for mechanistic follow-up.

### Study system and objectives

The Eastern Massasauga Rattlesnake (*Sistrurus catenatus*) and Western Massasauga (*S. tergeminus*) provide a compelling vertebrate system for investigating the genomic consequences of historical introgression. A distinct Iowa lineage of *S. catenatus* shows genomic signatures consistent with historical admixture with *S. tergeminus* (Sovic et al. 2016). However, the local genomic architecture of ancestry in Iowa individuals—particularly the locations and sizes of ancestry blocks that resist or permit introgression—remains incompletely characterized.

Here, we quantify **sliding-window hybrid index** (the proportion of ancestry derived from *S. tergeminus*) across the genome using a diagnostic-marker likelihood framework (Buerkle 2005). We then identify candidate **introgression deserts**—genomic regions retaining *S. catenatus* ancestry in otherwise admixed genomes—using confidence-interval–based classification and interval-calling procedures. These results provide a foundation for evolutionary inference (barrier loci, heterogeneous gene flow) and for prioritizing candidate genomic regions for downstream functional genomic analyses.

---

## Materials and Methods

### Variant filtering and sample definition

Hybrid ancestry was estimated from SNPs in VCF format using a shell workflow (`hybrid_index_sliding_window.sh`) and a Python estimator implementing a diagnostic-marker hybrid index model (`calculate_hybrid_index_ml.py`). Analyses used three sample sets:

- **Focal Iowa individuals** (putative hybrids; analyzed for local ancestry),
- **Parent 1**: *S. catenatus* reference samples (**excluding** Iowa individuals),
- **Parent 2**: *S. tergeminus* reference samples.

VCFs were filtered to retain **biallelic SNPs** and exclude indels. Sites were filtered for missingness using a maximum missingness threshold (e.g., requiring ≥80% non-missing genotypes across retained samples). A global MAF filter was not applied by default, because alleles fixed (or near-fixed) in the smaller parental group can be rare in the combined dataset and would be removed by a global MAF threshold.

### Identification of ancestry-informative markers (AIMs)

We identified ancestry-informative markers using Weir and Cockerham’s $F_{ST}$ between the two parental reference populations, computed with VCFtools (`--weir-fst-pop`). SNPs with high differentiation were retained as candidate AIMs using a stringent threshold (pipeline parameter `--fst-threshold`). The filtered VCF was then subset to this AIM set for all retained samples (focal individuals + both parental references).

### Diagnostic marker definition from parental allele frequencies

Within the AIM set, loci were further restricted to **near-diagnostic** markers based on parental allele frequencies computed directly from the VCF within `calculate_hybrid_index_ml.py`. For each AIM, the ALT allele frequency was estimated separately in parent1 and parent2 using diploid, non-missing genotypes only. A locus was considered diagnostic if it approximated a fixed difference between taxa, defined by a near-fixation threshold $\tau$:

- parent1 ALT frequency $\le \tau$ and parent2 ALT frequency $\ge 1-\tau$, **or**
- parent1 ALT frequency $\ge 1-\tau$ and parent2 ALT frequency $\le \tau$,

with $\tau = 0.05$ by default (`--diagnostic-threshold 0.05`). Loci were additionally required to have at least a minimum number of called parental genotypes per population (`--min-parent-called`).

For each diagnostic locus, genotypes in focal individuals were oriented to represent the **count of parent2 alleles** ($0,1,2$) regardless of whether the ALT allele was fixed in parent2 or fixed in parent1. Sites with missing or partially missing genotypes were excluded for that individual and window.

### Sliding-window hybrid index estimation

Local ancestry was summarized in sliding windows across each chromosome. Windows were defined by size $W$ (default 50 kb) and step $S$ (default 10 kb) (`--window`, `--step`). Only windows containing at least $M_{\min}$ diagnostic loci with called diploid genotypes for a given individual were retained (default $M_{\min}=10$; `--min-aims`) to reduce variance in low-information regions.

Hybrid index $h$ was defined as the proportion of ancestry derived from parent2 (*S. tergeminus*), with $h=0$ corresponding to pure parent1 ancestry and $h=1$ corresponding to pure parent2 ancestry. Under the diagnostic-marker (fixed-difference) model, each locus contributes a binomial likelihood for the observed parent2 allele count. Assuming loci are independent (for simplicity), the maximum-likelihood estimate within a window reduces to the observed parent2-allele fraction:

$$
\hat{h} \;=\; \frac{\sum_{\ell=1}^{M} n_{\ell,\mathrm{p2}}}{2M},
$$

where $n_{\ell,\mathrm{p2}}\in\{0,1,2\}$ is the parent2 allele count at diagnostic locus $\ell$ and $M$ is the number of diagnostic loci with called genotypes in the window for that individual.

### Mitigating parental sample-size imbalance via subsampling permutations

Because parental sample sizes were unbalanced (many parent1 individuals vs. fewer parent2 individuals), we incorporated a permutation-based subsampling procedure. For each replicate, a random subset of parent1 individuals of size $k$ (e.g., $k=10$) was drawn without replacement (`--parent1-subsample 10`) and used to compute parent1 allele frequencies and diagnostic status. Parent2 allele frequencies were computed from all available parent2 samples. Hybrid index was then estimated per window for each focal individual. This procedure was repeated $R$ times (default $R=100$; `--n-reps 100`).

For each sample × window, the final hybrid index estimate was summarized as:
- **Mean hybrid index** across replicates: $\overline{h}$,
- **Empirical 95% interval** from the 2.5th and 97.5th percentiles of $\hat{h}$ across replicates,
- **Number of replicates used** (windows can fail in some replicates if they fall below $M_{\min}$ due to replicate-specific diagnostic filtering).

The output table (`hybrid_index_windows.tsv`) contains, per sample and window: chromosome, window start/end, mean number of loci, mean hybrid index, CI bounds, a categorical ancestry label, and number of replicates contributing to the summary.

### CI-based ancestry classification for windows

For downstream region calling and visualization, each window was classified using CI bounds:
- **parent1** if $\mathrm{CI}_{\mathrm{high}} \le t_1$,
- **parent2** if $\mathrm{CI}_{\mathrm{low}} \ge t_2$,
- **admixed** otherwise,

with $t_1=0.1$ and $t_2=0.9$ by default. This conservative rule requires that the uncertainty interval lies entirely within a parental range before labeling a window as “pure-parent” ancestry.

### Identification of introgression deserts

Candidate **introgression deserts** (regions resistant to parent2 ancestry) were called from `hybrid_index_windows.tsv` using a dedicated Python script (`identify_introgression_deserts.py`) that operates on the CI-based window classifications.

**Per-sample deserts.** For each focal Iowa individual and chromosome:
1. Select windows classified as **parent1** (i.e., $\mathrm{CI}_{\mathrm{high}} \le t_1$).
2. Merge consecutive parent1 windows into intervals if their genomic coordinates are adjacent or separated by at most a small gap (default gap $\le S$).
3. Retain merged intervals exceeding a minimum physical length $L_{\min}$ (user-defined; e.g., 100–250 kb).

**Shared deserts across Iowa individuals.** To identify deserts plausibly shaped by selection at the population level, the script additionally computes, for each window, the proportion of Iowa individuals classified as parent1. Intervals are merged and retained if:
- the parent1 proportion exceeds a threshold (e.g., $\ge 0.8$ of individuals), and
- the merged interval exceeds $L_{\min}$.

These “shared deserts” highlight regions where Iowa genomes consistently retain *S. catenatus* ancestry in an otherwise admixed background.

### Gene annotation of deserts

To prioritize candidate genes for follow-up, desert intervals were intersected with the genome annotation (GFF). For each desert, the script reports:
- genes overlapping the interval, and
- the nearest upstream and downstream genes within a configurable flank distance (e.g., ±50 kb).

This yields a candidate gene set for desert regions suitable for downstream evolutionary and functional genomic analyses.

### Visualization

Hybrid index profiles were visualized per chromosome using `plot_hybrid_index.py`, which generates: (i) a mean profile across focal individuals (with uncertainty bands) and (ii) individual-level profiles. Plots include background shading for parent1/admixed/parent2 ranges and annotations for candidate extreme segments for integration with gene models.

---

## References

Abbott, R., Albach, D., Ansell, S., Arntzen, J. W., Baird, S. J. E., Bierne, N., et al. (2013). Hybridization and speciation. *Journal of Evolutionary Biology*, 26, 229–246.

Baiz, M. D., Tucker, P. K., Mueller, J. L., & Cortés-Ortiz, L. (2020). X-Linked Signature of Reproductive Isolation in Humans is Mirrored in a Howler Monkey Hybrid Zone. *Journal of Heredity*, 111, 419–428.

Buerkle, C. A. (2005). Maximum-likelihood estimation of a hybrid index based on molecular markers. *Molecular Ecology Notes*, 5, 684–687.

Harrison, R. G., & Larson, E. L. (2016). Heterogeneous genome divergence, differential introgression, and the origin and structure of hybrid zones. *Molecular Ecology*, 25, 2454–2466.

Payseur, B. A., Krenz, J. G., & Nachman, M. W. (2004). Differential patterns of introgression across the X chromosome in a hybrid zone between two species of house mice. *Evolution*, 58, 2064–2078.

Poelstra, J. W., Vijay, N., Hoeppner, M. P., & Wolf, J. B. W. (2014). The genomic landscape underlying phenotypic integrity in the face of gene flow in crows. *Science*, 344, 1410–1414.

Ravinet, M., Faria, R., Butlin, R. K., Galindo, J., Bierne, N., Rafajlović, M., et al. (2017). Interpreting the genomic landscape of speciation: a road map for finding barriers to gene flow. *Journal of Evolutionary Biology*, 30, 1450–1477.

Sankararaman, S., Mallick, S., Dannemann, M., Prüfer, K., Kelso, J., Pääbo, S., Reich, D., & Patterson, N. (2014). The genomic landscape of Neanderthal ancestry in present-day humans. *Nature*, 507, 354–357.

Scott, P. A., Glenn, T. C., & Rissler, L. J. (2019). Formation of a recent hybrid zone offers insight into the geographic puzzle and maintenance of species boundaries in musk turtles. *Molecular Ecology*, 28, 761–771.

Sovic, M. G., Fries, A. C., & Gibbs, H. L. (2016). Origin of a cryptic lineage in a threatened reptile through isolation and historical hybridization. *Heredity*, 117, 358–366.

Väli, Ü., Saag, L., Dombrovski, V., Meyburg, B.-U., Maciorowski, G., Mizera, T., et al. (2011). Sex- and species-biased gene flow in a spotted eagle hybrid zone. *BMC Evolutionary Biology*, 11, 100.

Walsh, J., Shriver, W. G., Olsen, B. J., & Kovach, A. I. (2016). Differential introgression and the maintenance of species boundaries in an advanced generation avian hybrid zone. *BMC Evolutionary Biology*, 16, 65.

Yang, W., Feiner, N., & Uller, T. (2020). Spatial variation in gene flow across a hybrid zone reveals causes of reproductive isolation and asymmetric introgression in wall lizards. *Evolution*, 74, 1289–1300.

Zancolli, G., Baker, T. G., Barlow, A., Bradley, R. K., Calvete, J. J., Carter, K. C., et al. (2016). Is Hybridization a Source of Adaptive Venom Variation in Rattlesnakes? A Test, Using a *Crotalus scutulatus* × *viridis* Hybrid Zone in Southwestern New Mexico. *Toxins*, 8, 188.

| Scatenatus gene_id | Ortholog / best match (mammal)             | Protein function (summary)                                                                                                | Source       |
| ------------------ | ------------------------------------------ | ------------------------------------------------------------------------------------------------------------------------- | ------------ |
| Scate06229         | OPRK1                                      | Kappa opioid receptor (GPCR); mediates opioid peptide signaling in excitable tissues                                      |              |
| Scate06230         | RB1CC1 (FIP200)                            | Scaffold regulating cell growth/proliferation, apoptosis, autophagy, and migration                                        | ([NCBI][1])  |
| Scate06149         | HNF4G                                      | Nuclear receptor transcription factor; activates RNA Pol II transcription programs                                        | ([NCBI][2])  |
| Scate10612         | CLSTN2 (mouse evidence)                    | Predicted Ca2+-binding / synapse-related roles (calsyntenin family; neuronal processes)                                   |              |
| Scate10613         | (unresolved; “L345_15459”)                 | insufficient data (uncharacterized in provided annotation)                                                                | —            |
| Scate11662         | ATF7IP                                     | Heterochromatin-associated nuclear regulator; can coactivate or corepress transcription depending on partners             |              |
| Scate11663         | PLB1                                       | Membrane phospholipase with lysophospholipase/phospholipase A2 activity; broader lipid hydrolase activity                 |              |
| Scate11596         | NUDT4 (DIPP2)                              | Hydrolase controlling turnover of diphosphoinositol polyphosphates (high-energy signaling metabolites)                    | ([NCBI][3])  |
| Scate11597         | UBE2N                                      | E2 ubiquitin-conjugating enzyme; implicated in DNA postreplication repair (mouse evidence cited in record)                |              |
| Scate11598         | MRPL42                                     | Nuclear-encoded mitochondrial ribosomal protein; supports mitochondrial translation                                       |              |
| Scate11599         | SOCS2                                      | Cytokine-inducible negative regulator of JAK/STAT signaling; interacts with signaling complexes                           | ([NCBI][4])  |
| Scate11562         | SYT1                                       | Synaptic vesicle Ca2+ sensor triggering neurotransmitter release during exocytosis                                        | ([NCBI][5])  |
| Scate11546         | KCNC2                                      | Voltage-gated K+ channel (delayed rectifier); enables rapid repolarization in excitable membranes                         | ([NCBI][6])  |
| Scate11547         | CAPS2 (calcyphosine 2)                     | Ca2+-binding protein with EF-hand motifs                                                                                  | ([NCBI][7])  |
| Scate11548         | CAPS2 (calcyphosine 2)                     | Ca2+-binding protein with EF-hand motifs                                                                                  | ([NCBI][7])  |
| Scate11549         | LOC101067630 (uncharacterized)             | insufficient data (uncharacterized in provided annotation)                                                                | —            |
| Scate11550         | GLIPR1                                     | CAP/PR-1 superfamily protein; implicated in immune defense/cancer biology; reported pro-apoptotic in some contexts        | ([NCBI][8])  |
| Scate11551         | (unresolved; “uncharacterized protein”)    | insufficient data (uncharacterized in provided annotation)                                                                | —            |
| Scate11552         | KRR1 (mouse/rat evidence)                  | SSU processome / ribosome biogenesis factor; rRNA processing / small-subunit biogenesis                                   |              |
| Scate11553         | PHLDA1                                     | Conserved nuclear protein; may contribute to anti-apoptotic effects of IGF-1                                              |              |
| Scate11554         | NAP1L1                                     | Histone chaperone (NAP family); involved in DNA replication/chromatin formation and proliferation regulation              | ([NCBI][9])  |
| Scate11555         | OSBPL8 (ORP8)                              | Lipid-binding/transfer protein; exchanges PS/PI4P/oxysterols between ER and plasma membrane                               | ([NCBI][10]) |
| Scate11556         | ZDHHC17 (HIP14)                            | Palmitoyltransferase (S-acylation); impacts protein palmitoylation and related signaling/cell-death programs              | ([NCBI][11]) |
| Scate11557         | ZDHHC17 (HIP14)                            | Palmitoyltransferase (S-acylation); impacts protein palmitoylation and related signaling/cell-death programs              | ([NCBI][11]) |
| Scate11558         | CSRP2 (mouse evidence)                     | LIM-domain protein; predicted actinin binding / muscle structural roles; muscle development/sarcomere organization        |              |
| Scate11559         | E2F7                                       | Transcriptional repressor in DNA damage/p53 signaling and transcription regulation; roles reported in angiogenesis        | ([NCBI][12]) |
| Scate11560         | NAV3 (mouse evidence)                      | Predicted microtubule-binding; roles in migration and microtubule dynamics regulation                                     |              |
| Scate11561         | NAV3 (mouse evidence)                      | Predicted microtubule-binding; roles in migration and microtubule dynamics regulation                                     |              |
| Scate11459         | ING-family (paralog ambiguous)             | ING proteins: chromatin-regulatory (PHD-finger) tumor-suppressor family; DNA repair/apoptosis regulation                  |              |
| Scate11460         | TSPAN12                                    | Tetraspanin family membrane protein; mediates signaling events affecting development/growth/motility                      | ([NCBI][13]) |
| Scate11461         | KCND2                                      | Voltage-gated K+ channel subunit; regulates excitability via voltage-dependent K+ permeability                            |              |
| Scate12346         | PCDH18                                     | Protocadherin (cadherin superfamily); implicated in neuronal cell-cell connectivity (specific function not fully defined) | ([NCBI][14]) |
| Scate12348         | (unresolved; “L345_05116”)                 | insufficient data (uncharacterized in provided annotation)                                                                | —            |
| Scate12349         | SCLT1                                      | Adaptor linking clathrin↔sodium channel components; also distal-appendage centriole protein needed for ciliogenesis       |              |
| Scate12350         | JADE1 (PHF17)                              | Transcriptional coactivator; part of histone acetyltransferase complexes; implicated in Wnt pathway/cell-cycle control    | ([NCBI][15]) |
| Scate12351         | PGRMC2                                     | Heme-binding membrane component; implicated in adipose development; ER/nuclear-envelope localization                      | ([NCBI][16]) |
| Scate12352         | LARP1B                                     | La-motif RNA-binding protein family member; likely RNA metabolism/regulation; multiple isoforms                           | ([NCBI][17]) |
| Scate12353         | (unresolved; “L345_05115” uncharacterized) | insufficient data (uncharacterized in provided annotation)                                                                | —            |
| Scate12354         | MFSD8 (mouse evidence)                     | Lysosomal membrane protein; implicated in lysosome/TORC1/autophagy-related processes; ortholog linked to CLN7             |              |
| Scate12355         | PLK4                                       | Ser/Thr kinase at centrioles; master regulator of centriole duplication in cell cycle                                     |              |
| Scate12356         | HSPA4L                                     | Heat-shock inducible Hsp70-family chaperone; protects against aggregated proteins under stress                            | ([NCBI][18]) |
| Scate12357         | SLC25A31 (ANT4)                            | Mitochondrial ADP/ATP carrier; implicated in spermatogenesis/flagellar cytoskeleton association                           |              |
| Scate12358         | INTU (aka PDZD6)                           | Ciliogenesis / planar cell polarity protein; involved in craniofacial/digit morphogenesis; basal-body/cilium localization | ([NCBI][19]) |

[1]: https://www.ncbi.nlm.nih.gov/gene/9821?utm_source=chatgpt.com "RB1CC1 RB1 inducible coiled-coil 1 [Homo sapiens (human)] - Gene - NCBI"
[2]: https://www.ncbi.nlm.nih.gov/gene/3174 "HNF4G hepatocyte nuclear factor 4 gamma [Homo sapiens (human)] - Gene - NCBI"
[3]: https://www.ncbi.nlm.nih.gov/gene/11163?utm_source=chatgpt.com "NUDT4 nudix hydrolase 4 [Homo sapiens (human)] - Gene - NCBI"
[4]: https://www.ncbi.nlm.nih.gov/gene/8835?utm_source=chatgpt.com "SOCS2 suppressor of cytokine signaling 2 [Homo sapiens (human)] - Gene - NCBI"
[5]: https://www.ncbi.nlm.nih.gov/gene/6857?utm_source=chatgpt.com "SYT1 synaptotagmin 1 [Homo sapiens (human)] - Gene - NCBI"
[6]: https://www.ncbi.nlm.nih.gov/gene/3747?utm_source=chatgpt.com "KCNC2 potassium voltage-gated channel subfamily C member 2 [Homo sapiens (human)] - Gene - NCBI"
[7]: https://www.ncbi.nlm.nih.gov/gene/84698?utm_source=chatgpt.com "CAPS2 calcyphosine 2 [Homo sapiens (human)] - Gene - NCBI"
[8]: https://www.ncbi.nlm.nih.gov/gene/11010?utm_source=chatgpt.com "GLIPR1 GLI pathogenesis related 1 [Homo sapiens (human)] - Gene - NCBI"
[9]: https://www.ncbi.nlm.nih.gov/gene/4673?utm_source=chatgpt.com "NAP1L1 nucleosome assembly protein 1 like 1 [Homo sapiens (human)] - Gene - NCBI"
[10]: https://www.ncbi.nlm.nih.gov/gtr/genes/114882/?utm_source=chatgpt.com "OSBPL8 oxysterol binding protein like 8 - NIH Genetic Testing Registry (GTR) - NCBI"
[11]: https://www.ncbi.nlm.nih.gov/gene/23390?utm_source=chatgpt.com "ZDHHC17 zDHHC palmitoyltransferase 17 [Homo sapiens (human)] - Gene - NCBI"
[12]: https://www.ncbi.nlm.nih.gov/gene/144455?utm_source=chatgpt.com "E2F7 E2F transcription factor 7 [Homo sapiens (human)] - Gene - NCBI"
[13]: https://www.ncbi.nlm.nih.gov/gene/23554?utm_source=chatgpt.com "TSPAN12 tetraspanin 12 [Homo sapiens (human)] - Gene - NCBI"
[14]: https://www.ncbi.nlm.nih.gov/gene/54510?utm_source=chatgpt.com "PCDH18 protocadherin 18 [Homo sapiens (human)] - Gene - NCBI"
[15]: https://www.ncbi.nlm.nih.gov/gene/79960?utm_source=chatgpt.com "JADE1 jade family PHD finger 1 [Homo sapiens (human)] - Gene - NCBI"
[16]: https://www.ncbi.nlm.nih.gov/gene/10424?utm_source=chatgpt.com "PGRMC2 progesterone receptor membrane component 2 [Homo sapiens (human)] - Gene - NCBI"
[17]: https://www.ncbi.nlm.nih.gov/gene/55132?utm_source=chatgpt.com "LARP1B La ribonucleoprotein 1B [Homo sapiens (human)] - Gene - NCBI"
[18]: https://www.ncbi.nlm.nih.gov/gene/22824?utm_source=chatgpt.com "HSPA4L heat shock protein family A (Hsp70) member 4 like [Homo sapiens (human)] - Gene - NCBI"
[19]: https://www.ncbi.nlm.nih.gov/gene/27152?utm_source=chatgpt.com "INTU inturned planar cell polarity protein [Homo sapiens (human)] - Gene - NCBI"
