## Introduction

Hybridization and introgression are pervasive evolutionary processes that can reshape genetic variation within and between species, influence adaptation, and alter the trajectory of speciation. Genome-wide data increasingly show that hybrid genomes are not uniform mixtures of parental lineages; instead, they are **mosaics** of ancestry blocks whose size and genomic distribution are determined by recombination, demographic history, and selection. Selection can remove introgressed ancestry near loci involved in reproductive isolation (producing ancestry “deserts”) while permitting or favoring introgressed haplotypes that are neutral or adaptive in the recipient background (producing ancestry “islands”). Consequently, *local ancestry* profiles across hybrid genomes provide a direct route to identifying candidate barrier regions, potential adaptive introgression, and genomic heterogeneity in gene flow.

Hybrid genomes are also powerful substrates for **functional genomics**. Introgressed haplotypes act as naturally occurring genetic perturbations segregating across many recombinant individuals. When ancestry is quantified along the genome, hybrid populations can be leveraged to: (i) associate ancestry tracts with phenotypes (admixture mapping), (ii) prioritize candidate genes within introgressed segments, and (iii) integrate ancestry with molecular traits (e.g., gene expression, splicing, chromatin accessibility) to reveal regulatory and coding effects that differ between parental lineages. In this framework, windowed ancestry estimates can be used to select individuals carrying contrasting ancestry at specific regions for targeted transcriptomic or regulatory assays, and to nominate genomic segments for mechanistic follow-up.

The Eastern Massasauga Rattlesnake (*Sistrurus catenatus*) and Western Massasauga (*S. tergeminus*) provide a compelling system for investigating the genomic consequences of historical introgression. Prior work has identified an Iowa lineage of *S. catenatus* with genomic signatures of admixture with *S. tergeminus*, consistent with secondary contact and introgression during post-glacial range shifts. However, the **local genomic architecture of ancestry** in Iowa individuals remains incompletely characterized. Key questions include: (1) how hybrid ancestry varies among individuals, (2) whether introgression is uniformly distributed or concentrated in particular genomic regions, and (3) whether ancestry patterns suggest selection for or against specific introgressed segments.

Here, we quantify **sliding-window hybrid index** (the proportion of ancestry derived from *S. tergeminus*) across the genome using a likelihood framework for diagnostic markers. We identify ancestry-informative markers (AIMs) using high differentiation between parental reference populations and estimate hybrid index in windows along each chromosome. This produces a genome-wide ancestry mosaic for each Iowa individual, providing a foundation for evolutionary inference and for prioritizing candidate genomic regions for downstream functional genomic analyses.

---

## Materials and Methods

### Variant filtering and sample definition

Hybrid ancestry was estimated from SNPs in VCF format using a custom pipeline (`hybrid_index_sliding_window.sh`) built around **VCFtools** and a Python estimator (`calculate_hybrid_index_ml.py`). Analyses were restricted to a predefined set of samples comprising:

- **Focal Iowa individuals** (putative hybrids; “input” samples),
- **Parent 1**: *S. catenatus* reference samples (excluding Iowa samples),
- **Parent 2**: *S. tergeminus* reference samples.

Prior to ancestry estimation, the VCF was filtered to retain **biallelic SNPs** only (explicit enforcement of min/max alleles = 2) and to remove indels. Sites were filtered for missingness using a maximum missingness threshold of 20% (i.e., requiring at least 80% genotypes called across retained samples). A **global MAF filter was not applied by default**, because alleles fixed (or near-fixed) in the smaller parental group can be rare in the combined dataset and would be removed by a global MAF threshold.

### Identification of ancestry-informative markers (AIMs)

To identify AIMs, we calculated Weir and Cockerham’s $F_{ST}$ between the two parental reference populations using VCFtools (`--weir-fst-pop`). SNPs with high interspecific differentiation were retained as candidate AIMs using a stringent threshold (parameterized as `--fst-threshold` in the pipeline). The filtered VCF was then subset to this AIM set for all retained samples (focal individuals + both parental references).

### Diagnostic marker definition

Within the AIM set, loci were further restricted to **near-diagnostic** markers based on parental allele frequencies computed directly from the VCF inside `calculate_hybrid_index_ml.py` (avoiding external frequency tables and eliminating coordinate misalignment risks). For each AIM, the ALT allele frequency was estimated separately in parent1 and parent2 using **diploid, non-missing genotypes only**. A locus was considered diagnostic if it approximated a fixed difference between taxa, defined by a near-fixation threshold:

- parent1 ALT frequency $\le \tau$ and parent2 ALT frequency $\ge 1-\tau$, **or**
- parent1 ALT frequency $\ge 1-\tau$ and parent2 ALT frequency $\le \tau$,

with $\tau = 0.05$ by default (`--diagnostic-threshold 0.05`). Loci were additionally required to have a minimum number of called parental genotypes per population (`--min-parent-called`; default $\ge 1$, adjustable for conservatism).

For each diagnostic locus, genotypes in focal individuals were re-oriented to represent the **count of parent2 alleles** ($0,1,2$) regardless of whether the ALT allele was fixed in parent2 or fixed in parent1. Sites with missing or partially missing genotypes (e.g., `0/.`, `./1`) were excluded for that individual.

### Sliding-window hybrid index estimation

Local ancestry was summarized in sliding windows across each chromosome. Windows were defined by size $W$ (default 50 kb) and step $S$ (default 10 kb) (`--window`, `--step`). Only windows containing at least $M_{\min}$ diagnostic loci with called diploid genotypes for a given individual were retained (default $M_{\min}=10$; `--min-aims`) to reduce variance in low-information regions.

Hybrid index $h$ was defined as the proportion of ancestry derived from parent2 (*S. tergeminus*), with $h=0$ corresponding to pure parent1 ancestry and $h=1$ corresponding to pure parent2 ancestry. Under the diagnostic-marker (fixed-difference) model, each diagnostic locus contributes a binomial likelihood for the parent2 allele count. Treating loci as independent (assumption used for simplicity), the maximum-likelihood estimate within a window reduces to the observed parent2-allele fraction:

$$
\hat{h} \;=\; \frac{\sum_{\ell=1}^{M} n_{\ell,\mathrm{p2}}}{2M},
$$

where $n_{\ell,\mathrm{p2}}\in\{0,1,2\}$ is the parent2 allele count at diagnostic locus $\ell$ and $M$ is the number of diagnostic loci with called genotypes in the window for that individual.

### Mitigating parental sample-size imbalance via subsampling permutations

To reduce sensitivity to unbalanced parental sample sizes (many parent1 individuals vs. fewer parent2 individuals), we incorporated a permutation-based subsampling procedure in the Python estimator. For each replicate, a random subset of parent1 individuals of size $k$ (e.g., $k=10$) was drawn without replacement (`--parent1-subsample 10`) and used to compute parent1 allele frequencies and diagnostic status. Parent2 allele frequencies were computed from all available parent2 samples. Hybrid index was then estimated per window for each focal individual. This procedure was repeated $R$ times (default $R=100$; `--n-reps 100`).

For each sample × window, the final hybrid index estimate was summarized as:

- **Mean hybrid index** across replicates: $\overline{h}$,
- **Permutation interval** from the empirical 2.5th and 97.5th percentiles of $\hat{h}$ across replicates,
- The **number of replicates used** (windows can fail in some replicates if they fall below $M_{\min}$ due to diagnostic filtering).

This permutation interval reflects variability induced by parental subsampling (and the induced diagnostic-site set), under the locus-independence simplification.

### Ancestry classification

For descriptive purposes, windows were classified into ancestry categories based on estimated hybrid index values:

- $\hat{h} < 0.1$: parent1 ancestry,
- $\hat{h} > 0.9$: parent2 ancestry,
- otherwise: admixed.

These thresholds were used for visualization and summary statistics only and were not part of the likelihood model.
