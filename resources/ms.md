Introduction
Hybridization and introgression are fundamental evolutionary processes that shape genetic diversity within and between species (Anderson, 1949; Barton & Hewitt, 1985; Harrison & Larson, 2016). While historically viewed as rare or evolutionarily insignificant, genomic analyses have revealed that introgressive hybridization is widespread across the tree of life and can have profound consequences for adaptation, speciation, and extinction (Mallet, 2005; Abbott et al., 2013; Suarez-Gonzalez et al., 2018). The genomic revolution has transformed our ability to detect historical gene flow and characterize its evolutionary outcomes, shifting focus from identifying whether hybridization occurred to understanding the genomic architecture and functional consequences of introgression (Gompert & Buerkle, 2009; Payseur & Rieseberg, 2016).
A central challenge in hybridization research is quantifying individual ancestry and identifying which genomic regions have experienced introgression. Hybrid index—the proportion of an individual's genome derived from a particular parental population—provides a fundamental metric for characterizing admixture levels (Buerkle, 2005). Traditional approaches estimated hybrid index from small numbers of genetic markers, limiting resolution and statistical power (Boecklen & Howard, 1997). Modern genomic methods enable genome-wide ancestry inference at fine spatial scales through local ancestry inference, revealing that hybrid genomes are often mosaics of parental ancestry rather than uniform blends (Gompert et al., 2012; Brandvain et al., 2014). These heterogeneous ancestry landscapes result from the interplay of selection, recombination, and drift acting on introgressed genomic regions (Harrison & Larson, 2016).
The Eastern Massasauga Rattlesnake (Sistrurus catenatus) and Western Massasauga (S. tergeminus) provide a compelling system for investigating the genomic consequences of historical hybridization. These sister species diverged approximately 2-3 million years ago and currently occupy largely parapatric distributions across central North America (Kubatko et al., 2011; Wooten & Gibbs, 2012). Sovic et al. (2016) identified a phylogenetically distinct lineage of S. catenatus in Iowa that exhibits genomic signatures of admixture with S. tergeminus. Demographic modeling supported a two-step process: initial divergence in allopatry approximately 34,000 years ago, followed by secondary contact and introgression approximately 11,000 years ago during post-glacial range shifts. This historical introgression occurred on evolutionary timescales and predates modern habitat fragmentation, representing natural hybridization rather than anthropogenic genetic swamping (Sovic et al., 2016).
Despite evidence for historical introgression, the genomic architecture of ancestry in Iowa rattlesnakes remains uncharacterized. Key questions include: (1) What is the extent of genomic variation in ancestry among individual hybrid snakes? (2) Are introgressed regions distributed uniformly across the genome or concentrated in specific chromosomal regions? (3) Do ancestry patterns suggest selection against or for particular introgressed genomic segments? Addressing these questions requires local ancestry inference methods that can identify parental origins of genomic segments at fine spatial scales.
Here, we apply sliding-window hybrid index analysis to characterize local ancestry across the genomes of 12 Iowa Sistrurus rattlesnakes. We use whole-genome resequencing data and a maximum likelihood framework to estimate hybrid index in 50 kb windows across all autosomes. We identify ancestry-informative markers (AIMs) through genome-wide FST scans between parental populations (S. catenatus, n=274; S. tergeminus, n=20) and convert high-FST loci to diagnostic markers with fixed allelic differences. We then quantify the proportion of S. tergeminus ancestry per window for each individual, enabling genome-wide visualization of ancestry mosaicism. Our analyses reveal substantial heterogeneity in ancestry both among individuals and across genomic regions, providing insights into the evolutionary dynamics of this ancient hybrid zone. These results have implications for conservation of this threatened species and contribute to broader understanding of hybrid genome evolution.

Materials & Methods
Variant filtering and sample definition
We estimated hybrid ancestry using whole-genome SNP data in VCF format. Analyses were conducted using a custom pipeline (hybrid_index_sliding_window.sh) built around VCFtools and Python scripts. Prior to analysis, the raw VCF was filtered to retain only high-quality biallelic SNPs. Variants were restricted to a predefined set of samples comprising focal individuals (putative hybrids) and two parental populations. SNPs were filtered by minor allele frequency (MAF ≥ 0.05), missing data (maximum missingness ≤ 20%), and indels were removed. Filtering was applied uniformly across all samples to avoid biases in downstream population genetic estimates.
Identification of ancestry-informative markers (AIMs)
To identify ancestry-informative markers, we calculated Weir and Cockerham’s F_STbetween the two parental populations using VCFtools. SNPs with high interspecific differentiation were retained as candidate AIMs using a stringent threshold (F_ST≥0.9). This conservative criterion was chosen to minimize the inclusion of loci with shared polymorphism and to focus on markers that approximate fixed differences between parental taxa.
For each candidate AIM, allele frequencies were estimated separately in each parental population. SNPs were classified as diagnostic when one allele was near fixation in one parental population (frequency ≥ 0.95) and rare in the other (frequency ≤ 0.05). Only these diagnostic loci were retained for hybrid index estimation. This approach corresponds to the “fixed differences” model described by Buerkle (2005).
Extraction of diagnostic genotypes
Genotypes at diagnostic AIMs were extracted from the filtered VCF for all focal (hybrid) individuals. Only these loci were used for downstream ancestry estimation, ensuring that hybrid index inference was based exclusively on markers with strong parental differentiation. 
Maximum-likelihood estimation of hybrid index
Hybrid index was estimated in sliding genomic windows to characterize local ancestry variation along the genome. Windows were defined with a size of 50 kb and a step size of 10 kb. Only windows containing at least 10 diagnostic AIMs were retained for analysis, reducing stochastic variance associated with small numbers of loci.
Within each window and for each individual, hybrid index (h) was estimated using a maximum-likelihood framework following Buerkle (2005). The hybrid index is defined as the proportion of ancestry derived from parental population 2, with h=0corresponding to pure parental population 1 ancestry and h=1corresponding to pure parental population 2 ancestry.
For each diagnostic locus, genotype likelihoods were modeled under a binomial framework assuming Mendelian segregation and no genotyping error. Under the fixed-differences model, the expected allele frequency in a hybrid individual is a linear mixture of parental allele frequencies:
ph=(1−h)p1+hp2,
where p_1and p_2are the allele frequencies in parental populations 1 and 2, respectively. The likelihood of observing each genotype was computed conditional on h, and log-likelihoods were summed across all loci within a window. The maximum-likelihood estimate of hwas obtained by numerical optimization over the interval h∈[0,1].
Approximate 95% confidence intervals for hwere calculated using a likelihood-ratio approach, defined as values of hwithin 2 log-likelihood units of the maximum.
Ancestry classification
For descriptive purposes, windows were classified into ancestry categories based on estimated hybrid index values: windows with h<0.1were classified as parental population 1, windows with h>0.9 as parental population 2, and intermediate values as admixed. These thresholds were used solely for visualization and summary statistics and were not part of the likelihood model.
This approach implements the core maximum-likelihood hybrid index framework described by Buerkle (2005), specifically the diagnostic-marker (fixed-differences) case. It also follows the conceptual foundation of genomic hybrid index estimation discussed by Gompert and Buerkle (2010) but represents a simplified implementation. In particular, the method does not incorporate genotype uncertainty, locus-specific weighting, explicit modeling of linkage disequilibrium, or Bayesian hierarchical structure. These simplifications were adopted to provide a conservative and computationally tractable assessment of local ancestry using strongly differentiated markers.


----------------

General figure description (hybrid index summary plots)

Each panel shows the mean hybrid index (h) across the 11 Iowa samples along a single chromosome/contig (x‑axis = genomic position; y‑axis = hybrid index, scaled 0–1).

iowa_samples

CM078115.hybrid_index_summary

CM078120.hybrid_index_summary

Hybrid index (h) is estimated by maximum likelihood (Buerkle-style hybrid index), where h ranges from 0 to 1 and reflects ancestry along the chromosome (conceptually: 0 = Parent1-like, 1 = Parent2-like).

calculate_hybrid_index_ml

calculate_hybrid_index_ml

Parent populations / colors (background shading):

Blue = Parent1 (S. catenatus)

Tan = Admixed

Pink = Parent2 (S. tergeminus)

hybrid_index_sliding_window

plot_hybrid_index

Black line (“Mean”): the mean h across the 11 Iowa individuals at each window midpoint (windowed across the chromosome).

CM078115.hybrid_index_summary

plot_hybrid_index

Uncertainty bands: gray bands in the legend indicate ±1 SE and ±1 SD around the mean (as displayed in the plot legend).

CM078115.hybrid_index_summary

CM078120.hybrid_index_summary

Horizontal dashed threshold lines: ancestry thresholds (used here as Parent1 threshold = 0.1 and Parent2 threshold = 0.9).

plot_hybrid_index

run_hybrid_script

How the background category is assigned (CI-aware): windows are classified using the window’s confidence interval, i.e. Parent1-like if upper CI ≤ 0.1; Parent2-like if lower CI ≥ 0.9; otherwise admixed.

plot_hybrid_index

Red vertical dotted lines (“Peak” markers): the genomic positions of the top 10 windows with the highest mean h (i.e., the most Parent2-shifted windows in the mean profile), marked with red dotted lines.

plot_hybrid_index

plot_hybrid_index

Yellow callout boxes: nearest annotated genes (typically top 1–2 shown, with “(+ more)” if additional genes fall within the flank distance) for each peak marker.

CM078120.hybrid_index_summary

CM078115.hybrid_index_summary

plot_hybrid_index

(Analysis settings commonly used for these plots: 50 kb windows, 10 kb step, minimum 10 AIMs per window.)

hybrid_index_sliding_window

run_hybrid_script

Plausible explanation for “highly admixed chr1” vs an “introgression desert”
Context: origin of the Iowa lineage

Demographic inference for this system supports a history where the Iowa lineage evolved in allopatry and then experienced a hybridization/introgression event with S. tergeminus (secondary contact), with the timing on the order of ~23 kya (model-based estimate).

hdy201656

hdy201656

hdy201656

Why you can see two very different genomic patterns under that history

Chromosome-wide mosaic admixture (e.g., chr1 example)
Under secondary contact followed by recombination over many generations, introgressed ancestry becomes broken into many tracts distributed across the chromosome. If introgression was not strongly selected against on that chromosome (i.e., relatively few/weak barrier loci), the result can look like a broadly “admixed” chromosome in a sliding-window mean.

A key technical/interpretive point: because you are plotting a mean across individuals, a “highly admixed” profile can also arise when different individuals carry different parental tracts at a given position (between-individual heterogeneity), even if no single individual is uniformly intermediate everywhere. (This is especially likely in older hybrid swarms.)

Introgression desert (e.g., chr6 example)
A long region where the mean hybrid index is consistently near Parent1 (and stays CI-classified as Parent1-like) is consistent with a genomic region that is resistant to introgression. Two biologically plausible mechanisms (not mutually exclusive) are:

Selection against introgressed alleles (Dobzhansky–Muller incompatibilities and/or maladaptation), which purges heterospecific ancestry in that region after secondary contact.

Reduced recombination (centromeric/low-recombination region) and/or structural variants (e.g., inversions) that maintain large linked blocks, making selection act on a broad region and producing an extended “desert.”

Why chr1 might be especially admixed (plausible hypotheses)

Higher effective recombination / fewer barrier loci on chr1 → introgressed tracts are more readily broken up and persist neutrally, yielding widespread mosaic ancestry in the population mean.

Adaptive introgression on chr1 → one or more S. tergeminus alleles (or haplotypes) on chr1 could have been favored in the Iowa environment, increasing retention and spread of Parent2 ancestry across multiple chr1 windows.

Among-individual heterogeneity (mean-of-11 effect) → chr1 might look “admixed” in the mean because different individuals carry different parental segments; this would be testable by checking the corresponding individual-level plots.

Confidence: 0.8 (mechanisms are standard for secondary contact/hybrid zones, but attributing chr1 specifically to recombination vs adaptive introgression vs among-individual heterogeneity requires the individual panels and/or recombination/LD/structural variant evidence).