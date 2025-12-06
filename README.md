# BINF-6210-Assignment-4
Differential Gene Expression: Impact of Normalization Choice (airway + edgeR)
Overview

This project investigates how different normalization methods in edgeR affect differential expression (DE) results for RNA-seq data.

Specifically, I compare three normalization strategies applied to the same dataset:

TMM (Trimmed Mean of M-values) – recommended default in edgeR.

Library-size–only scaling (method = "none") – uses raw library sizes as offsets.

Upper-quartile (UQ) normalization – scales based on the 75th percentile of counts.

I assess how these choices influence:

The set of genes called differentially expressed (DE) between treated and untreated samples.

Effect sizes (log2 fold-changes) for individual genes.

Precision (using within-group coefficient of variation, CoV).

Reliability (using within-group Pearson correlation).

The analysis follows the edgeR and limma RNA-seq workflows in Bioconductor 

Data
airways RNA-seq dataset

Source: Bioconductor airway package

Original study: Himes et al., 2014 (human airway smooth muscle cells with/without dexamethasone treatment).

Data type: Bulk RNA-seq counts (genes × samples).

Design:

8 samples total

4 cell lines

Each cell line has control (untrt) and dexamethasone-treated (trt) conditions

Key variables used:

counts: gene-level raw counts

dex: treatment status ("untrt" vs "trt")

Data are loaded directly in R via:

library(airway)
data("airway")


Date accessed: Nov 24, 2025

Research Question

How do different normalization methods in edgeR (TMM, simple library-size scaling, and upper-quartile) influence: the number and overlap of DE genes between treated and untreated samples, the estimated log2 fold-changes, and simple precision/reliability metrics?

This connects to the broader question: Are DE results robust to normalization choices in RNA-seq experiments?

Methods Summary

The analysis is implemented in a single R script:

R/Angelusm_script submission.R

Key steps:

Load data and metadata

Extract counts <- assay(airway) and coldata <- colData(airway)

Define group factor: group <- coldata$dex ("untrt" vs "trt")

Exploratory data analysis / QC

Compute library sizes per sample and plot a barplot.

Create boxplots of log2(counts+1) to inspect raw expression distributions.

Filtering lowly expressed genes

Use filterByExpr() from edgeR to remove genes with insufficient counts.

This reduces noise and focuses the DE analysis on reasonably expressed genes.

Compare three normalization methods

Construct DGEList and apply:

calcNormFactors(y, method = "TMM")

calcNormFactors(y, method = "none")

calcNormFactors(y, method = "upperquartile")

Compute logCPM matrices for downstream comparisons.

Differential expression (edgeR quasi-likelihood pipeline)

Design matrix: design <- model.matrix(~ group)

Estimate dispersion: estimateDisp()

Fit GLM QL model: glmQLFit()

Test treatment effect: glmQLFTest(coef = 2)

Extract full DE tables with topTags(..., n = Inf) for each normalization.

Compare results across normalizations

Count DE genes (FDR < 0.05) for each method.

Compute overlap of DE sets and visualize with a Venn diagram.

Compute Jaccard indices between method pairs and plot as a barplot.

Compare log2 fold-changes across methods with scatter plots:

TMM vs Simple

TMM vs UQ

Precision & reliability metrics

Precision proxy: median coefficient of variation (CoV) of logCPM within each group ("untrt" and "trt").

Reliability proxy: mean within-group Pearson correlation among samples.

MDS plots

Generate MDS plots on normalized logCPM to check whether sample clustering (trt vs untrt) is stable across normalization strategies.

Project Structure (files relevant to submission only)
project_root/
├── R/
│   └── airway_normalization_comparison.R   # main analysis script
├── figures/
│   ├── library_sizes_plot.png
│   ├── ggboxplot_rawED.png
│   ├── ggboxplot_normalized.png
│   ├── jaccard_DE_sets.png
│   ├── VennDiagram_degenesoverlap.png    
│   ├── logFC_TMM_vs_Simple.png
│   ├── logFC_TMM_vs_UQ.png
│   ├── MDS_TMM.png
│   └── MDS_Simple.png
├── doc/
│   └── Angelusm_assn4_storyboard.pdf           
└── README.md                    # this file

How to Run

Install required packages in R:

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("airway", "edgeR", "limma", "SummarizedExperiment"))
install.packages("ggplot2")


Ensure the directory structure exists:

dir.create("figures", showWarnings = FALSE)
dir.create("doc", showWarnings = FALSE)
dir.create("R", showWarnings = FALSE)


Place the main script in R/ and source it from the project root:

source("R/Angelusm_script submission.R")


The script should:

Run end-to-end without errors.

Save all figures into figures/.

Print summary outputs (numbers of DE genes, Jaccard indices, CoV and correlation summaries) to the console.

Key Outputs and brief Interpretation 

Figure: Library sizes
Shows total read counts per sample; here, library sizes show some difference in size, so normalization corrections are warranted.

Figure: Raw vs normalized expression distributions
Boxplots of log2(counts+1) and logCPM show that normalizations yield very similar distributions for this dataset.

Venn + Jaccard plots
DE gene sets (FDR < 0.05) are highly overlapping across TMM, Simple, and UQ normalization, suggesting that in this dataset DE calls are robust to normalization choice.

LogFC comparison plots (TMM vs Simple, TMM vs UQ)
Points lie very close to the y = x line, indicating almost identical log2 fold-change estimates across methods.

MDS plots
Treated vs untreated samples separate clearly under all methods, again indicating limited sensitivity to normalization choice in this relatively well-behaved dataset.

Precision / reliability metrics
Median CoV and mean within-group correlation are very similar across methods, supporting the conclusion that, for airway, changing normalization has only minor impact on precision and reliability.

Future directions

attempt to run script on different datasets: unofficially attempted on in class script example with not much difference-- script available in script folder.

attempt different pipelines on same dataset to see if there is any change

attempt more accurate precision, reliability and accuracy measures, like in referenced papers, with different, more complex data 

References

Chen, Y., Lun, A. T. L., & Smyth, G. K. (2016). From reads to genes to pathways: differential expression analysis of RNA-Seq experiments using Rsubread and the edgeR quasi-likelihood pipeline. F1000Research, 5, 1438. https://doi.org/10.12688/f1000research.8987.2
Chen, Y., Mccarthy, D., Ritchie, M., Robinson, M., Smyth, G., & Hall, E. (2008a). edgeR: differential analysis of sequence read count data User’s Guide. https://www.bioconductor.org/packages/devel/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf
Himes, B. E., Jiang, X., Wagner, P., Hu, R., Wang, Q., Klanderman, B., Whitaker, R. M., Duan, Q., Lasky-Su, J., Nikolos, C., Jester, W., Johnson, M., Panettieri, R. A., Tantisira, K. G., Weiss, S. T., & Lu, Q. (2014). RNA-Seq Transcriptome Profiling Identifies CRISPLD2 as a Glucocorticoid Responsive Gene that Modulates Cytokine Function in Airway Smooth Muscle Cells. PLoS ONE, 9(6), e99625. https://doi.org/10.1371/journal.pone.0099625
Law, C. W., Alhamdoosh, M., Su, S., Dong, X., Tian, L., Smyth, G. K., & Ritchie, M. E. (2018). RNA-seq analysis is easy as 1-2-3 with limma, Glimma and edgeR. F1000Research, 5. https://doi.org/10.12688/f1000research.9005.3
Love, M. I., Anders, S., Kim, V., & Huber, W. (2015). RNA-Seq workflow: gene-level exploratory analysis and differential expression. F1000Research, 4, 1070. https://doi.org/10.12688/f1000research.7035.1
Robinson, M. D., & Oshlack, A. (2010). A scaling normalization method for differential expression analysis of RNA-seq data. Genome Biology, 11(3), R25. https://doi.org/10.1186/gb-2010-11-3-r25
Tong, L., Wu, P.-Y., Phan, J. H., Hassazadeh, H. R., Tong, W., & Wang, M. D. (2020). Impact of RNA-seq data analysis algorithms on gene expression estimation and downstream prediction. Scientific Reports, 10(1), 17925. https://doi.org/10.1038/s41598-020-74567-y

