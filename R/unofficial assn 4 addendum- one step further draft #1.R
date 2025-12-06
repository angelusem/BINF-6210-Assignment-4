#steps adapted from in-class script exercise:
######SOFTWARE TOOLS - CLASS 19 - INTRO TO GENE EXPRESSION ANALYSIS TUTORIAL BY LAW ET AL 2018----

#Authors: This tutorial on RNA-Seq analysis using Bioconductor packages is by: Charity Law, Monther Alhamdoosh, Shian Su, Xueyi Dong, Luyi Tian, Gordon K. Smyth and Matthew E. Ritchie. December 17, 2018. Some commenting added by Sally Adamowicz, last updated November 13, 2023.

#load packages
library(limma)
library(Glimma)
library(edgeR)
library(Mus.musculus)

#If needed, install package R.utils and gplots from CRAN. Load packages.
#install.packages("R.utils")
library(R.utils)
#install.packages("gplots")
library(gplots)

####DATA----

#obtaining data. See online tutorial as well as the source paper (Sheridan et al. 2015).
url <- "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE63310&format=file"
utils::download.file(url, destfile="GSE63310_RAW.tar", mode="wb") 
utils::untar("GSE63310_RAW.tar", exdir = ".")
files <- c("GSM1545535_10_6_5_11.txt", "GSM1545536_9_6_5_11.txt", "GSM1545538_purep53.txt",
           "GSM1545539_JMS8-2.txt", "GSM1545540_JMS8-3.txt", "GSM1545541_JMS8-4.txt",
           "GSM1545542_JMS8-5.txt", "GSM1545544_JMS9-P7c.txt", "GSM1545545_JMS9-P8c.txt")
for(i in paste(files, ".gz", sep="")) {
  R.utils::gunzip(i, overwrite=TRUE)}

#Let's have a look. Let's look at the first 5 rows of file #1.
read.delim(files[1], nrow=5)

#So, this line is combining the 9 files, taking only the first column (the gene IDs) and the third column (the counts) for each file. The data consist of three different types of cells, with three replicates each, hence 9 files. The function readDGE() will combine the 9 samples into a single object for analysis.
x <- readDGE(files, columns=c(1,3))

#what is the class?
class(x)

#Exploring dimensions
dim(x)
length(x)
names(x)

#So, our first element is a dataframe, and the second element is a matrix.
class(x[[1]])
dim(x[[1]])
class(x[[2]])
dim(x[[2]])

#TIP: I always like to have a look at a subset of the data, and I recommend that you do so too. It helps us to understand the contents of a particular data object.

#Let's have a look at the first element of our DGEList. This tells us the information about the samples and the total library size for each sample.
x[[1]]
class(x[[1]])

#This first data frame is called "samples". We can refer to the first element by name:
x$samples

#Now, let's have a look at the first fifty rows of the second element of our DGEList.
x[[2]][1:50, ]

#Checking length of one column of this dataframe against our dimension information above for the DGEList object.
length(x[[2]][, 1])

#The second dataframe is called "counts". So, again we can refer to the second element of our DGEList by name if we wish, e.g.:
x$counts[1:10, ]

#Removing the GEO (Gene Expression Omnibus: https://www.ncbi.nlm.nih.gov/geo/) sample IDs, for simplicity for downstream analysis.

#First, let's look at colnames, so that we can clearly see what the below step is doing.
colnames(x)

#Authors of the tutorial wanted to clean up the sample names by omitting the GEO sample ids. This can be done by starting the sample names at character 12.
samplenames <- substring(text = colnames(x), first = 12, last = nchar(colnames(x)))
samplenames

#EXAMPLE: I am doing this a different way, as an example of using regular expressions to do this. This would be more flexible if there were variable numbers of characters.
testsamplenames <- sub(pattern = "^[A-Za-z]+[0-9]+_", x = colnames(x), replacement = "")

testsamplenames

all.equal(samplenames, testsamplenames)

#setting up group information
colnames(x) <- samplenames
group <- as.factor(c("LP", "ML", "Basal", "Basal", "ML", "LP", "Basal", "ML", "LP"))

#creating a new variable for the groups (i.e. the three types of cells in this analysis)
x$samples$group <- group

#creating a new variable for the lane (on the Illumina sequencer) - it is good to keep information that could have influenced the experiment.
lane <- as.factor(rep(c("L004","L006","L008"), c(3,4,2)))
x$samples$lane <- lane
x$samples

#retrieve annotations for mouse.
geneid <- rownames(x)

#exploring what these gene id's look like: These are Entrez gene ids and can be used to retrieve other information of interest.
head(geneid)

#retrieving gene symbol and chromosome information, using the Entrez gene ids.
genes <- select(x = Mus.musculus, keys = geneid, columns = c("SYMBOL", "TXCHROM"), keytype = "ENTREZID")
head(genes)

genes <- genes[!duplicated(genes$ENTREZID), ]

#adding gene annotations to our main data object
x$genes <- genes

#let's have a look
x

####DATA PREPROCESSING----

#Note that this analysis is looking for DIFFERENCES in gene expression among cell types for each gene, not absolute levels.The below metrics do not consider differences in gene length (which would impact the number of reads mapping to a particular gene).

#converting counts to counts per million and log counts per million
cpm <- cpm(x)
lcpm <- cpm(x, log=TRUE)

#to avoid creating spurious appearance of large fold changes in gene expression for those genes with very low counts. Shrinks inter-sample log-fold changes towards zero.
L <- mean(x$samples$lib.size) * 1e-6
M <- median(x$samples$lib.size) * 1e-6
c(L, M)

summary(lcpm)

#removing genes lowly expressed

#first, summarizing unexpressed genes. 19% of genes in this dataset have a zero count across all 9 samples. (By contrast, we would wish to include genes that are expressed in one condition but not another.)
table(rowSums(x$counts==0)==9)

#filtering unexpressed and lowly expressed genes as well as genes expressed in few samples. The number of reads is considered as well as library size and group size.
keep.exprs <- filterByExpr(x, group = group)
x <- x[keep.exprs,, keep.lib.sizes=FALSE]
dim(x)

#producing a figure to explore the effects of filtering
lcpm.cutoff <- log2(10/M + 2/L)
library(RColorBrewer)
nsamples <- ncol(x)
col <- brewer.pal(nsamples, "Paired")
par(mfrow=c(1,2))
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.26), las=2, main="", xlab="")
title(main="A. Raw data", xlab="Log-cpm")
abline(v=lcpm.cutoff, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", samplenames, text.col=col, bty="n")
lcpm <- cpm(x, log=TRUE)
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.26), las=2, main="", xlab="")
title(main="B. Filtered data", xlab="Log-cpm")
abline(v=lcpm.cutoff, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", samplenames, text.col=col, bty="n")

# segue into assn 4 pipeline:
y_mouse <- x                       #renamING for clarity
group_mouse <- y_mouse$samples$group
group_mouse                        # LP/ML/Basal

#creating normalization variants
y_mouse_TMM    <- calcNormFactors(y_mouse, method = "TMM")
y_mouse_simple <- calcNormFactors(y_mouse, method = "none")
y_mouse_UQ     <- calcNormFactors(y_mouse, method = "upperquartile")

#computing logCPM
logcpm_mouse_TMM    <- cpm(y_mouse_TMM,    log = TRUE)  
logcpm_mouse_simple <- cpm(y_mouse_simple, log = TRUE)
logcpm_mouse_UQ     <- cpm(y_mouse_UQ,     log = TRUE)

#the design matrix – e.g. compare Basal vs LP, Basal vs ML, LP vs ML
design_mouse <- model.matrix(~ 0 + group_mouse)
colnames(design_mouse) <- levels(group_mouse)
design_mouse

#considering the Basal vs LP contrast
contr_mouse <- makeContrasts(
  Basal_vs_LP = Basal - LP,
  levels = design_mouse
)

##DE for each normalization method for the basal vs lp group

#TMM
y_mouse_TMM <- estimateDisp(y_mouse_TMM, design_mouse)
fit_mouse_TMM  <- glmQLFit(y_mouse_TMM, design_mouse)
qlf_mouse_TMM  <- glmQLFTest(fit_mouse_TMM, contrast = contr_mouse[,"Basal_vs_LP"])
tab_mouse_TMM  <- topTags(qlf_mouse_TMM, n = Inf)$table

#Simple
y_mouse_simple <- estimateDisp(y_mouse_simple, design_mouse)
fit_mouse_simple <- glmQLFit(y_mouse_simple, design_mouse)
qlf_mouse_simple <- glmQLFTest(fit_mouse_simple, contrast = contr_mouse[,"Basal_vs_LP"])
tab_mouse_simple <- topTags(qlf_mouse_simple, n = Inf)$table

#UQ
y_mouse_UQ <- estimateDisp(y_mouse_UQ, design_mouse)
fit_mouse_UQ <- glmQLFit(y_mouse_UQ, design_mouse)
qlf_mouse_UQ <- glmQLFTest(fit_mouse_UQ, contrast = contr_mouse[,"Basal_vs_LP"])
tab_mouse_UQ <- topTags(qlf_mouse_UQ, n = Inf)$table

#figure to evaulate whether normalization method varying had an effect in this dataset:
##LOGFC COMPARISON PLOT:(Basal vs LP)

#taking just the common genes
common_genes_mouse <- intersect(rownames(tab_mouse_TMM),
                                rownames(tab_mouse_UQ))

#Building df of logFCs
fc_mouse_compare <- data.frame(
  gene        = common_genes_mouse,
  logFC_TMM   = tab_mouse_TMM[common_genes_mouse, "logFC"],
  logFC_UQ    = tab_mouse_UQ[common_genes_mouse, "logFC"],
  DE_TMM      = tab_mouse_TMM[common_genes_mouse, "FDR"] < 0.05
)

#plotting
library(ggplot2)

p_mouse_fc <- ggplot(fc_mouse_compare,
                     aes(x = logFC_TMM, y = logFC_UQ, colour = DE_TMM)) +
  geom_point(alpha = 0.5, size = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  scale_colour_manual(values = c("FALSE" = "grey70", "TRUE" = "red")) +
  theme_bw() +
  labs(
    title  = "Sheridan mouse data: log2 FC (Basal vs LP)\nTMM vs UQ normalization",
    x      = "log2 fold-change (TMM)",
    y      = "log2 fold-change (UQ)",
    colour = "DE under TMM\n(FDR < 0.05)"
  )

dir.create("figures", showWarnings = FALSE)
ggsave("figures/mouse_logFC_TMM_vs_UQ.png", p_mouse_fc, width = 6, height = 5)


#Jaccard index
sign_mouse_TMM    <- tab_mouse_TMM$FDR    < 0.05
sign_mouse_simple <- tab_mouse_simple$FDR < 0.05
sign_mouse_UQ     <- tab_mouse_UQ$FDR     < 0.05

genes_mouse_TMM    <- rownames(tab_mouse_TMM)[sign_mouse_TMM]
genes_mouse_simple <- rownames(tab_mouse_simple)[sign_mouse_simple]
genes_mouse_UQ     <- rownames(tab_mouse_UQ)[sign_mouse_UQ]

jaccard <- function(a, b) length(intersect(a, b)) / length(union(a, b))

jac_mouse_TMM_simple <- jaccard(genes_mouse_TMM,  genes_mouse_simple)
jac_mouse_TMM_UQ     <- jaccard(genes_mouse_TMM,  genes_mouse_UQ)
jac_mouse_simple_UQ  <- jaccard(genes_mouse_simple, genes_mouse_UQ)

#computing precision & reliability proxies reusing the same helper functions already defined:
#cov_by_group() and mean_within_cor(), just applied to logcpm_mouse_* and group_mouse.

cov_mouse_TMM_LP    <- cov_by_group(logcpm_mouse_TMM,    group_mouse, "LP")
cov_mouse_TMM_Basal <- cov_by_group(logcpm_mouse_TMM,    group_mouse, "Basal")
# reliability metrics for simple/UQ and other groups

rel_mouse_TMM_Basal <- mean_within_cor(logcpm_mouse_TMM, group_mouse, "Basal")
rel_mouse_simple_Basal  <- mean_within_cor(logcpm_mouse_simple, group_mouse, "Basal")
rel_mouse_UQ_Basal      <- mean_within_cor(logcpm_mouse_UQ, group_mouse, "Basal")

