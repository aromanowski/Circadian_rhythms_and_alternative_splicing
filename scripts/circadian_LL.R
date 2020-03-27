############################################################################
# Circadian_LL.R                                                           #
# --------------                                                           #
# This script analyses circadian gene expression and alternative splicing  #
# It was used to analyse the data described in Romanowski A. et al., 2020  #
#                                                                          #
############################################################################


### Notes
# Many lines of codes are commented. Uncomment them as needed for your usage


###############################################
#             Requires                        #
###############################################
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
  BiocManager::install()}
if(!require(ASpli)) BiocManager::install("ASpli")
if(!require(GenomicFeatures)) BiocManager::install("GenomicFeatures")
if(!require(GenomicRanges)) BiocManager::install("GenomicRanges")
if(!require(AnnotationDbi)) BiocManager::install("AnnotationDbi")
if(!require(BiocParallel)) BiocManager::install("BiocParallel")
if(!require(devtools)) {
   install.packages("devtools")
   devtools::install_github("kassambara/ggpubr")}
if(!require(RColorBrewer)) install.packages("RcolorBrewer")
if(!require(gplots)) install.packages("gplots")


###############################################
#             Includes                        #
###############################################
library(ASpli)
library(BiocParallel)
library(GenomicFeatures)
library(GenomicRanges)
library(AnnotationDbi)
library(gplots)
library(ggpubr)
library(RColorBrewer)

###############################################
#             Begin Analysis                  #
###############################################

# set base directory  ### modify accordingly ###
basedir <- "c:/circadianLL/comp"
setwd(basedir)

# load auxiliary functions
source("../scripts/functions.R")

# Make the file genome sqlite from the transcriptome, using AtRTD2 transcriptome genes.gtf
# This GTF was edited so that it had the correct chromosome naming
# genome<-makeTxDbFromGFF("/data/BioData/references/Arabidopsis_thaliana/AtRTD2/AtRTD2_edited.gtf", format = "gtf") 
# saveDb(genome, file="genomeAtRTD2.sqlite")
# features<-binGenome(genome)
# save(features, file="featuresAtRTD2.sqlite")

# Load the genome
# genome <- loadDb(file="genomeAtRTD2.sqlite")
# load(file="featuresAtRTD2.RData")

###########################################################
# Summary of AtRTD2 features:                             #
#                                                         #
# * Number of extracted Genes = 34212                     #
# * Number of extracted Exon Bins = 238662                #
# * Number of extracted intron bins = 178027              #
# * Number of extracted trascripts = 82190                #
# * Number of extracted junctions = 151944                #
# * Number of AS bins (not include external) = 41863      #
# * Number of AS bins (include external) = 41941          #
# * Classified as:                                        #
#   ES bins = 1686      (4%)                              #
#   IR bins = 13033     (31%)                             #
#   Alt5'ss bins = 4244	(10%)                             #
#   Alt3'ss bins = 7683	(18%)                             #
#   Multiple AS bins = 15217    (36%)                     #
#   classified as:                                        #
#               ES bins = 1627	(11%)                     #
#               IR bins = 5060	(33%)                     #
#               Alt5'ss bins = 2941 (19%)                 #
#               Alt3'ss bins = 5001 (33%)                 #
###########################################################

############################################################
# Create a target file with description of the experiment  #
#                                                          #
# Example                                                  #
# sample  bam condition                                    #
# Col_LL_A  Col_LL_A.bam  ctrl                             #
# Col_LL_B  Col_LL_B.bam  ctrl                             #
# Col_LL_C  Col_LL_C.bam  ctrl                             #
# SMN_LL_A  SMN_LL_A.bam  spf30-1                          #
# SMN_LL_B  SMN_LL_B.bam  spf30-1                          #
# SMN_LL_C  SMN_LL_C.bam  spf30-1                          #
############################################################

# Load the circadian targets and bam files
targets <- read.table("targets.txt", stringsAsFactors = FALSE, header = TRUE, row.names = 1)
# bam <- loadBAM(targets,cores=6) # load BAMS
# save(bam, file="bam.RData") # Save BAM data as an RData file
# load(file="bam.RData") # Load BAM data from RData file

##############################################
# OUTPUT module                              #
# Count tables, PSI, PIR (AS with junctions) #
##############################################
# Get raw counts from the BAM data
# l is read length from the experiment. Here it is set to 100bp
# counts <- readCounts(features, bam, cores=7, readLength = 100, targets=targets, maxISize = 5000)
# save(counts, file="counts.RData") # Save RAW count data
load("counts.RData")

# Get splicing events counts with ASpli
as <- AsDiscover(counts, targets, features, bam, readLength = 100, threshold = 5, cores = 1)
save(as, file="as.RData") # Save AS counts data

#######################################################
# OUTPUT module                                       #
# Differential expression and differential bin usage  #
#######################################################
# Obtain and save Differential Usage of genes and splicing events with ASpli
du <-DUreport(counts = counts, targets = targets, forceGLM = TRUE)
save(du, file="du_GLM.RData")

######################
#    Print results   #
######################
writeAll(counts = counts, du = du, as = as, output.dir = "ASpli_Results") # Generate ASpli reports

#######################################
#         LL data (Raw reads)         #
#######################################
# set base directory
setwd(basedir)

# Load target files
targets <- read.table("targets.txt", stringsAsFactors = FALSE, header = TRUE, row.names = 1) # Raw reads

# Load Counts data
load(file ="counts.RData")

# Get raw counts and define groups for EdgeR
cg <- countsg(counts)
lev = c("24",
        "28",
        "32",
        "36",
        "40",
        "44",
        "48",
        "52",
        "56",
        "60",
        "64",
        "68")
group <- factor(targets$condition, levels = lev)

# Create the design matrix
design <- model.matrix(~group) # Test everything against LL24
colnames(design) <- levels(group)

# Filter and keep genes with rd > 0.05
# Genes that have very low counts across all the libraries should be removed prior to downstream analysis. This is justified on both biological
# and statistical grounds. From biological point of view, a gene must be expressed at some minimal level before it is likely to be translated 
# into a protein or to be considered biologically important. From a statistical point of view, genes with consistently low counts are very unlikely 
# be assessed as significantly DE because low counts do not provide enough statistical evidence for a reliable judgement to be made. Such genes can 
# therefore be removed from the analysis without any loss of information.
df<-filterByRd(cg, targets = targets, min = 0.05, type="any") #18053 genes in total pass the rd filter
save(df, file = "df.RData")

# Create DGEList element with raw counts
# This object is easy to use as it can be manipulated like an ordinary list in R, and it can also be subsetted like a matrix. The main components of
# a DGEList object are a matrix of read counts, sample information in the data.frame format and optional gene annotation. We enter the counts into a
# DGEList object using the function DGEList in edgeR:
y <- DGEList(counts= df[,setdiff(colnames(df), c("symbol", "locus_overlap", "gene_coordinates","start","end","length","effective_length"))] ,
             group=group)
colnames(y)=rownames(targets) # Set Column names
y$samples$lib.size = colSums(y$counts) # Recalculate library sizes

# Plot the library sizes
# The names argument tells the barplot to use the sample names on the x-axis
# The las argument rotates the axis names
par(mar=c(9,4.1,4.1,2.1))
barplot(y$samples$lib.size, names=colnames(y), las=2)
# Add a title to the plot
title("Barplot of library sizes")

# Save the library sizes plot
png(file="library_sizes.png",    # create PNG for the MDS        
    width = 8*600,        # 10 x 600 pixels
    height = 5*600,
    res = 600,            # 300 pixels per inch
    pointsize = 8)        # smaller font size))
par(mar=c(9,4.1,4.1,2.1))
barplot(y$samples$lib.size, names=colnames(y), las=2)
# Add a title to the plot
title("Barplot of library sizes")
dev.off()


# Count data is not normally distributed, so if we want to examine the distributions of the raw counts we need to log the counts. 
# We will use box plots to check the distribution of the read counts on the log2 scale. 
# We can use the cpm function to get log2 counts per million, which are corrected for the different library sizes. 
# The cpm function also adds a small offset to avoid taking log of zero.
logcounts <- cpm(y,log=TRUE)
save(logcounts, file = "logcounts.RData")

# Check distributions of samples using boxplots
boxplot(logcounts, xlab="", ylab="Log2 counts per million",las=2)
# Let's add a blue horizontal line that corresponds to the median logCPM
abline(h=median(logcounts),col="blue")
title("Boxplots of logCPMs (unnormalised)")

png(filename = (paste0("logCPM_Boxplot",".png")),    # create PNG for the MDS        
    width = 10*600,        # 10 x 600 pixels
    height = 5*600,
    res = 600,            # 300 pixels per inch
    pointsize = 8)
par(mar=c(9,4.1,4.1,2.1))
boxplot(logcounts, xlab="", ylab="Log2 counts per million",las=2)
# Let's add a blue horizontal line that corresponds to the median logCPM
abline(h=median(logcounts),col="blue")
title("Boxplots of logCPMs (unnormalised)")
dev.off()

# Correlation analysis of replicate samples
library(ggpubr)
dir.create("correlation")
df_logcounts <- data.frame(logcounts)
for (i in (0:((length(colnames(df_logcounts))/2)-1)))
{
  print(2*i+1)
  print(colnames(df_logcounts)[2*i+1])
  print(2*i+1+1)
  print(colnames(df_logcounts)[2*i+1+1])
  print(" - ")
  gg <- ggscatter(data.frame(df_logcounts[,(2*i+1):(2*i+2)]),
                  x = colnames(df_logcounts)[2*i+1],
                  y = colnames(df_logcounts)[2*i+2],
                  color = "red",
                  add = "reg.line",
                  add.params = list(color = "black", fill = "lightgray"),
                  conf.int = TRUE,
                  cor.coef = TRUE,
                  cor.coeff.args = list(method = "pearson", label.sep = "\n"),
                  xlab = colnames(df_logcounts)[2*i+1],
                  ylab = colnames(df_logcounts)[2*i+2]
  )
  ggarrange(gg) %>% 
    ggexport(filename = paste0("correlation/", i, " - log2 ", colnames(df_logcounts)[2*i+1], " vs ", colnames(df_logcounts)[2*i+2], ".png"),
             pointsize = 12,
             height = 4.5 * 600,
             width = 4.5 * 600,
             res = 600)
}
dev.off()

# Create a Multi Dimensional Scaling Plot
# The RNA samples can be clustered in two dimensions using multi-dimensional scaling (MDS) plots. This is both an analysis step and a quality
# control step to explore the overall differences between the expression profiles of the different samples.
# In the MDS plot, the distance between each pair of samples can be interpreted as the leading log-fold change between the samples for the genes
# that best distinguish that pair of samples. By default, leading fold-change is defined as the root-mean-square of the largest 500 log2-fold changes
# between that pair of samples.

# We specify the option to let us plot only one plot
par(mfrow=c(1,1))
# Check number and order of the samples
levels(factor(targets$condition))
## Choose different colours for each sample type
col.cell <- c("violet", "purple", "blue", "cyan", "green", "dark gray","orange", "pink", "magenta", "red", "dark red", "brown")[factor(targets$condition)]
data.frame(targets$condition,col.cell)

# MDS with sample type colouring
plotMDS(y,col=col.cell)
dev.off()

png(filename = (paste0("MDS_plot",".png")),    # create PNG for the MDS        
    width = 10*600,        # 10 x 600 pixels
    height = 5*600,
    res = 600,            # 300 pixels per inch
    pointsize = 8)        # smaller font size)
plotMDS(y,col=col.cell)
title("MDS by Sample")
dev.off()

# Hierarchical clustering with heatmaps
# We estimate the variance for each row in the logcounts matrix
var_genes <- apply(logcounts, 1, var)
head(var_genes)
# Get the gene names for the top 500 most variable genes
select_var <- names(sort(var_genes, decreasing=TRUE))[1:500]
head(select_var)
# Subset logcounts matrix
highly_variable_lcpm <- logcounts[select_var,]
dim(highly_variable_lcpm)
head(highly_variable_lcpm)

## Get some nicer colours
library(RColorBrewer)
mypalette <- brewer.pal(7,"Spectral")
morecols <- colorRampPalette(mypalette)
# Set up colour vector for sample variable
col.cell <- c("violet", "purple", "blue", "cyan", "green", "dark gray","orange", "pink", "magenta", "red", "dark red", "brown")[factor(targets$condition)]

# Plot the heatmap
library(gplots)
par(mar=c(9,4.1,4.1,2.1))
heatmap.2(highly_variable_lcpm, 
          col=rev(morecols(50)),
          trace="column", 
          main="Variability of genes across samples",
          ColSideColors=col.cell,
          scale="row",
          cexRow=1,
          cexCol=1,
          margins=c(7,6),
          srtCol=45,
          labRow=FALSE)

# Save the heatmap
png(file="Top500_var_genes.heatmap.png",    # create PNG for the MDS        
    width = 10*600,        # 10 x 600 pixels
    height = 5*600,
    res = 600,            # 300 pixels per inch
    pointsize = 8)        # smaller font size))
par(mar=c(9,4.1,4.1,2.1))
heatmap.2(highly_variable_lcpm, 
          col=rev(morecols(50)),
          trace="column", 
          main="Variability of genes across samples",
          ColSideColors=col.cell,scale="row",
          cexRow=1,
          cexCol=1,
          margins=c(7,6),
          srtCol=45,
          labRow=FALSE)
dev.off()

# Hierarchical clustering with heatmaps
# We estimate the variance for each row in the logcounts matrix
var_genes <- apply(logcounts, 1, var)
head(var_genes)
# Get the gene names for all genes
select_var <- names(sort(var_genes, decreasing=TRUE))
head(select_var)
# Subset logcounts matrix
highly_variable_lcpm <- logcounts[select_var,]
dim(highly_variable_lcpm)
head(highly_variable_lcpm)

## Get some nicer colours
library(RColorBrewer)
mypalette <- brewer.pal(7,"Spectral")
morecols <- colorRampPalette(mypalette)
# Set up colour vector for sample variable
col.cell <- c("violet", "purple", "blue", "cyan", "green", "dark gray","orange", "pink", "magenta", "red", "dark red", "brown")[factor(targets$condition)]

# Plot the heatmap
library(gplots)
par(mar=c(9,4.1,4.1,2.1))
heatmap.2(highly_variable_lcpm, 
          col=rev(morecols(50)),
          trace="column", 
          main="Variability of genes across samples",
          ColSideColors=col.cell,scale="row",
          cexRow=1,
          cexCol=1,
          margins=c(7,6),
          srtCol=45,
          labRow=FALSE)
dev.off()

# Save the heatmap
png(file="Full_var_genes.heatmap.png",    # create PNG for the MDS        
    width = 10*600,        # 10 x 600 pixels
    height = 5*600,
    res = 600,            # 300 pixels per inch
    pointsize = 8)        # smaller font size))
par(mar=c(9,4.1,4.1,2.1))
heatmap.2(highly_variable_lcpm, 
          col=rev(morecols(50)),
          trace="column", 
          main="Variability of genes across samples",
          ColSideColors=col.cell,scale="row",
          cexRow=1,
          cexCol=1,
          margins=c(7,6),
          srtCol=45,
          labRow=FALSE)
dev.off()

########################################################################
#         Leaf Dev (Fitted values and normalized expression values     #
########################################################################

# Create design model to define comparison groups
# Linear modeling and differential expression analysis in edgeR requires a design matrix to be specified. The design matrix records which
# treatment conditions were applied to each samples, and it also defines how the experimental effects are parametrized in the linear models.
#
# Similar example from Limma (see user guide), with a timecourse:
# lev <- c("wt.0hr","wt.6hr","wt.24hr","mu.0hr","mu.6hr","mu.24hr")
# f <- factor(targets$Target, levels=lev)
# design <- model.matrix(~0+f)
# colnames(design) <- lev
# fit <- lmFit(eset, design)

lev = c("24",
        "28",
        "32",
        "36",
        "40",
        "44",
        "48",
        "52",
        "56",
        "60",
        "64",
        "68")
group <- factor(targets$condition, levels = lev)
design <- model.matrix(~group, data = y$samples) # Test everythings against LL24
colnames(design) <- unique(targets$condition)
save(design, file = "design.RData")

# Calculate normalization factors by library
# Normalization by trimmed mean of M values (TMM) (Robinson and Oshlack 2010) is performed by using the calcNormFactors function, which returns
# the DGEList argument with only the  norm.factors changed. It calculates a set of normalization factors, one for each sample, to eliminate 
# composition biases between libraries. The product of these factors and the library sizes defines the effective library size, which replaces 
# the original library size in all downstream analyses.
y<-calcNormFactors(y)

# Check before and after TMM normalization effect
# The expression profiles of individual samples can be explored more closely with mean-difference (MD) plots. An MD plot visualizes the library
# size-adjusted log-fold change between two libraries (the difference) against the average log-expression across those libraries (the mean).
par(mfrow=c(2,2)) # plot 2 by 2 graphs in the same layout
plotMD(logcounts,column = 7) # before normalization
abline(h=0,col="blue")
plotMD(logcounts,column = 8) # before normalization
abline(h=0,col="blue")
plotMD(y,column = 7) # after TMM normalization
abline(h=0,col="red")
plotMD(y,column = 8) # after TMM normalization
abline(h=0,col="red")

# Normalization graphs (one per sample side by side)
pdf("circadian_samples_before_and_after_TMM_normalization.pdf",
    paper = "a4r")
for (i in seq(from = 0, to = length(rownames(targets)), by=6)) {
  par(mfrow=c(2,4), mar=c(9,4.1,4.1,2.1))
  
  print(i)
  for (j in seq(c(0:7))) {
    if ((i+j) <= length(rownames(targets))) {  
      print(paste0("j =", j, " and i=", i+j))
      plotMD(logcounts,column = i+j) # before normalization
      abline(h=0,col="blue")
      plotMD(y,column = i+j) # after TMM normalization
      abline(h=0,col="red")
    }
  }
  if (i > length(rownames(targets))) { i = length(rownames(targets))}
  
}
dev.off()

# Estimate dispersions:
# EdgeR estimates an empirical Bayes moderated dispersion for each individual gene. It also estimates a common dispersion, which is a
# global dispersion estimate averaged over all genes, and a trended dispersion where the dispersion of a gene is predicted from its abundance
y<-estimateDisp(y, design = design, robust = TRUE)
save(y, file = "y.RData")
# load(file = "y.RData")

# This returns a DGEList object with additional components (common.dispersion,  trended.dispersion and tagwise.dispersion) added to hold the 
# estimated dispersions. Here robust=TRUE has been used to protect the empirical Bayes estimates against the possibility of outlier genes with 
# exceptionally large or small individual dispersions (Phipson et al. 2016).

# Plot the dispersions:
# The vertical axis of the plotBCV plot shows square-root dispersion, also known as biological coefficient of variation (BCV) 
# (McCarthy, Chen, and Smyth 2012). For RNA-seq studies, the NB dispersions tend to be higher for genes with very low counts. 
# The dispersion trend tends to decrease smoothly with abundance and to asymptotic to a constant value for genes with larger counts.
par(mfrow=c(1,1)) # Plot only one graph
plotBCV(y) # Plot the dispersion of the data

# Save the dispersion plot
png(file="circadian_data_dispersions.png",    # create PNG for the MDS        
    width = 9*600,        # 10 x 600 pixels
    height = 5*600,
    res = 600,            # 300 pixels per inch
    pointsize = 8)        # smaller font size))
plotBCV(y) # Plot the dispersion of the data
dev.off()

# Calculate fitted.values counts (pseudo counts) and perform ANOVA like function (Quasi Likelihood F-test)

# The NB model can be extended with quasi-likelihood (QL) methods to account for gene-specific variability from both biological and technical
# sources (Lund et al. 2012; Lun, Chen, and Smyth 2016). Under the QL framework, the NB dispersion trend is used to describe the overall biological 
# variability across all genes, and gene-specific variability above and below the overall level is picked up by the QL dispersion. In the QL approach,
# the individual (tagwise) NB dispersions are not used.
# The estimation of QL dispersions is performed using the glmQLFit function.
# Setting robust=TRUE in glmQLFit is usually recommended (Phipson et al. 2016). This allows gene-specific prior df estimates, with lower values for
# outlier genes and higher values for the main body of genes. This reduces the chance of getting false positives from genes with extremely high or low
# raw dispersions, while at the same time increasing statistical power to detect differential expression for the main body of genes.
fit <- glmQLFit(y, design = design, robust = TRUE)
save(fit, file = "fit.RData")

# This returns a DGEGLM object with the estimated values of the GLM coefficients for each gene. It also contains a number of empirical Bayes (EB)
# statistics including the QL dispersion trend, the squeezed QL dispersion estimates and the prior degrees of freedom (df). The QL dispersions can be
# visualized by plotQLDisp function.
plotQLDisp(fit)


# Save the QL dispersion plot
png(file="QL_dispersions.png",    # create PNG for the MDS        
    width = 9*600,        # 10 x 600 pixels
    height = 5*600,
    res = 600,            # 300 pixels per inch
    pointsize = 8)        # smaller font size))
plotQLDisp(fit) # Plot the QL dispersion of the data
dev.off()

# We use QL F-tests instead of the more usual likelihood ratio tests (LRT) as they give stricter error rate control by accounting for the uncertainty
# in dispersion estimation
# See: https://support.bioconductor.org/p/84338/
qlf <- glmQLFTest(fit, coef = c(2:12))

fdr<-p.adjust(qlf$table$PValue, method="BH") # Calculate FDR
tt<-cbind(qlf$table, fdr) 
final<-list(qlf, tt)
names(final)=c("qlf", "full")


# Genes that pass read density filter and likelihood test
dim(final$qlf) # 18503 genes
genes_qlf <- dim(final$qlf)[1]
index_qlf <- which(final$full$PValue < 0.05 & final$full$fdr < 0.1) # 13256 / 18503 genes
save(index_qlf, file = "index_qlf.RData")
genes <- final
genesfv <- genes$qlf$fitted.values

# Calculate normalized logCPM data
logcounts2 <- cpm(y, normalized.lib.sizes = TRUE, log=TRUE)
save(logcounts2, file = "logcounts2.RData")

# Calculate normalized counts (cpm)
norm_counts <- cpm(y, normalized.lib.sizes = TRUE) 
save(norm_counts, file = "norm_counts.RData")

# Plot an example gene
par(mfrow=c(2,2), mar=c(9,4.1,4.1,2.1))
barplot(norm_counts["AT4G16780", ], las = 2, ylab = "Normalised counts (CPM)") # Example value
title("HXK1 - AT4G16780")
barplot(norm_counts["AT2G46830", ], las = 2, ylab = "Normalised counts (CPM)") # Example value
title("CCA1 - AT2G46830")
barplot(norm_counts["AT1G02340", ], las = 2, ylab = "Normalised counts (CPM)") # Example value
title("HFR1 - AT1G02340")
barplot(norm_counts["AT5G61380", ], las = 2, ylab = "Normalised counts (CPM)") # Example value
title("TOC1 - AT5G61380")

# Plot an example gene
png(file="gene_controls_clock_and_shade_responders.png",    # create PNG for the MDS        
    width = 9*600,        # 10 x 600 pixels
    height = 5*600,
    res = 600,            # 300 pixels per inch
    pointsize = 8)        # smaller font size))
par(mfrow=c(2,2), mar=c(9,4.1,4.1,2.1))
barplot(norm_counts["AT4G16780", ], las = 2, ylab = "Normalised counts (CPM)") # Example value
title("HXK1 - AT4G16780")
barplot(norm_counts["AT2G46830", ], las = 2, ylab = "Normalised counts (CPM)") # Example value
title("CCA1 - AT2G46830")
barplot(norm_counts["AT1G02340", ], las = 2, ylab = "Normalised counts (CPM)") # Example value
title("HFR1 - AT1G02340")
barplot(norm_counts["AT5G61380", ], las = 2, ylab = "Normalised counts (CPM)") # Example value
title("TOC1 - AT5G61380")
dev.off()

par(mfrow = c(1,1))
# Check distributions of samples using boxplots
boxplot(logcounts2, xlab="", ylab="Log2 counts per million",las=2)
# Let's add a blue horizontal line that corresponds to the median logCPM
abline(h=median(logcounts2),col="blue")
title("Boxplots of logCPMs (normalised)")

png(filename = (paste0("logCPM_normalized_Boxplot",".png")),    # create PNG for the MDS        
    width = 10*600,        # 10 x 600 pixels
    height = 5*600,
    res = 600,            # 300 pixels per inch
    pointsize = 8)
par(mar=c(9,4.1,4.1,2.1))
boxplot(logcounts2, xlab="", ylab="Log2 counts per million",las=2)
# Add a blue horizontal line that corresponds to the median logCPM
abline(h=median(logcounts2),col="blue")
title("Boxplots of logCPMs (normalised)")
dev.off()

# Correlation analysis of replicate samples
library(ggpubr)
dir.create("correlation")
df_logcounts2 <- data.frame(logcounts2)
for (i in (0:((length(colnames(df_logcounts2))/2)-1)))
{
  print(2*i+1)
  print(colnames(df_logcounts2)[2*i+1])
  print(2*i+1+1)
  print(colnames(df_logcounts2)[2*i+1+1])
  print(" - ")
  gg <- ggscatter(data.frame(df_logcounts2[,(2*i+1):(2*i+2)]),
                  x = colnames(df_logcounts2)[2*i+1],
                  y = colnames(df_logcounts2)[2*i+2],
                  color = "red",
                  add = "reg.line",
                  add.params = list(color = "black", fill = "lightgray"),
                  conf.int = TRUE,
                  cor.coef = TRUE,
                  cor.coeff.args = list(method = "pearson", label.sep = "\n"),
                  xlab = colnames(df_logcounts2)[2*i+1],
                  ylab = colnames(df_logcounts2)[2*i+2]
  )
  ggarrange(gg) %>% 
    ggexport(filename = paste0("correlation/", i, " - normalised log2 ", colnames(df_logcounts2)[2*i+1], " vs ", colnames(df_logcounts2)[2*i+2], ".png"),
             pointsize = 12,
             height = 4.5 * 600,
             width = 4.5 * 600,
             res = 600)
}
dev.off()


# Create a Multi Dimensional Scaling Plot of the normalized data
## Let's choose colours for the different samples
col.cell <- c("violet", "purple", "blue", "cyan", "green", "dark gray","orange", "pink", "magenta", "red", "dark red", "brown")[factor(targets$condition)]
data.frame(targets$condition,col.cell)
plotMDS(y, col=col.cell) # Plot MDS

# Save the MDS plot to a png file
png(filename = (paste0("MDS_plot_normalized_data",".png")),    # create PNG for the MDS        
    width = 10*600,        # 10 x 600 pixels
    height = 5*600,
    res = 600,            # 300 pixels per inch
    pointsize = 8)        # smaller font size)
# We specify the option to let us plot only one plot
par(mfrow=c(1,1))
plotMDS(y,col=col.cell)
title("MDS by Sample (normalized data)")
dev.off()

# Hierarchical clustering with heatmaps
# We estimate the variance for each row in the logcounts matrix (normalized values)
var_genes <- apply(logcounts2[index_qlf,], 1, var)
head(var_genes)
# Get the gene names for the top 500 most variable genes
select_var <- names(sort(var_genes, decreasing=TRUE))[1:500]
head(select_var)
# Subset logcounts2 matrix
highly_variable_lcpm <- logcounts2[index_qlf,][select_var,]
dim(highly_variable_lcpm)
head(highly_variable_lcpm)

## Get some nicer colours
library(RColorBrewer)
mypalette <- brewer.pal(7,"Spectral")
morecols <- colorRampPalette(mypalette)
# Set up colour vector for sample variable
col.cell <- c("violet", "purple", "blue", "cyan", "green", "dark gray","orange", "pink", "magenta", "red", "dark red", "brown")[factor(targets$condition)]

# Plot the heatmap
library(gplots)
par(mar=c(9,4.1,4.1,2.1))
heatmap.2(highly_variable_lcpm, 
          col=rev(morecols(50)),
          trace="column", 
          main="Variability of genes across samples (normalized data)",
          ColSideColors=col.cell,scale="row",
          cexRow=1,
          cexCol=1,
          margins=c(7,6),
          srtCol=45,
          labRow=FALSE)
dev.off()

# Save the heatmap
png(file="Top500_var_genes_normalized_data.heatmap.png",    # create PNG for the MDS        
    width = 10*600,        # 10 x 600 pixels
    height = 5*600,
    res = 600,            # 300 pixels per inch
    pointsize = 8)        # smaller font size))
par(mar=c(9,4.1,4.1,2.1))
heatmap.2(highly_variable_lcpm, 
          col=rev(morecols(50)),
          trace="column", 
          main="Variability of genes across samples (normalized data)",
          ColSideColors=col.cell,scale="row",
          cexRow=1,
          cexCol=1,
          margins=c(7,6),
          srtCol=45,
          labRow=FALSE)
dev.off()

# We estimate the variance for each row in the logcounts matrix (normalized values)
var_genes <- apply(logcounts2[index_qlf,], 1, var)
head(var_genes)
# Get the gene names for the all variable genes
select_var <- names(sort(var_genes, decreasing=TRUE))
head(select_var)
# Subset logcounts2 matrix
highly_variable_lcpm <- logcounts2[index_qlf,][select_var,]
dim(highly_variable_lcpm)
head(highly_variable_lcpm)

## Get some nicer colours
library(RColorBrewer)
mypalette <- brewer.pal(7,"Spectral")
morecols <- colorRampPalette(mypalette)
# Set up colour vector for sample variable
col.cell <- c("violet", "purple", "blue", "cyan", "green", "dark gray","orange", "pink", "magenta", "red", "dark red", "brown")[factor(targets$condition)]

# Plot the heatmap
library(gplots)
par(mar=c(9,4.1,4.1,2.1))
heatmap.2(highly_variable_lcpm, 
          col=rev(morecols(50)),
          trace="column", 
          main="Variability of genes across samples (normalized data)",
          ColSideColors=col.cell,scale="row",
          cexRow=1,
          cexCol=1,
          margins=c(7,6),
          srtCol=45,
          labRow=FALSE)
dev.off()

# Save the heatmap
png(file="Full_var_genes_normalized_data.heatmap.png",    # create PNG for the MDS        
    width = 10*600,        # 10 x 600 pixels
    height = 5*600,
    res = 600,            # 300 pixels per inch
    pointsize = 8)        # smaller font size))
par(mar=c(9,4.1,4.1,2.1))
heatmap.2(highly_variable_lcpm, 
          col=rev(morecols(50)),
          trace="column", 
          main="Variability of genes across samples (normalized data)",
          ColSideColors=col.cell,scale="row",
          cexRow=1,
          cexCol=1,
          margins=c(7,6),
          srtCol=45,
          labRow=FALSE)
dev.off()


# Save all count data
write.table(qlf$fitted.values[index_qlf,], file="circadian.genes.fitted.values.tab", sep="\t", row.names = TRUE, col.names = TRUE, na = "NA")
write.table(cg, file="circadian.genes.raw.values.tab", sep="\t", row.names = TRUE, col.names = TRUE, na = "NA")
write.table(df, file="circadian.genes.filtered.values.tab", sep="\t", row.names = TRUE, col.names = TRUE, na = "NA")
write.table(logcounts[index_qlf,], file="circadian.genes.logCPM.unnormalised.values.tab", sep="\t", row.names = TRUE, col.names = TRUE, na = "NA")
write.table(logcounts2[index_qlf,], file="circadian.genes.logCPM.normalised.values.tab", sep="\t", row.names = TRUE, col.names = TRUE, na = "NA")
write.table(norm_counts[index_qlf,], file="circadian.genes.CPM.normalised.values.tab", sep="\t", row.names = TRUE, col.names = TRUE, na = "NA")
write.table(logcounts, file="full.circadian.genes.logCPM.unnormalised.values.tab", sep="\t", row.names = TRUE, col.names = TRUE, na = "NA")
write.table(logcounts2, file="full.circadian.genes.logCPM.normalised.values.tab", sep="\t", row.names = TRUE, col.names = TRUE, na = "NA")
write.table(norm_counts, file="full.circadian.genes.CPM.normalised.values.tab", sep="\t", row.names = TRUE, col.names = TRUE, na = "NA")

# Add annotation data to normalised counts and logCPM data
md <- read.csv("gene.metadata.tab", na.strings=c("", "NA"), sep = "\t", header = TRUE, stringsAsFactors = FALSE)
logcounts2.annot <- getAnnnot(logcounts2, md)
norm_counts.annot <- getAnnnot(norm_counts, md)
write.table(logcounts2.annot[index_qlf,], file="annot.circadian.genes.logCPM.normalised.values.tab", sep="\t", row.names = TRUE, col.names = TRUE, na = "NA")
write.table(norm_counts.annot[index_qlf,], file="annot.circadian.genes.CPM.normalised.values.tab", sep="\t", row.names = TRUE, col.names = TRUE, na = "NA")

############################################
#       Circadian Analysis - JTK Cycle     #
############################################

# Set basedir
if (get_os() != "windows") {basedir <- "/data/projects/circadianLL/comp" 
} else {basedir <- "c:/circadianLL/comp"}
setwd(basedir)

# Create the output directory within the project folder
dir.create("JTK")
  
# Load JTK Cycle functions
source("JTKversion3/JTK_CYCLEv3.1.R")

#### Gene analysis ####
project <- "genes"
options(stringsAsFactors=FALSE)
annot <- norm_counts.annot[index_qlf,]$Symbol
data <- norm_counts[index_qlf,] # genes that passed the rd filter, and the GLM + QL F-test (p<0.01 and q<0.01)
jtkdist(ncol(data))


jtkdist(12,2) # 12 total time points, 2 replicates per time point
periods <- 5:7  # looking for rhythms between 20-28 hours (i.e. between 5 and 7 time points per cycle).
jtk.init(periods, 4) # sampling rate in hours

cat("JTK analysis started on",date(),"\n")
flush.console()

st <- system.time({
  res <- apply(data,1,function(z) {
    jtkx(z)
    c(JTK.ADJP,JTK.PERIOD,JTK.LAG,JTK.AMP)
  })
  res <- as.data.frame(t(res))
  bhq <- p.adjust(unlist(res[,1]),"BH")
  res <- cbind(bhq,res)
  colnames(res) <- c("BH.Q","ADJ.P","PER","LAG","AMP")
  results <- cbind(annot,res,data)
  results <- results[order(res$ADJ.P,-res$AMP),]
})
print(st)

gene.coord.md <- read.table("gene.metadata.tab")
save(results,file=paste("JTK/JTK",project,"rda",sep="."))
write.table(results,file=paste("JTK/JTK",project,"txt",sep="."),row.names=T,col.names=T,quote=F,sep="\t")
write.table(getAnnnot(results,gene.coord.md),file=paste("JTK/annot.JTK",project,"txt",sep="."),row.names=T,col.names=T,quote=F,sep="\t")

# recalculate phases to be circadian  
results2 <-data.frame(cbind(results[,1:5],
                            recalculated.phase = results$LAG/(results$PER/24),
                            results[,6:30]))
save(results2, file = paste("JTK/JTK_recalculated_",project,"rda",sep="."))
write.table(results,file=paste("JTK/JTK_recalculated_",project,"txt",sep="."),row.names=T,col.names=T,quote=F,sep="\t")
write.table(getAnnnot(results2,gene.coord.md),file=paste("JTK/annot.JTK_recalculated",project,"txt",sep="."),row.names=T,col.names=T,quote=F,sep="\t")

################################
# Plot Circadian gene Heatmaps #
################################
setwd(basedir)
dir.create("heatmaps")

load(file = "JTK/JTK.genes.rda")
load(file = "JTK/JTK_recalculated_.genes.rda")

length(which(results$ADJ.P<0.001 & results$BH.Q < 0.001)) #6944
length(which(results$ADJ.P<0.001 & results$BH.Q < 0.01)) #7418    
length(which(results$ADJ.P<0.01 & results$BH.Q < 0.01)) #9127
length(which(results$ADJ.P<0.05 & results$BH.Q < 0.01)) #9127
length(which(results$ADJ.P<0.05 & results$BH.Q < 0.05)) #10608

results <- results[which(results$ADJ.P<0.01 & results$BH.Q < 0.01),] # 9127 genes with p < 0.01 and q < 0.01
results2 <- results2[which(results2$ADJ.P<0.01 & results2$BH.Q < 0.01),] # 9127 genes with p < 0.01 and q < 0.01

results <- results[order(results$LAG, results$AMP),] # Tau phases
results2 <- results2[order(results2$recalculated.phase, results$AMP),] # Circadian phases
indexGenes <- match(rownames(results), rownames(logcounts2))

# Create a scaled logcounts to use with the plot.clust function
scaled_logcounts2 <- t(scale(t(as.matrix(logcounts2[indexGenes,]))))

# Create a cluster type heatmap of gene expression
# create the auxiliary data needed
cluster <- unique(results$LAG)
j <- NULL
for (i in seq_along(cluster)){
  print(i)
  j[i] <- length(which(results$LAG == cluster[i]))
} 
size <- j

# Plot the cluster
plot.clust(list(data = logcounts2[indexGenes,],
                       cluster = results$LAG,
                       size = size),
                standardize.gene = TRUE,
                gene_labels = FALSE)
dev.off()

# Save the cluster heatmap
png(file="heatmaps/circadian.genes.heatmap.png",    # create PNG for the MDS        
    width = 5*1200,        # 6.5 x 1200 pixels  
    height = 5*1200,
    res = 1200,            # 1200 pixels per inch
    pointsize = 10)        # smaller font size
plot.clust(list(data = logcounts2[indexGenes,],
                      cluster = results$LAG,
                      size = size),
                 standardize.gene = TRUE,
                 gene_labels = FALSE)
dev.off()

# Create a scale for the plot
png(file="heatmaps/circadian.genes.heatmap.scale.png",    # create PNG for the MDS        
    width = 6*600,        # 6.5 x 1200 pixels
    height = 2*600,
    res = 600,            # 1200 pixels per inch
    pointsize = 10)        # smaller font size
image.scale(z = scaled_logcounts2, zlim = c(min(scaled_logcounts2), max(scaled_logcounts2)), col = colorRampPalette(c("blue", "gray","yellow"))(64))
dev.off()

###############################################
#       Splicing events analysis              #
###############################################
# set base directory
if (get_os() != "windows") {basedir <- "/data/projects/circadianLL/comp" 
} else {basedir <- "c:/circadianLL/comp"}
setwd(basedir)


load(file = "counts.RData")

lev = c("24",
        "28",
        "32",
        "36",
        "40",
        "44",
        "48",
        "52",
        "56",
        "60",
        "64",
        "68")
group <- factor(targets$condition, levels = lev)
coef <- c(2:12)

# Create the design matrix
design <- model.matrix(~group) # Contrast against LL24
colnames(design) <- levels(group)

cg <- countsg(counts)
exon.intron.counts <- countsb(counts)

norm = TRUE
ignoreExtremes = TRUE

dfGen0=cg
dfBin0=exon.intron.counts #374748
dfGen<-filterByRd(dfGen0, targets, min = 0.05, type="all") # 15529 genes
dfBin=dfBin0[dfBin0[,"locus"]%in%row.names(dfGen),] # 269290 bins

if (ignoreExtremes==TRUE) {
  dfBin=dfBin[dfBin$event!="external",] 
}
# 205091 bins (not considering external bins)

########### We need to pass the filter rd bin / rd gen > 5% ###############
df<-filterByRdBin(dfGen=dfGen,
                  dfBin=dfBin,
                  targets=targets, 
                  min=0.05) 
# 177619 bins have an rd bin / rd gen > 0.05

  y <- DGEList(counts = df[,setdiff(colnames(df), 
                               c("feature","event", "locus", "locus_overlap", "symbol", "gene_coordinates", "start", "end", "length"))] ,
               group = group)
colnames(y)=rownames(targets)
y$samples$lib.size = colSums(y$counts)
design=model.matrix(~group, data=y$samples) #test against LL24
colnames(design)<-levels(y$samples$group)
y<-calcNormFactors(y)

# MDS with sample type colouring
plotMDS(y,col=col.cell)
  
y<-estimateDisp(y, design = design, robust = TRUE)
fit <- glmQLFit(y, design = design, robust = TRUE)

###########################################################################
if  (norm==TRUE) {
  start=ncol(df) - nrow(targets) +1
  end=ncol(df)
  
  bin2gen<-match(df$locus, row.names(genesfv))
  nmatrix<-genesfv[bin2gen,]
  gene.mean=rowMeans(nmatrix)
  counts.N<-round(fit$fitted.values  /  nmatrix * gene.mean)
  counts.N[is.na(counts.N)] <- 0 # we have to remove NAs
  counts.final.N<-cbind(df[,1:start -1], counts.N) # have a df equal to original df counts (without norm)
  df=counts.final.N
  head(df)
  y <- DGEList(counts=
                 df[,setdiff(colnames(df), 
                             c("feature","event", "locus", "locus_overlap", "symbol", "gene_coordinates", "start", "end", "length"))] ,
               group=group)
  colnames(y)=rownames(targets)
  y$samples$lib.size = colSums(y$counts)
  design=model.matrix(~group, data=y$samples) # Test everything against LL24
  colnames(design)<-levels(y$samples$group)
  y<-calcNormFactors(y)
  y<-estimateDisp(y, design = design, robust = TRUE)
  fit <- glmQLFit(y, design = design, robust = TRUE)
  #end norm
}

qlf <- glmQLFTest(fit, coef = coef)
fdr<-p.adjust(qlf$table$PValue, method="BH") # Calculate FDR
tt<-cbind(qlf$table, fdr)
final<-list(qlf, tt)
names(final)=c("qlf", "full")
dim(final$qlf)

exon.intron_qlf <- dim(final$qlf)[1]
exon.intron.index_qlf <- which(final$full$PValue < 0.05 & final$full$fdr < 0.1) # 84155 / 177619 bins
save(exon.intron.index_qlf , file = "exon.intron.index_qlf.RData")
exon.intron <- final
exon.intron.fv <- exon.intron$qlf$fitted.values

# Calculate normalized logCPM data
exon.intron.logcounts2 <- cpm(y, normalized.lib.sizes = TRUE, log=TRUE)
save(exon.intron.logcounts2, file = "exon.intron.logcounts2.RData")

# Calculate normalized counts (cpm)
exon.intron.norm_counts <- cpm(y, normalized.lib.sizes = TRUE) 
save(exon.intron.norm_counts, file = "exon.intron.norm_counts.RData")

# Save all exon.intron count data
write.table(qlf$fitted.values[exon.intron.index_qlf,], file="circadian.bins.fitted.values.tab", sep="\t", row.names = TRUE, col.names = TRUE, na = "NA")
write.table(exon.intron.counts , file="circadian.bins.raw.values.tab", sep="\t", row.names = TRUE, col.names = TRUE, na = "NA")
write.table(df, file="circadian.bins.filtered.values.tab", sep="\t", row.names = TRUE, col.names = TRUE, na = "NA")
write.table(exon.intron.logcounts2[exon.intron.index_qlf,], file="circadian.bins.logCPM.normalised.values.tab", sep="\t", row.names = TRUE, col.names = TRUE, na = "NA")
write.table(exon.intron.norm_counts[exon.intron.index_qlf,], file="circadian.bins.CPM.normalised.values.tab", sep="\t", row.names = TRUE, col.names = TRUE, na = "NA")
write.table(exon.intron.logcounts2, file="full.circadian.bins.logCPM.normalised.values.tab", sep="\t", row.names = TRUE, col.names = TRUE, na = "NA")
write.table(exon.intron.norm_counts, file="full.circadian.bins.CPM.normalised.values.tab", sep="\t", row.names = TRUE, col.names = TRUE, na = "NA")

# Add annotation data to normalised counts and logCPM data
event.md <- read.csv("event.metadata.tab", na.strings=c("", "NA"), sep = "\t", header = TRUE, stringsAsFactors = FALSE)
exon.intron.logcounts2.annot <- getEventInfo(exon.intron.logcounts2, event.md)
exon.intron.norm_counts.annot <- getEventInfo(exon.intron.norm_counts, event.md)
write.table(logcounts2.annot[exon.intron.index_qlf,], file="annot.circadian.bins.logCPM.normalised.values.tab", sep="\t", row.names = TRUE, col.names = TRUE, na = "NA")
write.table(norm_counts.annot[exon.intron.index_qlf,], file="annot.circadian.bins.CPM.normalised.values.tab", sep="\t", row.names = TRUE, col.names = TRUE, na = "NA")

#####################################################
#       Splicing Circadian Analysis - JTK Cycle     #
#####################################################

# Set basedir
if (get_os() != "windows") {basedir <- "/data/projects/circadianLL/comp" 
} else {basedir <- "c:/circadianLL/comp/"}
setwd(basedir)

# Load JTK Cycle functions
source("JTKversion3/JTK_CYCLEv3.1.R")
# load("exon.intron.norm_counts.RData")
# load("exon.intron.index_qlf.RData")

#### Gene analysis ####
project <- "exon.intron.LL"
options(stringsAsFactors=FALSE)
annot <- rownames(exon.intron.norm_counts[exon.intron.index_qlf,])
data <- exon.intron.norm_counts[exon.intron.index_qlf,] # 84,155 bins that passed the rd filter, and the GLM + QL F-test (p<0.05 and q<0.10)
jtkdist(ncol(data))


jtkdist(12,2) # 12 total time points, 2 replicates per time point
periods <- 5:7  # looking for rhythms between 20-28 hours (i.e. between 5 and 7 time points per cycle).
jtk.init(periods, 4) # sampling rate in hours

cat("JTK analysis started on",date(),"\n")
flush.console()

st <- system.time({
  res <- apply(data,1,function(z) {
    jtkx(z)
    c(JTK.ADJP,JTK.PERIOD,JTK.LAG,JTK.AMP)
  })
  res <- as.data.frame(t(res))
  bhq <- p.adjust(unlist(res[,1]),"BH")
  res <- cbind(bhq,res)
  colnames(res) <- c("BH.Q","ADJ.P","PER","LAG","AMP")
  results <- cbind(annot,res,data)
  results <- results[order(res$ADJ.P,-res$AMP),]
})
print(st)
save(results,file=paste("JTK/JTK",project,"rda",sep="."))
write.table(results,file=paste("JTK/JTK",project,"txt",sep="."),row.names=T,col.names=T,quote=F,sep="\t")
write.table(getEventInfo(results,event.md),file=paste("JTK/annot.JTK",project,"txt",sep="."),row.names=T,col.names=T,quote=F,sep="\t")

# recalculate phases to be circadian  
results2 <-data.frame(cbind(results[,1:5],
                            recalculated.phase = results$LAG/(results$PER/24),
                            results[,6:30]))
save(results2, file = paste("JTK/JTK_recalculated_",project,"rda",sep="."))
write.table(results2,file=paste("JTK/JTK_recalculated_",project,"txt",sep="."),row.names=T,col.names=T,quote=F,sep="\t")
write.table(getEventInfo(results2,event.md),file=paste("JTK/annot.JTK_recalculated_",project,"txt",sep="."),row.names=T,col.names=T,quote=F,sep="\t")

###########################################################
# Plot Circadian Splicing evenets heatmaps and histograms #
###########################################################
if (get_os() != "windows") {basedir <- "/data/projects/circadianLL/comp/" 
} else {basedir <- "c:/circadianLL/comp"}
setwd(basedir)

load(file = "JTK/JTK.exon.intron.LL.rda")
load(file = "JTK/JTK_recalculated_.exon.intron.LL.rda")

# check p values cut offs
length(which(results$ADJ.P<0.001 & results$BH.Q < 0.001 & results$PER != 0)) #1227
length(which(results$ADJ.P<0.001 & results$BH.Q < 0.01 & results$PER != 0)) #4460
length(which(results$ADJ.P<0.001 & results$BH.Q < 0.05 & results$PER != 0)) #5306
length(which(results$ADJ.P<0.01 & results$BH.Q < 0.01 & results$PER != 0)) #4460
length(which(results$ADJ.P<0.01 & results$BH.Q < 0.05 & results$PER != 0)) #10786
length(which(results$ADJ.P<0.05 & results$BH.Q < 0.05 & results$PER != 0)) #10786

# filter clock controlled events
results <- results[which(results$ADJ.P<0.01 & results$BH.Q < 0.01 & results$PER != 0),] # 4460 events with p < 0.01 and q < 0.01
results2 <- results2[which(results2$ADJ.P<0.01 & results2$BH.Q < 0.01 & results2$PER != 0),] # 4460 events with p < 0.01 and q < 0.01

results <- results[order(results$LAG, results$AMP),] # Tau phases
results2 <- results2[order(results2$recalculated.phase, results$AMP),] # Circadian phases
indexEvents <- match(rownames(results), rownames(exon.intron.logcounts2))

# Create a scaled logcounts to use with the image.scale function
scaled_exon.intron.logcounts2 <- t(scale(t(as.matrix(exon.intron.logcounts2[indexEvents,]))))
range(scaled_exon.intron.logcounts2)

# Create a cluster type heatmap of gene expression
# create the auxiliary data needed
cluster <- unique(results$LAG)
j <- NULL
for (i in seq_along(cluster)){
  print(i)
  j[i] <- length(which(results$LAG == cluster[i]))
} 
size <- j

# Plot the cluster
plot.clust(list(data = exon.intron.logcounts2[indexEvents,],
                      cluster = results$LAG,
                      size = size),
                 standardize.gene = TRUE,
                 gene_labels = FALSE)
dev.off()

# Save the cluster heatmap
png(file="heatmaps/circadian.bins.heatmap.png",    # create PNG for the MDS        
    width = 6.5*1200,        # 6.5 x 1200 pixels
    height = 5*1200,
    res = 1200,            # 1200 pixels per inch
    pointsize = 10)        # smaller font size
plot.clust(list(data = exon.intron.logcounts2[indexEvents,],
                      cluster = results$LAG,
                      size = size),
                 standardize.gene = TRUE,
                 gene_labels = FALSE)
dev.off()

# Create a scale for the plot
png(file="heatmaps/circadian.bins.heatmap.scale.png",    # create PNG for the MDS        
    width = 6*600,        # 6.5 x 1200 pixels
    height = 2*600,
    res = 600,            # 1200 pixels per inch
    pointsize = 10)        # smaller font size
image.scale(z = scaled_exon.intron.logcounts2, zlim = range(scaled_exon.intron.logcounts2), col = colorRampPalette(c("blue", "gray","yellow"))(64))
dev.off()

# Load annotated exon.intron file
annot.exon.intron <- read.csv(file = "JTK/annot.JTK_recalculated_.exon.intron.LL.txt", stringsAsFactors = FALSE, header = TRUE, row.names = 1, sep = "\t")

save.image(file = "myEnvironment.RData")
# load(file = "myEnvironment.RData")
sessionInfo()

