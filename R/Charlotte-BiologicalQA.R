#' ---
#' title: "Biological QA (mod'ed from Charlotte Soneson's)"
#' author: "Nicolas Delhomme, mod'ed from Charlotte Soneson's"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#' ---
#' # Introduction
#' In this tutorial we walk through a gene-level RNA-seq differential expression analysis using Bioconductor packages. 
#' We start from the gene-vs-sample count matrix, and thus assume that the raw reads have already been quality controlled 
#' and that the gene expression has been quantified (either using alignment and counting, or by applying an alignment-free 
#' quantification tool). We perform exploratory data analysis (EDA) for quality assessment and to explore the relationship 
#' between samples, then perform differential gene expression analysis, and visually explore the results.
#' 
#' Bioconductor has many packages supporting analysis of high-throughput sequence data, including RNA-seq. 
#' The packages that we will use in this tutorial include core packages maintained by the [Bioconductor core team](https://www.bioconductor.org/about/core-team/) 
#' for importing and processing raw sequencing data and loading gene annotations. We will also use contributed packages 
#' for statistical analysis and visualization of sequencing data. Through scheduled releases every 6 months, the Bioconductor 
#' project ensures that all the packages within a release will work together in harmony (hence the “conductor” metaphor). 
#' The packages used in this tutorial are loaded with the library function and can be installed by following the [Bioconductor 
#' package installation instructions](http://bioconductor.org/install/#install-bioconductor-packages).
#' 
#' Many parts of this tutorial are based on parts of a published RNA-seq workflow available via 
#' [F1000Research](http://f1000research.com/articles/4-1070) (Love et al. 2015) 
#' and as a [Bioconductor package](https://www.bioconductor.org/packages/release/workflows/html/rnaseqGene.html).
#' 
#' ## Citing scientific research software
#' 
#' If you use the results from an R package in published research, you can find the proper citation for the software 
#' by typing `citation("pkgName")`, where you would substitute the name of the package for `pkgName`. Citing methods papers
#'  helps to support and reward the individuals who put time into open source software for genomic data analysis.
#' 
#' ## Experimental data
#' 
#' The data used in this workflow comes from an RNA-seq experiment conducted in _Norway spruce_ investigating
#' the zygotic embryo development whereby a seed constituting of two tissues, the embryo and it's supporting megagametophyte, 
#' develops in to a plantlet. Out of the eight time points available in the original (unpublished as of yet, but hopefully soon) studies,
#' we have selected two: **B4** and **B8**. Hence, we have 12 samples, three biological replicates per time point and tissue 
#' (**ZE** or **FMG**, zygotic and female megagametophyte, respectively). The rationale for the studies is to compare it to the
#' process used by the forestry industry, called _somatic embryogenesis_ that is used to regenerate an embryo from a somatic tissue.
#' 
#' We start by setting the path to the folder containing the data that will be used in the workflow. 
datadir <- "/home/ubuntu/raw_data/exploratoryDataAnalysis"
#' ```{r dir,echo=FALSE}
#' datadir <- "/mnt/picea/home/delhomme/Git/RNA-Seq-workflow-Differential-Expression/data"
#' ```

#' That directory contains gene- and transcript-level quantifications for the samples in this experiment, 
#' as well as a metadata table indicating the identity of the samples. All analyses are based on the Norway spruce reference genome (v1.0) 
#' and the corresponding annotation. The types of available quantifications are:
#'   
#' * featureCounts: gene-level read counts obtained by featureCounts (Liao, Smyth, and Shi 2014) 
#' from the [Rsubread](https://bioconductor.org/packages/3.12/Rsubread) package, following alignment by STAR (Dobin et al. 2013).
#' * salmon: transcript-level quantifications (read counts and TPMs) obtained by Salmon (Patro et al. 2017).
#' 
#' The metadata.txt file contains the sample information.
#' 
#' ## Goal of this tutorial
#' 
#' Our goal in this tutorial is to bring a summary of the RNA-seq experiment into R/Bioconductor for visualization and statistical testing. 
#' We want to visualize the relationships between the samples (within and across the treatment), and then we want to perform statistical
#'  tests to find which genes are changing their expression due to treatment.
#' 
#' # Setup
#' 
#' Here we load the needed libraries we will need for the analysis and set some graphical parameters.
#' * Libraries
suppressPackageStartupMessages({
  library(DESeq2)
  library(edgeR)
  library(LSD)
  library(RColorBrewer)
  library(SummarizedExperiment)
  library(tidyverse)
  library(tximport)
  library(vsn)
})

#' * Graphics
pal <- brewer.pal(8,"Dark2")

#' # Reading the metadata
#' 
#' First, we will read the metadata for the experiment. The two annotations of primary interest for this tutorial is `Tissue`, which 
#' represents the tissue, and `Time`, which indicates the developmental time point at which the sample was collected. 
#' There are 3 biological replicates per `Tissue` and `Time`. The sample identifier is given by the names column, and will be used to match the metadata table to the quantifications.
meta <- read_tsv(file.path(datadir, "metadata.txt"),col_types=cols(.default=col_factor()))
meta

#' What is a factor?
#' 
#' The representation in R of a categorical variable.
meta$Tissue
levels(meta$Tissue)
as.integer(meta$Tissue)

#' # Summarizing an RNA-seq experiment as a count matrix

#' Count-based statistical methods such as [DESeq2](https://bioconductor.org/packages/3.12/DESeq2) (Love, Huber, and Anders 2014), 
#' [edgeR](https://bioconductor.org/packages/3.12/edgeR) (Robinson, McCarthy, and Smyth 2009), 
#' [limma](https://bioconductor.org/packages/3.12/limma) with the voom method (Law et al. 2014), 
#' [DSS](https://bioconductor.org/packages/3.12/DSS) (Wu, Wang, and Wu 2013), [EBSeq](https://bioconductor.org/packages/3.12/EBSeq) (Leng et al. 2013), 
#' [BaySeq](https://bioconductor.org/packages/3.12/BaySeq) (Hardcastle and Kelly 2010) and [DEXSeq](https://bioconductor.org/packages/3.12/DEXSeq)
#'  (Anders, Reyes, and Huber 2012) expect input data as obtained, e.g., from RNA-seq or another high-throughput
#'  sequencing experiment in the form of a matrix of integer values, or “counts”. The value in the i-th row and the j-th column 
#'  of the matrix tells how many reads (or fragments, for paired-end RNA-seq) have been assigned to feature i in sample j. 
#'  For RNA-seq, a feature is typically a gene, a transcript or an exon. Analogously, for other types of assays, the rows of the 
#'  matrix might correspond e.g., to binding regions (with ChIP-Seq), species of bacteria (with metagenomic datasets), or 
#'  peptide sequences (with quantitative mass spectrometry).
#'
#' The fact that the values in the matrix are counts of sequencing reads (in the case of single-end sequencing) or fragments
#'  (for paired-end sequencing) is important for the count-based statistical models, e.g. DESeq2 or edgeR, as only the counts 
#'  allow assessing the measurement precision correctly. It is important to never provide counts that have been normalized 
#'  for sequencing depth/library size to these packages, as the statistical model is most powerful when applied to counts, 
#'  and is designed to account for library size differences internally.
#'
#' An alternative to using actual counts of reads or fragments aligned to the genome is to use estimated counts from software 
#' that use pseudo-alignment to the transcriptome. Since these represent expected counts rather than observed counts they are 
#' not necessarily integers, and thus may need to be rounded before they are fed to the count-based pipelines.
#'
#' In the sections below, we will show how to generate gene-level count matrices in R from the output of two of 
#' the most common quantification pipelines (featureCounts and Salmon).

#' ## featureCounts
#' The `featureCounts` function (from the [Rsubread](https://bioconductor.org/packages/3.12/Rsubread) package) takes as input 
#' bam files resulting from read alignment to the genome, and counts the number of reads overlapping each genomic feature 
#' (here, each gene). For the purposes of this tutorial, `featureCounts` was run as follows (**DON’T RUN THIS**):
#' ```{r subread, eval=FALSE}  
#' library(Rsubread)
#' fc <- featureCounts(files = files, 
#'                    annot.ext = gtf, 
#'                    isGTFAnnotationFile = TRUE,
#'                    GTF.featureType = "exon", 
#'                    GTF.attrType = "gene_id", 
#'                    useMetaFeatures = TRUE, 
#'                    strandSpecific = 0, 
#'                    isPairedEnd = TRUE, 
#'                    nthreads = 6)
#' ```
#'
#' where `files` is a vector of file names pointing to the bam files for the different samples, and `gtf` points to a gtf 
#' file with the genomic regions corresponding to each gene. `featureCounts` returns a list with several element, one of 
#' which is the estimated count matrix. For simplicity, we saved the list output from the command above to a file and 
#' here we just load it back into R.

fc <- readRDS(file.path(datadir, "featureCounts/star_featureCounts.rds"))
names(fc)
counts_featurecounts <- fc$counts
colnames(counts_featurecounts) <- sub("_S.*","",colnames(counts_featurecounts))
head(counts_featurecounts)
dim(counts_featurecounts)
fc$stat

#' Let us visualise it
stats <- as.matrix(fc$stat %>% column_to_rownames("Status"))
colnames(stats) <- sub("_S.*","",colnames(stats))
barplot(stats[rowSums(stats)>0,],
        beside = TRUE,col = pal[1:3],las=2,
        cex.axis = .8,cex.names = .8,
        legend.text = rownames(stats)[rowSums(stats)>0],
        args.legend = list(horiz=TRUE,cex=.8,x="top",bty="n"))

#' ## Alignment-free quantification
#'
#' Alignment-free transcript quantification software such as kallisto (Bray et al. 2016), Salmon (Patro et al. 2017) and 
#' Sailfish (Patro, Mount, and Kingsford 2014), as well as other transcript quantification methods like Cufflinks 
#' (Trapnell et al. 2010, 2013) and RSEM (Li and Dewey 2011), differ from the counting methods covered above in that they 
#' provide quantifications (usually both as counts and as TPMs) for each transcript. These can then be summarized on the 
#' gene level by adding all values for transcripts from the same gene. A simple way to import results from these packages 
#' into R is provided by the [tximport](https://bioconductor.org/packages/3.12/tximport) and 
#' [tximeta](https://bioconductor.org/packages/3.12/tximeta) packages. Here, `tximport` reads the quantifications into a list
#'  of matrices, and `tximeta` aggregates the information into a `SummarizedExperiment` object, and also automatically adds
#'   additional annotations for the features. Both packages can return quantifications on the transcript level or aggregate 
#'   them on the gene level. They also calculate average transcript lengths for each gene and each sample, which can be used 
#'   as offsets to improve the differential expression analysis by accounting for differential isoform usage across 
#'   samples (Soneson, Love, and Robinson 2015).

#' ### Salmon
#' #' The code below imports the Salmon quantifications into R using the `tximport` package. Using the `txmeta` package would require
#' us to work with an organism for which the annotation would be available from the bioconductor project (unlikely for non-model organism).

# List all quant.sf output files from Salmon
salmonfiles <- file.path(datadir, "salmon", meta$names,"quant.sf")
names(salmonfiles) <- meta$names
stopifnot(all(file.exists(salmonfiles)))

#' The Norway spruce annotation have only one transcript per gene, as such there is not need to summarise the expression at the gene
#' level, which would require the use of the `tx2gene` argument instead of `txOut`. The `tx2gene` expects a tab delimited file with the
#' transcript ID and the corresponding gene ID it two columns.
sg <- tximport(files=salmonfiles,type="salmon",txOut=TRUE)

#' Note that _Salmon_ returns _estimated_ counts, which are not necessarily integers. They may need to be rounded before they are passed 
#' to count-based statistical methods (e.g. _DESeq2_). To obtain consistent results with different pipelines, we round the estimated counts 
#' here and use the resulting matrix as input also to edgeR.

counts_salmon <- round(sg$counts)
rownames(counts_salmon) <- sub("\\.1$","",rownames(counts_salmon))

#' ## Comparison of counts
#' For illustration, we compare the counts obtained by the two quantification approaches for the first sample. As we can see, there is an overall 
#' good correspondence between the different methods, especially for the genes with high counts.
spl <- "P11562_110"

#' for comparing, we need to make sure, both counts table are
#' sorted in the same way. We create a data.frame containing 
#' both sorted according to their gene IDs
gns <- sort(intersect(rownames(counts_salmon),rownames(counts_featurecounts)))
quants <- data.frame(featureCounts = counts_featurecounts[gns, spl],
                     salmon = counts_salmon[gns, spl])

head(quants[order(quants$salmon,decreasing=TRUE),])

#' Let us plot the raw counts
heatscatter(quants$featureCounts,quants$salmon,
            xlab="featureCounts",ylab="salmon",
            main="scatterplot (raw counts)")
abline(0,1,lty=2,col="black")

#' Let us go onto a log scale
heatscatter(log10(quants$featureCounts+1),log10(quants$salmon+1),
            xlab="featureCounts",ylab="salmon",
            main="scatterplot (log10 counts + 1)")
abline(0,1,lty=2,col="black")

#' # Representing counts for differential expression packages
#' 
#' At this point, we have a gene-level count matrix. In the rest of this tutorial, we will work with the counts generated by `Salmon` and imported into R with `tximport`. 
#' This is a branching point where we could use a variety of Bioconductor packages for exploration and differential expression of the count matrix, including 
#' [edgeR](https://bioconductor.org/packages/3.12/edgeR) (Robinson, McCarthy, and Smyth 2009), [DESeq2](https://bioconductor.org/packages/3.12/DESeq2) (Love, Huber, and Anders 2014), 
#' [limma](https://bioconductor.org/packages/3.12/limma) with the voom method (Law et al. 2014), [DSS](https://bioconductor.org/packages/3.12/DSS) (Wu, Wang, and Wu 2013), 
#' [EBSeq](https://bioconductor.org/packages/3.12/EBSeq) (Leng et al. 2013) and [BaySeq](https://bioconductor.org/packages/3.12/BaySeq) (Hardcastle and Kelly 2010). 
#' We will continue using `DESeq2` and `edgeR`. #' Why? Because you are likely to encounter both in the literature.
#' 
#' Differences? Little - mostly taste. `DESeq2` abstracts a lot 
#' of the decisions from the user, while `edgeR` leave these all
#' to the user.

#' 
#' Bioconductor software packages often define and use a custom class for storing data that makes sure that all the needed data slots are consistently provided 
#' and fulfill any requirements. In addition, Bioconductor has general data classes (such as the `SummarizedExperiment`) that can be used to move data between packages. 
#' The [DEFormats](https://bioconductor.org/packages/3.12/DEFormats) package can be useful for converting between different classes. The core Bioconductor classes also 
#' provide useful functionality: for example, subsetting or reordering the rows or columns of a `SummarizedExperiment` automatically subsets or reorders the associated 
#' `rowRanges` and `colData`, which can help to prevent accidental sample swaps that would otherwise lead to spurious results. With `SummarizedExperiment` this is all 
#' taken care of behind the scenes.
#' 
#' Each of the packages we will use for differential expression has a specific class of object used to store the summarization of the RNA-seq experiment and the intermediate 
#' quantities that are calculated during the statistical analysis of the data. `DESeq2` uses a `DESeqDataSet` and `edgeR` uses a `DGEList`.
#' 
#' ## The DESeqDataSet, sample information, and the design formula
#' 
#' In `DESeq2`, the custom class is called `DESeqDataSet`. It is built on top of the `SummarizedExperiment` class, and it is easy to convert `SummarizedExperiment`
#'  objects into `DESeqDataSet` objects. One of the two main differences compared to a `SummarizedExperiment` object is that the `assay` slot is instead accessed 
#'  using the `counts` accessor function, and the `DESeqDataSet` class enforces that the values in this matrix are non-negative integers.
#'  
#'  A second difference is that the `DESeqDataSet` has an associated _design formula_. The experimental design is specified at the beginning of the analysis, 
#'  as it will inform many of the `DESeq2` functions how to treat the samples in the analysis (one exception is the size factor estimation, i.e., the adjustment 
#'  for differing library sizes, which does not depend on the design formula). The design formula tells which columns in the sample information table (`colData`)
#'   specify the experimental design and how these factors should be used in the analysis.
#'   
#'   Let’s remind ourselves of the design of our experiment:

#' Let us take a look at the information we have about
#' our samples to build our model
meta <- meta %>% column_to_rownames("names")
meta

#' Time and Tissue are the two categorical variables of interest
meta$Time
meta$Tissue

#' Both are factors (categorical variables), and based on the 
#' order of the levels, B4 and FMG will be the reference we 
#' compare to (_i.e._ B4 FMG will be our expression baseline)

#' Some sanity checking
stopifnot(all(colnames(counts_salmon) == rownames(meta)))

#' Let us create the DESeq objects. Using the method `DESeqDataSetFromTximport` make sure to take
#' into account the observed transcript length as a representation of a potential differential isoform 
#' usage (_i.e._ setting an offset), which is the point of the tximport methods (Soneson, Love, and Robinson 2015) for gene-level analysis. 
#' See [there](https://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html#downstream-dge-in-bioconductor).
#' 
ds_se <- DESeqDataSetFromTximport(txi=sg, 
                                colData=meta ,
                                design = ~ Tissue * Time)

#' Before proceeding, we know that sample `P11562_112` is undersampled. to prevent it from biasing our further processing we remove it, 
#' as well as the corresponding metadata.
sel <- rownames(meta) != "P11562_112"
ds_se <- ds_se[,sel]
counts_salmon <- counts_salmon[,sel]
meta <- meta[sel,]

#' ## The DGEList
#' 
#' As mentioned above, the edgeR package uses another type of data container, namely a `DGEList` object. 
#' It is just as easy to create a `DGEList` object using a count matrix and a table with information 
#' about samples. We can additionally add information about the genes. Here we illustrate how to 
#' generate a `DGEList` object from the Salmon count matrix.

#' Create the object
genetable <- data.frame(gene.id = rownames(counts_salmon),
                        stringsAsFactors = FALSE)
stopifnot(all(rownames(meta) == colnames(counts_salmon)))
dge <- DGEList(counts = counts_salmon, 
               samples = meta, 
               genes = genetable)
names(dge)

#' Just like the `SummarizedExperiment` and the `DESeqDataSet` the `DGEList` contains all the information we need: 
#' the count matrix, information about the samples (the columns of the count matrix), and information about the 
#' genes (the rows of the count matrix). One difference compared to the `DESeqDataSet` is that the experimental 
#' design is not defined when creating the `DGEList`, but later in the workflow.
#' 
#' To include information about the average transcript lengths (offsets) estimated by `tximport` in `edgeR`, 
#' the offsets must be added manually. The `tximport` vignette shows how this is done, and we repeat it here.
#' 
#' Here is an illustration of what I meant earlier, `edgeR` leaves a lot for the user to do
#' 
#' Here we retrieve the average transcript length, use
#' it to correct the salmon expression values, before
#' calculating the library size factor (_i.e._ the difference
#' in sequencing depth between samples). 
avetxlengths <- sg$length[,sel]
rownames(avetxlengths) <- sub("\\.1$","",rownames(avetxlengths))
stopifnot(all(rownames(avetxlengths) == rownames(counts_salmon)))
stopifnot(all(colnames(avetxlengths) == colnames(counts_salmon)))
avetxlengths <- avetxlengths/exp(rowMeans(log(avetxlengths)))
offsets <- log(calcNormFactors(counts_salmon/avetxlengths)) + 
  log(colSums(counts_salmon/avetxlengths))
dge <- scaleOffset(dge, t(t(log(avetxlengths)) + offsets))
names(dge)

#' Once a DGEList has been created, we calculate between-sample (TMM) normalization factors, using the `calcNormFactors` function in edgeR.
dge <- edgeR::calcNormFactors(dge)
dge$samples

#' Let us visualise the library size factor differences
#' 
#' It is very similar across all samples varying from -25% to +15% 
boxplot(dge$samples$norm.factors)

#' In the remainder of this tutorial, we will use the `ds_se` and `dge` objects for `DESeq2` and `edgeR` analyses, respectively.
#'
#' # Exploratory analysis and visualization
#' There are two separate analysis paths in this tutorial:
#' 
#' 1. visual exploration of sample relationships, in which we will discuss transformation of the counts for computing distances or making plots
#' 2. statistical testing for differences attributable to treatment, controlling for cell line effects
#' 
#' Importantly, **the statistical testing methods rely on original count data (not scaled or transformed)** for calculating the precision of measurements. 
#' However, for visualization and exploratory analysis, transformed counts are typically more suitable. Thus, it is critical to separate the two workflows 
#' and use the appropriate input data for each of them.

#' ## Transformations
#' Many common statistical methods for exploratory analysis of multidimensional data, for example clustering and _principal components analysis_ 
#' (PCA), work best for data that generally has the same range of variance at different ranges of the mean values. When the expected amount 
#' of variance is approximately the same across different mean values, the data is said to be _homoskedastic_. For RNA-seq raw counts, however, 
#' the variance grows with the mean. For example, if one performs PCA directly on a matrix of size-factor-normalized read counts, the result 
#' typically depends only on the few most strongly expressed genes because they show the largest absolute differences between samples. A simple 
#' and often used strategy to avoid this is to take the logarithm of the normalized count values plus a small pseudocount; however, now the 
#' genes with the very lowest counts will tend to dominate the results because, due to the strong Poisson noise inherent to small count values, 
#' and the fact that the logarithm amplifies differences for the smallest values, these low count genes will show the strongest relative differences between samples.
#' 
#' As a solution, `DESeq2` offers transformations for count data that stabilize the variance across the mean: the `regularized logarithm` (rlog) 
#' and the `variance stabilizing transformation` (VST). These have slightly different implementations, discussed a bit in the `DESeq2` paper and 
#' in the vignette, but a similar goal of stabilizing the variance across the range of values. Both produce log2-like values for high counts. 
#' Here we will use the variance stabilizing transformation implemented with the `vst` function.
#' 
#' Neither `DESeq2` nor `edgeR` modify the count data. 
#' They calculate the parameters (such as the library size factor,
#' the dispersion, etc.) to use in their model. If we want to 
#' visualise the data in a reduced dimensionality (not 50+ thousands genes 
#' and x samples as a matrix), such as a PCA (principal component analysis),
#' or MDS (Multi-dimensional scaling), we need to normalise the data.
#' 
#' In addition to a difference in sequencing depth between samples, RNA-Seq
#' data suffers from another problem, a mean-variance relationship, _i.e._ the
#' data is heteroskedastic.
meanSdPlot(log(assay(ds_se)[rowSums(assay(ds_se))>30,]))

#' To counteract this, one can use a variance stabilising transformation
#' (VST), a heuristic that will transform the variance so it becomes 
#' independent of the mean
vsd <- DESeq2::vst(ds_se)

#' Et voila! Not perfect but good enough for visualisation - Actually a terrible plot here, as we have so few data points. I'll show case a better one.
meanSdPlot(log(assay(vsd)[rowSums(assay(vsd))>0,]))

#' ## PCA plot

#' One way to visualize sample-to-sample distances is a principal components analysis (PCA). 
#' In this ordination method, the data points (here, the samples) are projected onto the 2D 
#' plane such that they spread out in the two directions that explain most of the differences (Figure below). 
#' The x-axis (the first principal component, or PC1) is the direction that separates the data points the most 
#' (_i.e._, the direction with the largest variance). The y-axis (the second principal component, or PC2) 
#' represents the direction with largest variance subject to the constraint that it must be orthogonal to the 
#' first direction. The percent of the total variance that is contained in the direction is printed in the axis
#'  label. Note that these percentages do not sum to 100%, because there are more dimensions that contain the 
#'  remaining variance (although each of these remaining dimensions will explain less than the two that we see).
#' 
DESeq2::plotPCA(vsd, intgroup = "Tissue")
DESeq2::plotPCA(vsd, intgroup = "Time")

#' ## MDS plot

#' Another way to reduce dimensionality, which is in many ways similar to PCA, is multidimensional scaling (MDS). 
#' For MDS, we first have to calculate all pairwise distances between our objects (samples in this case), and then 
#' create a (typically) two-dimensional representation where these pre-calculated distances are represented as 
#' accurately as possible. This means that depending on how the pairwise sample distances are defined, the two-dimensional 
#' plot can be very different, and it is important to choose a distance that is suitable for the type of data at hand.
#'
#' `edgeR` contains a function `plotMDS`, which operates on a `DGEList` object and generates a two-dimensional MDS 
#' representation of the samples. The default distance between two samples can be interpreted as the “typical” 
#' log fold change between the two samples, for the genes that are most different between them (by default, the 
#' top 500 genes, but this can be modified). We generate an MDS plot from the DGEList object dge, coloring by 
#' the treatment and using different plot symbols for different cell lines.
plotMDS(dge, top = 500, 
        labels = NULL, 
        col = as.numeric(dge$samples$Time), 
        pch = as.numeric(dge$samples$Tissue), 
        cex = 2, gene.selection = "common")

#' # Session Info
#' ```{r, session info, echo=FALSE}
#' sessionInfo()
#' ```
#' # References
#' Anders, Simon, Alejandro Reyes, and Wolfgang Huber. 2012. “Detecting Differential Usage of Exons from RNA-seq Data.” Genome Res. 22 (10): 2008–17.
#' 
#' Benjamini, Yoav, and Yosef Hochberg. 1995. “Controlling the False Discovery Rate: A Practical and Powerful Approach to Multiple Testing.” Journal of the Royal Statistical Society. Series B (Methodological) 57 (1): 289–300. [http://www.jstor.org/stable/2346101](http://www.jstor.org/stable/2346101).
#'
#' Bray, Nicolas L, Harold Pimentel, Páll Melsted, and Lior Pachter. 2016. “Near-Optimal RNA-Seq Quantification.” Nat. Biotechnol.
#'
#' Dobin, Alexander, Carrie A. Davis, Felix Schlesinger, Jorg Drenkow, Chris Zaleski, Sonali Jha, Philippe Batut, Mark Chaisson, and Thomas R. Gingeras. 2013. “STAR: ultrafast universal RNA-seq aligner.” Bioinformatics 29 (1). Oxford University Press: 15–21. [https://doi.org/10.1093/bioinformatics/bts635](https://doi.org/10.1093/bioinformatics/bts635).
#'
#' Dudoit, Sandrine, Yee H. Yang, Matthew J. Callow, and Terence P. Speed. 2002. “Statistical methods for identifying differentially expressed genes in replicated cDNA microarray experiments.” In Statistica Sinica, 111–39. [http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.117.9702](http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.117.9702).
#'
#' Hardcastle, Thomas, and Krystyna Kelly. 2010. “baySeq: Empirical Bayesian methods for identifying differential expression in sequence count data.” BMC Bioinformatics 11 (1): 422+. [https://doi.org/10.1186/1471-2105-11-422](https://doi.org/10.1186/1471-2105-11-422).
#'
#' Himes, Blanca E., Xiaofeng Jiang, Peter Wagner, Ruoxi Hu, Qiyu Wang, Barbara Klanderman, Reid M. Whitaker, et al. 2014. “RNA-Seq transcriptome profiling identifies CRISPLD2 as a glucocorticoid responsive gene that modulates cytokine function in airway smooth muscle cells.” PloS One 9 (6). [https://doi.org/10.1371/journal.pone.0099625](https://doi.org/10.1371/journal.pone.0099625).
#' 
#' Law, Charity W., Yunshun Chen, Wei Shi, and Gordon K. Smyth. 2014. “Voom: precision weights unlock linear model analysis tools for RNA-seq read counts.” Genome Biology 15 (2). BioMed Central Ltd: R29+. [https://doi.org/10.1186/gb-2014-15-2-r29](https://doi.org/10.1186/gb-2014-15-2-r29).
#'
#' Leng, N., J. A. Dawson, J. A. Thomson, V. Ruotti, A. I. Rissman, B. M. G. Smits, J. D. Haag, M. N. Gould, R. M. Stewart, and C. Kendziorski. 2013. “EBSeq: an empirical Bayes hierarchical model for inference in RNA-seq experiments.” Bioinformatics 29 (8). Oxford University Press: 1035–43. [https://doi.org/10.1093/bioinformatics/btt087](https://doi.org/10.1093/bioinformatics/btt087).
#'
#' Li, Bo, and Colin N. Dewey. 2011. “RSEM: accurate transcript quantification from RNA-Seq data with or without a reference genome.” BMC Bioinformatics 12: 323+. [https://doi.org/10.1186/1471-2105-12-3231](https://doi.org/10.1186/1471-2105-12-3231).
#'
#' Liao, Y., G. K. Smyth, and W. Shi. 2014. “featureCounts: an efficient general purpose program for assigning sequence reads to genomic features.” Bioinformatics 30 (7). Oxford University Press: 923–30. [https://doi.org/10.1093/bioinformatics/btt656](https://doi.org/10.1093/bioinformatics/btt656).
#'
#' Love, Michael I., Simon Anders, Vladislav Kim, and Wolfgang Huber. 2015. “RNA-Seq Workflow: Gene-Level Exploratory Analysis and Differential Expression.” F1000Research, October. [https://doi.org/10.12688/f1000research.7035.1](https://doi.org/10.12688/f1000research.7035.1).
#'
#' Love, Michael I., Wolfgang Huber, and Simon Anders. 2014. “Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2.” Genome Biology 15 (12). BioMed Central Ltd: 550+. [https://doi.org/10.1186/s13059-014-0550-8](https://doi.org/10.1186/s13059-014-0550-8).
#'
#' Patro, Rob, Geet Duggal, Michael I Love, Rafael A Irizarry, and Carl Kingsford. 2017. “Salmon Provides Fast and Bias-Aware Quantification of Transcript Expression.” Nat. Methods.
#'
#' Patro, Rob, Stephen M. Mount, and Carl Kingsford. 2014. “Sailfish enables alignment-free isoform quantification from RNA-seq reads using lightweight algorithms.” Nature Biotechnology 32: 462–64. [https://doi.org/10.1038/nbt.2862](https://doi.org/10.1038/nbt.2862).
#'
#' Robinson, M. D., D. J. McCarthy, and G. K. Smyth. 2009. “edgeR: a Bioconductor package for differential expression analysis of digital gene expression data.” Bioinformatics 26 (1). Oxford University Press: 139–40. [https://doi.org/10.1093/bioinformatics/btp616](https://doi.org/10.1093/bioinformatics/btp616).
#'
#' Soneson, Charlotte, Michael I. Love, and Mark Robinson. 2015. “Differential analyses for RNA-seq: transcript-level estimates improve gene-level inferences.” F1000Research 4 (1521). [https://doi.org/10.12688/f1000research.7563.1](https://doi.org/10.12688/f1000research.7563.1).
#'
#' Trapnell, Cole, David G Hendrickson, Martin Sauvageau, Loyal Goff, John L Rinn, and Lior Pachter. 2013. “Differential Analysis of Gene Regulation at Transcript Resolution with RNA-seq.” Nat. Biotechnol. 31 (1): 46–53.
#'
#' Trapnell, Cole, Brian a Williams, Geo Pertea, Ali Mortazavi, Gordon Kwan, Marijke J van Baren, Steven L Salzberg, Barbara J Wold, and Lior Pachter. 2010. “Transcript Assembly and Quantification by RNA-Seq Reveals Unannotated Transcripts and Isoform Switching During Cell Differentiation.” Nat. Biotechnol. 28 (5): 511–15.
#'
#' Wu, Hao, Chi Wang, and Zhijin Wu. 2013. “A new shrinkage estimator for dispersion improves differential expression detection in RNA-seq data.” Biostatistics 14 (2). Oxford University Press: 232–43. [https://doi.org/10.1093/biostatistics/kxs033](https://doi.org/10.1093/biostatistics/kxs033).
#'
#' Zeeberg, Barry R, Joseph Riss, David W Kane, Kimberly J Bussey, Edward Uchio, W Marston Linehan, J Carl Barrett, and John N Weinstein. 2004. “Mistaken Identifiers: Gene Name Errors Can Be Introduced Inadvertently When Using Excel in Bioinformatics.” BMC Bioinformatics 5: 80.
#'
#' Ziemann, Mark, Yotam Eren, and Assam El-Osta. 2016. “Gene Name Errors Are Widespread in the Scientific Literature.” Genome Biol. 17 (1): 1–3.
#'
#'