#' ---
#' title: "Biological QA (mod'ed from Charlotte Soneson's)"
#' author: "Nicolas Delhomme, mod'ed from Charlotte Soneson's"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#' ---
#' # Setup
#' 
#' Here we load the needed libraries we will need for the analysis and set some graphical parameters.
#' * Libraries
suppressPackageStartupMessages({
  library(DESeq2)
  library(edgeR)
  library(here)
  library(pheatmap)
  library(RColorBrewer)
  library(VennDiagram)
})

#' * Helpers
source(here("R/volcanoPlot.R"))

#' * Graphics
pal <- brewer.pal(8,"Dark2")  

#' * Load the data
load(here("data/DE/differentialExpressionObjects.rda"))

#' # Differential expression analysis
#' ## Performing differential expression testing with DESeq2
#' 
#' As we have already specified an experimental design when we created the _DESeqDataSet_, 
#' we can run the differential expression pipeline on the raw counts with a single call to the function `DESeq`.
#'We can also plot the estimated dispersions.

#' Very easy, just one command
ds_se <- DESeq2::DESeq(ds_se)

#' Here, we plot the dispersion estimation to ensure that the assumption that most genes are not differentially expressed holds
DESeq2::plotDispEsts(ds_se)

#' This function will print out a message for the various steps it performs. These are described in more detail in the manual page
#' for `DESeq`, which can be accessed by typing `?DESeq`. Briefly these are: the estimation of size factors (controlling for 
#' differences in the sequencing depth of the samples), the estimation of dispersion values for each gene, and fitting a generalized linear model.
#' 
#' A DESeqDataSet is returned that contains all the fitted parameters within it, and the following section describes 
#' how to extract out results tables of interest from this object.
#' 
#' But wait! Before we proceed, what was our design again?
design(ds_se)
colData(ds_se)[,c("Tissue","Time")]

#' We have multiple possible contrast (_i.e_ possible comparisons), which is the one R returns by default?
resultsNames(ds_se)

#' By default the last one from that list (which is sorted alphabetically). 
#' Here we have the Tissue contribution of _ZE_, the Time contribution of _B8_ and
#' as last the contribution specific to the interaction _ZE_:_B8_
#' 
#' Hence, running the command to retrieve the results will by default return the latter. Let's rather start by looking 
#' at the Tissue ZE _vs._ FMG
res <- results(ds_se,name="Tissue_ZE_vs_FMG")

#' OK, but what was the baseline? FMG and B4
levels(ds_se$Time)
levels(ds_se$Tissue)

#' So we will be looking at the effect of the tissue during the developmental time B4
head(res)

#' As res is a DataFrame object, it carries metadata with information on the meaning of the columns:
mcols(res, use.names = TRUE)

#' The first column, baseMean, is a just the average of the normalized count values, dividing by size factors, 
#' taken over all samples in the DESeqDataSet. The remaining four columns refer to a specific contrast, 
#' namely the comparison of the ZE level over the FMG level for the factor variable Tissue, at the developmental time B4 (the baseline).
#' 
#' The column log2FoldChange is the effect size estimate. It tells us how much the gene’s expression seems to have changed in between the Tissues during the Time B4. 
#' This value is reported on a logarithmic scale to base 2: for example, a log2 fold change of 1.5 means that the gene’s expression is increased by a multiplicative factor of \(2^{1.5} \approx 2.82\).
#' 
#' Of course, this estimate has an uncertainty associated with it, which is available in the column lfcSE, the standard error estimate for the log2 fold change estimate. 
#' We can also express the uncertainty of a particular effect size estimate as the result of a statistical test. The purpose of a test for differential expression is to 
#' test whether the data provides sufficient evidence to conclude that this value is really different from zero. DESeq2 performs for each gene a hypothesis test to see 
#' whether evidence is sufficient to decide against the null hypothesis that there is zero effect of the tissue type on the gene and that the observed difference between tissues
#'  was merely caused by experimental variability (i.e., the type of variability that you can expect between different samples in the same tissue group). 
#'  As usual in statistics, the result of this test is reported as a p value, and it is found in the column pvalue. Remember that a p value indicates the probability that 
#'  an effect as strong as the observed one, or even stronger, would be seen under the situation described by the null hypothesis.
#'  
#'  We can also summarize the results with the following line of code, which reports some additional information, that will be covered in later sections.
summary(res)
hist(res$pvalue,breaks=seq(0,1,.01))

#' Note that there are few genes with differential expression between Tissues at the FDR level of 10%. 
#' Nonetheless, there are two ways to be more strict about which set of genes are considered significant:
#' 
#' * lower the false discovery rate threshold (the threshold on padj in the results table)
#' * raise the log2 fold change threshold from 0 using the lfcThreshold argument of results

#' If we lower the false discovery rate threshold, we should also tell this value to `results()`, 
#' so that the function will use an alternative threshold for the optimal independent filtering step:
res.05 <- results(ds_se, name="Tissue_ZE_vs_FMG", alpha = 0.05)
table(res.05$padj < 0.05)

#' If we want to raise the log2 fold change threshold, so that we test for genes that show more substantial 
#' changes between tissues, we simply supply a value on the log2 scale. For example, by specifying lfcThreshold = 1, 
#' we test for genes that show significant effects between tissues on gene counts more than doubling or less than halving, because `2^1 = 2`.
resLFC1 <- results(ds_se, name="Tissue_ZE_vs_FMG", lfcThreshold = 1)
summary(resLFC1)
table(resLFC1$padj < 0.1)

#' Following Schurch _et al._, (RNA, 2016) recommandations, we would have for 3 replicates per conditions:
resSchurch <- results(ds_se, name="Tissue_ZE_vs_FMG", lfcThreshold = 0.5, alpha = 0.01)
summary(resSchurch)

#' Sometimes a subset of the p values in res will be NA (“not available”). This is DESeq’s way of reporting 
#' that all counts for this gene were zero, and hence no test was applied. In addition, p values can be assigned `NA`
#' if the gene was excluded from analysis because it contained an extreme count outlier. 
#' For more information, see the outlier detection section of the `DESeq2` vignette.
#' 
#' ### Visualisation
#' #### Assessment
#' A volcano plot (the naming is probably obvious) is very useful to check
#' whether the assumption that most genes are not DE holds.
#' 
#' The plot has a representation of density as a color gradient 
#' from sparse to dense (gray -> blue -> red -> yellow). Note that the light blue, 
#' confusingly shows DE genes.
#' 
#' Here clearly most genes are located at 0,0 (the yellow dot)
volcanoPlot(resSchurch)

#' #### Individual genes
#' With DESeq2, there is also an easy way to plot the (normalized, transformed) 
#' counts for specific genes, using the plotCounts function:
plotCounts(ds_se, gene = "G10D_2", intgroup = "Tissue", 
           normalized = TRUE, transform = FALSE)

#' This is actually an insteresting case as G10D_2 is a Long Terminal Repeat (LTR). 
#' Fascinating to see that they are deregulated in the ZE, as compared to the FMG, 
#' probably related to the modified epigenome of the embryo.
#' 
#' ## Performing differential expression testing with edgeR
#' 
#' Next we will show how to perform differential expression analysis with edgeR. Recall that we have a _DGEList_ `dge`, containing all the necessary information:

names(dge)

#' We first define a design matrix, using the same formula syntax as for _DESeq2_ above.

design <- model.matrix(~ Tissue * Time, data = dge$samples)

#' While `DESeq2` performs independent filtering of lowly expressed genes internally, this is done by the user before applying edgeR. 
#' Here, we filter out lowly expressed genes using the filterByExpr() function, and then estimate the dispersion for each gene. 
#' Note that it is important that we specify the design in the dispersion calculation. Afterwards, we plot the estimated dispersions.
#' 
#' First, we perform the independent filtering (based on the model, filter the genes that have no 
#' power to be detected as differentially expressed. This aims at removing genes so lowly expressed 
#' that they represent "noise")
keep <- edgeR::filterByExpr(dge, design)
table(keep)
dge <- dge[keep, ]

#' Then we estimate the dispersion
dge <- edgeR::estimateDisp(dge, design)

#' And we can take a look at the dispersion estimation (dispersion == biological coefficient of variation)
edgeR::plotBCV(dge)

#' Finally, we fit the generalized linear model and perform the test. In the `glmQLFTest` function, we indicate which coefficient (which column in the design matrix) that we would like to test for. 
#' It is possible to test more general contrasts as well, and the user guide contains many examples on how to do this. The topTags function extracts the top-ranked genes. 
#' You can indicate the adjusted p-value cutoff, and/or the number of genes to keep.
fit <- edgeR::glmQLFit(dge, design)

#' What is the design?
design

# To test as we did for DESeq2, we select the second column
qlf <- edgeR::glmQLFTest(fit, coef = 2)

#' Finally looking at the results
tt.all <- edgeR::topTags(qlf, n = nrow(dge), sort.by = "none") # all genes
hist(tt.all$table$PValue,breaks=seq(0,1,.01))
tt <- edgeR::topTags(qlf, n = nrow(dge), p.value = 0.01) # genes with adj.p<0.01
tt10 <- edgeR::topTags(qlf) # just the top 10 by default
tt10

#' The columns in the edgeR result data frame are similar to the ones output by DESeq2. edgeR 
#' represents the overall expression level on the log-CPM scale rather than on the normalized 
#' count scale that DESeq2 uses. The F column contains the test statistic, and the FDR column 
#' contains the Benjamini-Hochberg adjusted p-values.
#' 
#' We can compare the sets of significantly differentially expressed genes to see how the results from the two packages overlap:
shared <- intersect(rownames(res), tt.all$table$gene.id)
table(DESeq2 = res$padj[match(shared, rownames(res))] < 0.1, 
      edgeR = tt.all$table$FDR[match(shared, tt.all$table$gene.id)] < 0.1)

#' We can do the same visually
grid.newpage()
d.sel <- match(shared, rownames(res))
e.sel <- match(shared, tt.all$table$gene.id)
grid.draw(venn.diagram(list(DESeq2 = rownames(res)[d.sel][res$padj[d.sel] < 0.1 & !is.na(res$padj[d.sel])], 
                            edgeR = tt.all$table$gene.id[e.sel][tt.all$table$FDR[e.sel] < 0.1 & !is.na(tt.all$table$FDR[e.sel])]),
                       NULL,fill=pal[1:2]))

#' We can also compare the two result lists by the ranks:
plot(rank(res$pvalue[match(shared, rownames(res))]), 
     rank(tt.all$table$PValue[match(shared, tt.all$table$gene.id)]), 
     cex = 0.1, xlab = "DESeq2", ylab = "edgeR")

#' Also with edgeR we can test for significance relative to a fold-change threshold, 
#' using the function glmTreat. Below we set the log fold-change threshold to 1 (i.e., fold change threshold equal to 2), as for DESeq2 above.
treatres <- edgeR::glmTreat(fit, coef = ncol(design), lfc = 1)
tt.treat <- edgeR::topTags(treatres, n = nrow(dge), sort.by = "none")

#' ## Multiple testing
#' In high-throughput biology, we are careful to not use the p values directly as evidence against the null, but to correct for multiple testing. 
#' What would happen if we were to simply threshold the p values at a low value, say 0.05?
sum(res$pvalue < 0.05, na.rm = TRUE)

#' and that many genes not NA
sum(!is.na(res$pvalue))

#' Now, assume for a moment that the null hypothesis is true for all genes, i.e., no gene is affected in either Tissue. 
#' Suppose we are interesting in a significance level of 0.05. Then, by the definition of the p value, we expect up to 5% of the 
#' genes to have a p value below 0.05. This amounts to
round(sum(!is.na(res$pvalue)) * 0.05)

#' If we just considered the list of genes with a p value below 0.05 as differentially expressed, this list should therefore be 
#' expected to contain many false positives:
round(sum(!is.na(res$pvalue))*0.05 / sum(res$pvalue < 0.05, na.rm = TRUE), 2)

#' DESeq2 and edgeR use the Benjamini-Hochberg (BH) adjustment (Benjamini and Hochberg 1995) as implemented in the base R p.adjust 
#' function; in brief, this method calculates for each gene an adjusted p value that answers the following question: if one called 
#' significant all genes with an adjusted p value less than or equal to this gene’s adjusted p value threshold, what would be the
#' fraction of false positives (the false discovery rate, FDR) among them, in the sense of the calculation outlined above? These 
#' values, called the BH-adjusted p values, are given in the column padj of the res object from DESeq2, and in the FDR column in 
#' the TopTags object from edgeR.
#'  
#' The FDR is a useful statistic for many high-throughput experiments, as we are often interested in reporting or focusing on a
#' set of interesting genes, and we would like to put an upper bound on the percent of false positives in this set.
#'   
#' Hence, if we consider a fraction of 10% false positives acceptable, we can consider all genes with an adjusted p value 
#' below 10% = 0.1 as significant. How many such genes are there?
sum(res$padj < 0.1, na.rm = TRUE)

#' # Plotting results
#' ## MA plot with DESeq2
#' An MA-plot (Dudoit et al. 2002) provides a useful overview for an experiment with a two-group comparison (Figure below). 
#' The log2 fold change for a particular comparison is plotted on the y-axis and the average of the counts normalized by size 
#' factor is shown on the x-axis (“M” for minus, because a log ratio is equal to log minus log, and “A” for average). Each gene 
#' is represented with a dot. Genes with an adjusted p value below a threshold (here 0.1, the default with DESeq2) are shown in blue.
DESeq2::plotMA(res, ylim = c(-5, 5))

#' ## MA / Smear plot with edgeR
#' In edgeR, the MA plot is obtained via the plotSmear function. The genes below the adjusted p value threshold are shown in red.
edgeR::plotSmear(qlf, de.tags = tt$table$gene.id)

#' ## Heatmap of the most significant genes
#' Another way of representing the results of a differential expression analysis is to construct a heatmap of the top differentially 
#' expressed genes. Here, we would expect the contrasted sample groups to cluster separately. A heatmap is a “color coded expression matrix”, 
#' where the rows and columns are clustered using hierarchical clustering. Typically, it should not be applied to counts, but works better 
#' with transformed values. Here we show how it can be applied to the variance-stabilized values generated above. We choose the top 30 differentially expressed genes.
#'  There are many functions in R that can generate heatmaps, here we show the one from the pheatmap package.
mat <- assay(vsd)[head(order(res$padj), 30), ]
mat <- mat - rowMeans(mat)
df <- as.data.frame(colData(vsd)[, c("Tissue", "Time")])
pheatmap(mat, annotation_col = df)

#' # Session Info
#'  ```{r, session info, echo=FALSE}
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


