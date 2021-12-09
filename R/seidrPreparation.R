#' ---
#' title: "ZE Network data preparation"
#' author: "Nicolas Delhomme"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#'    code_folding: hide
#' ---
#' # Setup
suppressPackageStartupMessages({
  library(here)
  library(readr)
})

#' Helper
source(here("R/featureSelection.R"))

#' # Data
#' ```{r CHANGEME2,eval=FALSE,echo=FALSE}
#'  CHANGEME is the variance stabilised data, where the transformation was done taking the
#'  model into account (_i.e._ `blind=FALSE`)
#' ```
load(here("data/DE/differentialExpressionObjects.rda"))

#' # Filter
conds <- factor(paste(ds_se$Time,ds_se$Tissue))
sels <- rangeFeatureSelect(counts=as.matrix(assay(vsd)),
                           conditions=conds,
                           nrep=2)

#' ```{r CHANGEME3,eval=FALSE,echo=FALSE}
#'  CHANGEME is the vst cutoff devised from the plot above. The goal is to remove / reduce
#'  the signal to noise. Typically, this means trimming the data after the first sharp decrease 
#'  on the y axis, most visible in the non logarithmic version of the plot
#' ```
vst.cutoff <- 1

#' # Export
dir.create(here("seidr/data"),showWarnings=FALSE,recursive=TRUE)

#' * gene by column, without names matrix
write.table(t(assay(vsd)[sels[[vst.cutoff+1]],]),
            file=here("seidr/data/headless.tsv"),
            col.names=FALSE,
            row.names=FALSE,
            sep="\t",quote=FALSE)

#' * gene names, one row
write.table(t(sub("\\.[0-9]+$","",rownames(vsd)[sels[[vst.cutoff+1]]])),
            file=here("seidr/data/genes.tsv"),
            col.names=FALSE,
            row.names=FALSE,
            sep="\t",quote=FALSE)

#' # Session Info
#' ```{r session info, echo=FALSE}
#' sessionInfo()
#' ```
