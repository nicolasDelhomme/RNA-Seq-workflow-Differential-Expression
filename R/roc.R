#' ---
#' title: "Zygotic Embryogenesis Seidr ROC"
#' author: "Nicolas Delhomme"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#' ---
#' # Setup
#' * Libraries
suppressPackageStartupMessages({
  library(here)
  library(pracma)
})

#' * function
"plotROC" <- function(f){
  
  dat <- read.delim(f,header=FALSE,skip=1,
                    col.names = c("TP","FP","PR"))
  
  auc <- round(trapz(dat[,2],dat[,1]),digits=3)
  
  plot(dat[,2],dat[,1],type="l",main=sprintf("%s (AUC = %s)",sub("\\.roc","",basename(f)),auc),
       xlab="False Positive Rate",ylab="True Positive Rate")
  
  abline(0,1,lty=2)

}

#' # Plot
sapply(dir(here("data/seidr/roc"),pattern="*\\.roc.gz",full.names = TRUE),plotROC)

#' # Session Info
#'  ```{r session info, echo=FALSE}
#'  sessionInfo()
#'  ```

