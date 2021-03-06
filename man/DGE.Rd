% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DGE.R
\name{DGE}
\alias{DGE}
\title{limma voom method}
\usage{
DGE(
  count.table,
  design.matrix,
  method = c("EdgeR", "EdgeRLTR", "voom", "DESeq2")
)
}
\arguments{
\item{count.table}{count table}

\item{design.matrix}{design matrix}

\item{method}{DGEs methods EdgeR , EdgeRLTR, voom , DESeq2. User can choose
any of the listed methods or all methods will be run by default.}
}
\value{
list of data frame of differential gene expression quantification by
various method
}
\description{
limma voom method
}
\examples{

data("count.table")
design <- data.frame("trt" = colnames(count.table))
rownames(design) <- design$trt
design$trt <- as.integer(grepl("T",design[,1]))

DGE.list <- DGE(count.table = count.table , design.matrix = design)
}
