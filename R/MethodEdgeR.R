#' EdgeR method
#'
#' @param count.table count table
#' @param design.matrix design matrix
#' @importFrom edgeR DGEList
#' @importFrom edgeR calcNormFactors
#' @importFrom edgeR estimateDisp
#' @importFrom edgeR glmFit
#' @importFrom edgeR topTags
#' @return data.frame of differntial gene expression quantification by EdgeR
#'
#' @examples
#'
#' data("count.table")
#' ss <- data.frame("trt" = colnames(count.table))
#' rownames(ss) <- ss$trt
#' ss$trt <- as.integer(grepl("T",ss[,1]))
#' 
#' design <- model.matrix(~ss$trt)
#' 
#' 
#' EdgeR <- MethodEdgeR(count.table = count.table, design.matrix = design)
#' @export


MethodEdgeR <- function(count.table , design.matrix)
{
  #stopifnot(is.null(setdiff(colnames(count.table), row.names(design.matrix)))))

z <- edgeR::DGEList(counts=count.table)
z <- edgeR::calcNormFactors(z)
z <- edgeR::estimateDisp(z, design.matrix,robust=TRUE,prior.df=1)
fit <- edgeR::glmFit(z, design.matrix)
lrt <- edgeR::glmLRT(fit)
dge<-as.data.frame(edgeR::topTags(lrt,n=Inf))
dge$dispersion<-lrt$dispersion
dge<-merge(dge,lrt$fitted.values,by='row.names')
rownames(dge)=dge$Row.names
dge$Row.names=NULL
return(dge)

}
