#' limma voom method
#'
#' @param count.table count table
#' @param design.matrix design matrix
#' @importFrom limma voom
#' @importFrom limma lmFit
#' @importFrom limma eBayes
#' @importFrom limma topTable
#' @return data.frame of differential gene expression quantification by Voom
#'
#' @examples
#'
#' data("count.table")
#' ss <- data.frame("trt" = colnames(count.table))
#' rownames(ss) <- ss$trt
#' ss$trt <- as.integer(grepl("T",ss[,1]))
#' design <- model.matrix(~ss$trt)
#' 
#' VoomDGE <- MethodVoom(count.table = count.table , design.matrix = design)
#' @export

MethodVoom <- function(count.table , design.matrix )
{

v <- limma::voom(count.table, design.matrix, plot=FALSE)
vfit <- limma::lmFit(v, design.matrix)
efit <- limma::eBayes(vfit)
dge <- limma::topTable(efit,n=Inf)
return(dge)
}
