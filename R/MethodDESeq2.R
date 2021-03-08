#' DESeq2 method
#'
#' @param count.table count table
#' @param design.matrix design matrix
#' @importFrom DESeq2 DESeqDataSetFromMatrix
#' @importFrom DESeq2 DESeq
#' @importFrom DESeq2 results
#' @importFrom DESeq2 vst
#' @importFrom SummarizedExperiment assay
#' @return data.frame of differential gene expression quantification by Voom
#'
#' @examples
#'
#' data("count.table")
#' design <- data.frame("trt" = colnames(count.table))
#' rownames(design) <- design$trt
#' design$trt <- as.integer(grepl("T",design[,1]))
#' design$trt <- paste0("C",design$trt)
#'
#' DESeq2.DGE <- MethodDESeq2(count.table = count.table, design.matrix = design)
#' @export

MethodDESeq2 <- function(count.table , design.matrix )
{

  
  condition <- factor(design.matrix[,1])
dds <- DESeq2::DESeqDataSetFromMatrix(countData = count.table , colData = 
                                data.frame(condition), design = ~ condition )
res <- DESeq2::DESeq(dds)
z<- DESeq2::results(res)
vsd <- DESeq2::vst(dds, blind=FALSE)
zz<-cbind(as.data.frame(z),SummarizedExperiment::assay(vsd))
dge<-as.data.frame(zz[order(zz$pvalue),])
return(dge)
}
