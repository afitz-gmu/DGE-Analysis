#' limma voom method
#'
#' @param count.table count table
#' @param design.matrix design matrix
#' @param method DGEs methods EdgeR , EdgeRLTR, voom , DESeq2. User can choose
#' any of the listed methods or all methods will be run by default.
#' @return list of data frame of differential gene expression quantification by
#' various method
#'
#' @examples
#'
#' data("count.table")
#' design <- data.frame("trt" = colnames(count.table))
#' rownames(design) <- design$trt
#' design$trt <- as.integer(grepl("T",design[,1]))
#'
#' DGE.list <- DGE(count.table = count.table , design.matrix = design)
#' @export

DGE <- function(count.table , design.matrix  , method = c("EdgeR", "EdgeRLTR",
                                                          "voom", "DESeq2"))
{


  DEG.list <- list()
  for(m in method)
  {
    if(m == "EdgeR" )
    {
      DEG.list[[m]]<-MethodEdgeR(count.table, design.matrix =
                                    model.matrix(~design.matrix[,1]))
    }
    if(m == "EdgeRLTR" )
    {
      DEG.list[[m]]<-MethodEdgeRLRT(count.table, design.matrix =
                                   model.matrix(~design.matrix[,1]))
    }
    if(m == "voom" )
    {
      DEG.list[[m]]<-MethodVoom(count.table, design.matrix =
                                      model.matrix(~design.matrix[,1]))
    }
    if(m == "DESeq2" )
    {
      DEG.list[[m]]<-MethodDESeq2(count.table, design.matrix = design.matrix)
    }


  }
  return(DEG.list)
}
