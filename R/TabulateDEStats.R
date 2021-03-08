#' plot Prapotional venn diagram for upregulated and downregulated genes
#'
#' @param DE list of Differentially Expressed Genes from Various methods.
#' @param FoldChange is cuoff to mark genes as Differentially expressed
#' @param cutoff is either pvalue of FDR cutoff for filtering
#' @param DE.Only Logical, show unchanged genes or not default TRUE
#' @param type.sig vector \code{c('p', 'FDR')} default \code{'p'}
#' @importFrom ggplot2 ggplot
#' @importFrom reshape2 melt
#' @importFrom gridExtra grid.arrange
#' @return data frame of Deferentially expressed genes counts in various methods
#'
#' @examples
#'
#' data("DEG")
#' DE.list<-list("edger" =dge_edger, "edgerql" = dge_edgerql,
#' "deseq2" = dge_deseq2, "voom" = dge_voom )
#' TabulateStats(DE=DE.list , FoldChange=1.5 , cutoff =0.01 , type.sig ="FDR")
#'
#'
#' @export

TabulateStats <- function(DE  , FoldChange = 0 , cutoff=0.05 , DE.Only= TRUE ,
                          type.sig='p')
    {

  stopifnot(is.logical(DE.Only) , is.numeric(FoldChange) , is.numeric(cutoff))

  FC <- c("log2FoldChange", "logFC")
  pv <-c( "pvalue", "P.Value", "PValue")

  if (type.sig=="FDR"){
     pv <-c( "FDR", "padj", "adj.P.Val")
  }
      names<-names(DE)
      tableStats<-list()
      for(name in names)
      {

        DERes<-DE[[name]]
        DERes$regulation<-"No"
        DERes$regulation[DERes[intersect(FC , colnames(DERes))] > FoldChange &
          DERes[intersect(pv , colnames(DERes))] < cutoff] <- "UP"

        DERes$regulation[DERes[intersect(FC , colnames(DERes))] < -1* FoldChange
                  & DERes[intersect(pv , colnames(DERes))] < cutoff] <- "DOWN"
        tableStats[[name]] <- as.matrix(table(DERes$regulation))
      }
  result<-as.data.frame(tableStats)
  result$id<-row.names(result)
  result.m<-reshape2::melt(result, id.vars = "id")
  if(DE.Only)
  {
    result.m <- result.m[result.m$id != "No",]
  }
  colnames(result.m) <- c("Condition" , "Methods" ,"Count")
  plt <- ggplot2::ggplot(result.m , mapping=ggplot2::aes(
    fill = Condition, x = Methods, y = Count)) +
    ggplot2::geom_bar(position = "dodge" , stat="identity")  +
    ggplot2::theme_classic() +
    ggplot2::ggtitle("Barplot #DE status in various methods")
  print(plt)
  result$id<-NULL
  return(result)
}

