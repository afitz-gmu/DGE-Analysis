#' plot Venn diagram for upregulated and downregulated genes
#'
#' @param DE list of Differentially Expressed Genes from Various methods.
#' @param FoldChange is cuoff to mark genes as Differentially expressed
#' @param cutoff is either pvalue of FDR cutoff for filtering
#' @param type.sig vector \code{c('p', 'FDR')} default \code{'p'}
#' @importFrom ggvenn ggvenn
#' @importFrom ggplot2 ggtitle
#' @importFrom gridExtra grid.arrange
#' @return \code{NULL}
#'
#' @examples
#'
#' data("DEG")
#' DE.list<-list("edger" =dge_edger, "edgerql" = dge_edgerql,
#' "deseq2" = dge_deseq2, "voom" = dge_voom )
#' VennDig(DE=DE.list , FoldChange=1.5 , cutoff =0.01 , type.sig ="FDR")
#'
#'
#' @export



VennDig <- function(DE, FoldChange = 0, cutoff=0.05, type.sig='p')
{

  stopifnot(is.list(DE), is.numeric(FoldChange) , is.numeric(cutoff))

  FC <- c("log2FoldChange", "logFC")
  pv <-c( "pvalue", "P.Value", "PValue")
  if (type.sig=="FDR"){
    pv <-c( "FDR", "padj", "adj.P.Val")
  }
  names<-names(DE)
  up<-list()
  down<-list()
  for(name in names)
  {

    DERes<-DE[[name]]
    DERes$regulation<-"No"
    DERes$regulation[DERes[intersect(FC , colnames(DERes))] > FoldChange &
                      DERes[intersect(pv , colnames(DERes))] < cutoff] <- "UP"
    up[[name]]<-row.names(DERes[DERes$regulation == "UP",])
    DERes$regulation[DERes[intersect(FC , colnames(DERes))] < -1* FoldChange
                  & DERes[intersect(pv , colnames(DERes))] < cutoff] <- "DOWN"
    down[[name]]<-row.names(DERes[DERes$regulation == "DOWN",])

  }
  gridExtra::grid.arrange( ggvenn::ggvenn(up, set_name_size=4 )+
                ggplot2::ggtitle("Upregulated Genes") ,
                ggvenn::ggvenn(down, set_name_size=4) +
                ggplot2::ggtitle("Downregulated Genes"),
                nrow=1)
}
