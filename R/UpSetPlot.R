#' plot UpSetPlot diagram for upregulated and downregulated genes
#'
#' @param DE list of Differentially Expressed Genes from Various methods.
#' @param FoldChange is cuoff to mark genes as Differentially expressed
#' @param cutoff is either pvalue of FDR cutoff for filtering
#' @param type.sig vector \code{c('p', 'FDR')} default \code{'p'}
#' @param regulation 
#' @importFrom UpSetR upset
#' @importFrom UpSetR fromList
#' @importFrom grid grid.text
#' @importFrom grid grid.edit
#' @importFrom grid grid.grab
#' @importFrom grid gpar
#' @importFrom gridExtra grid.arrange
#' @return \code{NULL}
#' @examples
#'
#' data("DEG")
#' DE.list<-list("edger" =dge_edger, "edgerql" = dge_edgerql,
#' "deseq2" = dge_deseq2, "voom" = dge_voom )
#' UpSetPlot(DE=DE.list , FoldChange=0 , cutoff =0.01 , type.sig = "FDR")
#'
#'
#' @export

UpSetPlot <- function(DE  , FoldChange = 0 , cutoff=0.05, type.sig=NULL)
{
  type.sig <-type.sig
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

  plot.ls <- list(
    UpRegulated = UpSetR::upset(UpSetR::fromList(up), order.by = "freq"),
    DownRegulated = UpSetR::upset(UpSetR::fromList(down), order.by = "freq"))

  for(v in names(plot.ls))
  {
    print(plot.ls[[v]])
    grid::grid.text(v ,x=0.65, y=0.97 ,gp =grid::gpar(fontsize=8))
    grid::grid.edit('arrange', name=v)
    vp<-grid::grid.grab()
    plot.ls[[v]] <-vp
  }
  
  gridExtra::grid.arrange(grobs=plot.ls , ncol=2)

}
