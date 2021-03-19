#' plot and get summary
#'
#' @param DE list of Differentially Expressed Genes from Various methods.
#' @param FoldChange is cuoff to mark genes as Differentially expressed
#' @param cutoff is either pvalue of FDR cutoff for filtering
#' @param type.sig vector \code{c('p', 'FDR')} default \code{'p'}
#' @importFrom ggplot2 ggplot
#' @importFrom reshape2 acast
#' @return data frame of Deferentially expressed genes counts in various methods
#'
#' @examples
#'
#' data("DEG")
#' DE.list<-list("edger" =dge_edger, "edgerql" = dge_edgerql,
#' "voom" = dge_voom )
#' DGE.CI(DE=DE.list , FoldChange=1.2 , cutoff =0.01 , type.sig ="FDR")
#'
#'
#' @export(DGE.CI.Plot)
#'
DGE.CI.Plot <- function(DE  , FoldChange = 0 , cutoff=0.05 , type.sig='p') {

  stopifnot(is.numeric(FoldChange) , is.numeric(cutoff))

  FC <- c("log2FoldChange", "logFC")
  pv <-c( "pvalue", "P.Value", "PValue")

  if (type.sig=="FDR"){
    pv <-c( "FDR", "padj", "adj.P.Val")
  }
  names<-names(DE)
  DEG.df <- data.frame(FoldChange = numeric(), pvalue = numeric(),
                       Regulation = character() , Method = character() ,
                       genes = character())
  for(name in names)
  {

    DERes<-DE[[name]]
    DERes$Method <- name
    DERes$genes <- row.names(DERes)
    DERes$regulation<-"No"
    DERes$regulation[DERes[intersect(FC , colnames(DERes))] > FoldChange &
                       DERes[intersect(pv , colnames(DERes))] < cutoff] <- "UP"

    DERes$regulation[DERes[intersect(FC , colnames(DERes))] < -1* FoldChange
                     & DERes[intersect(pv , colnames(DERes))] < cutoff] <- "DOWN"
    df<-DERes[c(intersect(FC , colnames(DERes)), intersect(pv , colnames(DERes)),
                "regulation" , "Method" ,"genes")]
    colnames(df) <- colnames(DEG.df)
    DEG.df <-rbind(DEG.df , df)

  }
  gn<-unique(DEG.df[DEG.df$Regulation != "No",]$genes)
  df<-DEG.df[DEG.df$genes %in% gn ,]
  df$FoldChange=exp(df$FoldChange)
  
  
  fc <- reshape2::acast(data=df ,genes ~Method, value.var="FoldChange")
  colnames(fc) <-paste(colnames(fc),"FoldChange")
  
  ggplot2::ggplot(df , ggplot2::aes(x= as.numeric(as.factor(genes)) , y= -log10(pvalue) )) +
  ggplot2::geom_point(ggplot2::aes(pch = Regulation  , color=Method , size=FoldChange) ) +
  ggplot2::geom_smooth( method=loess , se=TRUE, color="red")  +
  ggplot2::scale_x_discrete(labels=gn , name="genes") +
  ggplot2::theme_minimal() +
  ggplot2::ggtitle("95 confidence interval") +
  ggplot2::geom_hline(yintercept=-log10(cutoff), col="cyan")
  
}
