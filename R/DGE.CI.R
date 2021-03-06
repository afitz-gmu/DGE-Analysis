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
#' @export
#'
DGE.CI <- function(DE  , FoldChange = 0 , cutoff=0.05 , type.sig='p') {

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
  
  Range= t(apply(fc,1,range))
  colnames(Range)<-c("Min","Max")
  fc<-cbind(fc, Range)
  pv <- reshape2::acast(data=df ,genes ~Method, value.var="pvalue")
  colnames(pv) <-paste(colnames(pv),"pvalue")
  Range= t(apply(pv,1,range))
  colnames(Range)<-c("Min","Max")
  pv<-cbind(pv, Range)
  res<-cbind(fc,pv)
  res <- as.data.frame(res)
  return(res)

}
