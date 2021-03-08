#' plot volcano for easy comparision
#'
#' @param DE list of Differentially Expressed Genes from Various methods.
#' @param FoldChange is cuoff to mark genes as Differentially expressed
#' @param cutoff is either pvalue of FDR cutoff for filtering
#' @param DE.Only Logical, show unchanged genes or not default TRUE
#' @param type.sig vector \code{c('p', 'FDR')} default \code{'p'}
#' @param show.genes logical weather to show genes as label
#' @param genes vector tor of genes for which only volcano plot will be drawn
#'
#' @importFrom ggplot2 ggplot
#' @importFrom ggrepel geom_text_repel
#' @return NULL
#'
#' @examples
#'
#' data("DEG")
#' DE.list<-list("edger" =dge_edger, "edgerql" = dge_edgerql,
#'                  "deseq2" = dge_deseq2, "voom" = dge_voom )
#' genes.list<-c("Gene7673", "Gene28034", "Gene38639")
#' multiVolcano(DE.list, FoldChange = 1.4, DE.Only = FALSE )
#' multiVolcano(DE.list, FoldChange = 1.4, DE.Only = FALSE , genes = genes.list)
#'
#'
#' @export


multiVolcano <- function(DE, FoldChange = 0 , cutoff=0.05 , DE.Only= TRUE ,
                          type.sig='p' , show.genes =FALSE , genes = NULL)
{

  stopifnot(is.logical(DE.Only) , is.numeric(FoldChange) , is.numeric(cutoff) ,
            is.logical(show.genes))


FC <- c("log2FoldChange", "logFC")
pv <-c( "pvalue", "P.Value", "PValue")

if (type.sig=="FDR"){
  pv <-c( "FDR", "padj", "adj.P.Val")
}
names <- names(DE)
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
if(is.null(genes))
{
DEG.df$genes[DEG.df$Regulation == "No"] <- NA
} else {
  DEG.df <- DEG.df[DEG.df$genes %in% genes,]
}

col.vec=c("blue", "grey", "red")
if(!show.genes)
{
  DEG.df$genes <- NA
}
if(DE.Only)
{
  DEG.df <- DEG.df[DEG.df$Regulation != "No",]
  col.vec=c("blue", "red")
}

ggplot2::ggplot(data=DEG.df, ggplot2::aes(x=FoldChange, y=-log10(pvalue),
                                          col=Regulation, label=genes)) +
  ggplot2::geom_point() +
  ggplot2::theme_minimal() +
  ggplot2::scale_color_manual(values=col.vec) +
  ggrepel::geom_text_repel() +
  ggplot2::facet_wrap(~Method)


}

