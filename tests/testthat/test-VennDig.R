context("VennDig")

test_that("VennDig", {
  data("DEG")
  DE.list<-list("edger" =dge_edger, "edgerql" = dge_edgerql,
                "deseq2" = dge_deseq2, "voom" = dge_voom )
  expect_error(VennDig(DE=DE.list , FoldChange=1.5 , cutoff =p, type.sig ="FDR"))  
})