context("TabulateStats")

test_that("TabulateStats", {
  data("DEG")
  DE.list<-list("edger" =dge_edger, "edgerql" = dge_edgerql,
                "deseq2" = dge_deseq2, "voom" = dge_voom )
  expect_error(TabulateStats(DE=DE.list , FoldChange=1.5 , cutoff =0.01 , DE.Only = "ABC"))  
})