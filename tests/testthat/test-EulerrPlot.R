context("EulerrPlot")

test_that("EulerrPlot", {
  data("DEG")
  DE.list<-list("edger" =dge_edger, "edgerql" = dge_edgerql,
                "deseq2" = dge_deseq2, "voom" = dge_voom )
  expect_error(EulerrPlot(DE=DE.list , FoldChange=20 , cutoff ="l" , type.sig ="p"))  
})