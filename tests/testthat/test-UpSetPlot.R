context("UpSetPlot")

test_that("mat", {
  data("DEG")
  DE.list<-list("edger" =dge_edger, "edgerql" = dge_edgerql,
                "deseq2" = dge_deseq2, "voom" = dge_voom )
  expect_error(UpSetPlot(DE=DE.list , FoldChange=10 , cutoff =0.01 , type.sig = "FDR"))
  
})
