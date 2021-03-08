context("multiVolcano")

test_that("multiVolcano", {
  data("DEG")
  DE.list<-list("edger" =dge_edger, "edgerql" = dge_edgerql,
                "deseq2" = dge_deseq2, "voom" = dge_voom )
  expect_error(multiVolcano(DE.list, FoldChange = 1.4, DE.Only = "ABC" ))  
})