test_tomodel_path <- testthat::test_path("test_data", "to_model_ls.Rds")

to.model.ls <- readRDS(test_tomodel_path)
to.model.df <- to.model.ls[["to.model"]]

tst.df <- to.model.df[to.model.df$rowname == "ENSG00000000460", ]

res.lm.fit <- kimma_lm(
  model.lm = "expression~virus+asthma",
  to.model.gene = tst.df,
  gene = "ENSG00000000460",
  use.weights = TRUE,
  metrics = FALSE
)

testthat::test_that("kmFit_contrast produces correct results", {

  res <- kmFit_contrast(
    fit = res.lm.fit[["fit"]],
    contrast.var = "virus",
    to.model.gene = tst.df,
    genotype.name = NULL
  )

  estimate.val <- res[["estimate"]]
  contrast.ref <- res[["contrast_ref"]]
  contrast.lvl <- res[["contrast_lvl"]]

  testthat::expect_true(contrast.ref == "none")
  testthat::expect_true(contrast.lvl == "HRV")
  testthat::expect_equal(estimate.val, 0.4630, tolerance = 0.001)
})
