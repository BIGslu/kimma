test_tomodel_path <- testthat::test_path("test_data", "to_model_ls.Rds")

to.model.ls <- readRDS(test_tomodel_path)
to.model.df <- to.model.ls[["to.model"]]

tst.df <- to.model.df[to.model.df$rowname == "ENSG00000000460", ]

res.lm.fit <- kimma_lm(
  model_lm = "expression~virus+asthma",
  to_model_gene = tst.df,
  gene = "ENSG00000000460",
  use_weights = TRUE,
  metrics = FALSE
)

testthat::test_that("kmFit_contrast produces correct results", {

  res <- kmFit_contrast(
    fit = res.lm.fit[["fit"]],
    contrast_var = "virus",
    to_model_gene = tst.df,
    genotype_name = NULL
  )

  estimate.val <- res[["estimate"]]
  contrast.ref <- res[["contrast_ref"]]
  contrast.lvl <- res[["contrast_lvl"]]

  testthat::expect_true(contrast.ref == "none")
  testthat::expect_true(contrast.lvl == "HRV")
  testthat::expect_equal(estimate.val, 0.4630, tolerance = 0.001)
})


testthat::test_that("kmFit_contrast produces correct results with interaction term", {

  res <- kmFit_contrast(
    fit = res.lm.fit[["fit"]],
    contrast_var = "virus:asthma",
    to_model_gene = tst.df,
    genotype_name = NULL
  )

  testthat::expect_equal(
    res[res$contrast_ref == "none healthy" & res$contrast_lvl == "HRV healthy", ][["estimate"]],
    0.4630,
    tolerance = 0.001
  )

  testthat::expect_equal(
    res[res$contrast_ref == "none healthy" & res$contrast_lvl == "none asthma", ][["estimate"]],
    0.7872,
    tolerance = 0.001
  )

  testthat::expect_equal(
    res[res$contrast_ref == "none healthy" & res$contrast_lvl == "HRV asthma", ][["estimate"]],
    1.2502,
    tolerance = 0.001
  )

  testthat::expect_equal(
    res[res$contrast_ref == "HRV healthy" & res$contrast_lvl == "none asthma", ][["estimate"]],
    0.3241,
    tolerance = 0.001
  )

  testthat::expect_equal(
    res[res$contrast_ref == "HRV healthy" & res$contrast_lvl == "HRV asthma", ][["estimate"]],
    0.7872,
    tolerance = 0.001
  )

  testthat::expect_equal(
    res[res$contrast_ref == "none asthma" & res$contrast_lvl == "HRV asthma", ][["estimate"]],
    0.4630,
    tolerance = 0.001
  )

})

testthat::test_that(
  "kmFit_contrast produces correct results when interaction term is numeric", {

    res.lm.fit <- kimma_lm(
      model_lm = "expression~virus + virus * median_cv_coverage",
      to_model_gene = tst.df,
      gene = "ENSG00000000460",
      use_weights = TRUE,
      metrics = FALSE
    )

    res <- kmFit_contrast(
      fit = res.lm.fit[["fit"]],
      contrast_var = "virus:median_cv_coverage",
      to_model_gene = tst.df,
      genotype_name = NULL
    )

    estimate.val <- res[["estimate"]]
    contr.var <- res[["variable"]]
    contrast.ref <- res[["contrast_ref"]]
    contrast.lvl <- res[["contrast_lvl"]]

    testthat::expect_true(contr.var == "virus*median_cv_coverage")
    testthat::expect_true(contrast.ref == "none")
    testthat::expect_true(contrast.lvl == "HRV")
    testthat::expect_equal(estimate.val, 3.5497, tolerance = 0.001)

  }
)
