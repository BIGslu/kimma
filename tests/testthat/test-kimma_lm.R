test_tomodel_path <- testthat::test_path("test_data", "to_model_ls.Rds")

to.model.ls <- readRDS(test_tomodel_path)
to.model.df <- to.model.ls[["to.model"]]


testthat::test_that("kimma_lm produces correct results with gene weights", {

  # check for 1st gene
  tst.df.gene1 <- to.model.df[to.model.df$rowname == "ENSG00000000460", ]

  res.gene1 <- kimma_lm(
    model.lm = "expression~virus+asthma",
    to.model.gene = tst.df.gene1,
    gene = "ENSG00000000460",
    use.weights = TRUE,
    metrics = FALSE
  )

  coef.tab <- res.gene1$fit$coefficients
  intercept.est <- coef.tab[names(coef.tab) == "(Intercept)"][[1]]
  virus.est <- coef.tab[names(coef.tab) == "virusHRV"][[1]]
  asthma.est <- coef.tab[names(coef.tab) == "asthmaasthma"][[1]]

  testthat::expect_equal(
    object = intercept.est,
    expected = 7.3243,
    tolerance = 0.001
  )

  testthat::expect_equal(
    object = virus.est,
    expected = 0.4630,
    tolerance = 0.001
  )

  testthat::expect_equal(
    object = asthma.est,
    expected = 0.7872,
    tolerance = 0.001
  )


  # check for 2nd gene
  tst.df.gene2 <- to.model.df[to.model.df$rowname == "ENSG00000001460", ]

  res.gene2 <- kimma_lm(
    model.lm = "expression~virus+asthma",
    to.model.gene = tst.df.gene2,
    gene = "ENSG00000001460",
    use.weights = TRUE,
    metrics = FALSE
  )

  coef.tab <- res.gene2$fit$coefficients
  intercept.est <- coef.tab[names(coef.tab) == "(Intercept)"][[1]]
  virus.est <- coef.tab[names(coef.tab) == "virusHRV"][[1]]
  asthma.est <- coef.tab[names(coef.tab) == "asthmaasthma"][[1]]

  testthat::expect_equal(intercept.est, 6.2902, tolerance = 0.001)
  testthat::expect_equal(virus.est, 0.9594, tolerance = 0.001)
  testthat::expect_equal(asthma.est, 0.7352, tolerance = 0.001)
})


testthat::test_that("kimma_lm produces correct results without gene weights", {

  tst.df <- to.model.df[to.model.df$rowname == "ENSG00000000460", ]

  res <- kimma_lm(
    model.lm = "expression~virus+asthma",
    to.model.gene = tst.df,
    gene = "ENSG00000000460",
    use.weights = FALSE,
    metrics = FALSE
  )

  coef.tab <- res$fit$coefficients
  intercept.est <- coef.tab[names(coef.tab) == "(Intercept)"][[1]]
  virus.est <- coef.tab[names(coef.tab) == "virusHRV"][[1]]
  asthma.est <- coef.tab[names(coef.tab) == "asthmaasthma"][[1]]

  testthat::expect_equal(intercept.est, 7.1501, tolerance = 0.001)
  testthat::expect_equal(virus.est, 0.5873, tolerance = 0.001)
  testthat::expect_equal(asthma.est, 0.9121, tolerance = 0.001)

})


testthat::test_that("kimma_lm adds fit metrics if set to TRUE", {

  tst.df <- to.model.df[to.model.df$rowname == "ENSG00000000460", ]

  res <- kimma_lm(
    model.lm = "expression~virus+asthma",
    to.model.gene = tst.df,
    gene = "ENSG00000000460",
    use.weights = TRUE,
    metrics = TRUE
  )

  # metrics data is added
  testthat::expect_true("metrics" %in% names(res))

  testthat::expect_equal(res$metrics[["sigma"]], 0.7763, tolerance = 0.001)
  testthat::expect_equal(res$metrics[["AIC"]], 28.9584, tolerance = 0.001)
  testthat::expect_equal(res$metrics[["BIC"]], 30.8980, tolerance = 0.001)
  testthat::expect_equal(res$metrics[["Rsq"]], 0.3854, tolerance = 0.001)
  testthat::expect_equal(res$metrics[["adj_Rsq"]], 0.2488, tolerance = 0.001)
})
