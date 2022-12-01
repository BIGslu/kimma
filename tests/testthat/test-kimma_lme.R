test_tomodel_path <- testthat::test_path("test_data", "to_model_ls.Rds")

to.model.ls <- readRDS(test_tomodel_path)
to.model.test.df <- to.model.ls[["to.model"]]


testthat::test_that("kimma_lme produces correct results with gene weights", {

  tst.df <- to.model.test.df[to.model.test.df$rowname == "ENSG00000000460", ]

  res <- kimma_lme(
    model.lm = "expression~virus+asthma+(1|ptID)",
    to.model.gene = tst.df,
    gene = "ENSG00000000460",
    use.weights = TRUE,
    metrics = FALSE
  )

  res.fit <- as.data.frame(
    summary(res$fit)[["coefficients"]]
  )

  intercept.est <- res.fit[rownames(res.fit) == "(Intercept)", ][["Estimate"]]
  virus.est <- res.fit[rownames(res.fit) == "virusHRV", ][["Estimate"]]
  asthma.est <- res.fit[rownames(res.fit) == "asthmaasthma", ][["Estimate"]]

  testthat::expect_equal(intercept.est, 7.2840, tolerance = 0.001)
  testthat::expect_equal(virus.est, 0.4744, tolerance = 0.001)
  testthat::expect_equal(asthma.est, 0.8369, tolerance = 0.001)

})


testthat::test_that("kimma_lme produces correct results without gene weights", {

  tst.df <- to.model.test.df[to.model.test.df$rowname == "ENSG00000000460", ]

  res <- kimma_lme(
    model.lm = "expression~virus+asthma+(1|ptID)",
    to.model.gene = tst.df,
    gene = "ENSG00000001460",
    use.weights = FALSE,
    metrics = FALSE
  )

  res.fit <- as.data.frame(
    summary(res$fit)[["coefficients"]]
  )

  intercept.est <- res.fit[rownames(res.fit) == "(Intercept)", ][["Estimate"]]
  virus.est <- res.fit[rownames(res.fit) == "virusHRV", ][["Estimate"]]
  asthma.est <- res.fit[rownames(res.fit) == "asthmaasthma", ][["Estimate"]]

  testthat::expect_equal(intercept.est, 7.1501, tolerance = 0.001)
  testthat::expect_equal(virus.est, 0.5873, tolerance = 0.001)
  testthat::expect_equal(asthma.est, 0.9121, tolerance = 0.001)

})


testthat::test_that("kimma_lme fails if no random effects terms specified in model", {

  testthat::expect_error(
    kimma_lme(
      model.lm = "expression~virus",
      to.model.gene = tst.df,
      gene = "ENSG00000000460",
      use.weights = TRUE,
      metrics = FALSE
    )
  )

})
