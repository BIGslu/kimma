test_tomodel_path <- testthat::test_path("test_data", "to_model_ls.Rds")

to.model.ls <- readRDS(test_tomodel_path)
to.model.df <- to.model.ls[["to.model"]]

testthat::test_that("kimma_lmerel produces correct results with gene weights", {

  tst.df <- to.model.df[to.model.df$rowname == "ENSG00000002587", ]

  lmerel.res <- kimma_lmerel(
    model.lme = "expression~virus*asthma+(1|ptID)",
    to.model.gene = tst.df,
    gene = "ENSG00000002587",
    kin.subset = kimma::example.kin,
    use.weights = TRUE,
    patientID = "ptID",
    metrics = TRUE
  )

  res.tab <- lmerel.res$results
  virus.est <- res.tab[res.tab$variable == "virus", "estimate"]
  asthma.est <- res.tab[res.tab$variable == "asthma", "estimate"]
  virus_asthma.est <- res.tab[res.tab$variable == "virus:asthma", "estimate"]
  rand.eff.est <- res.tab[res.tab$variable == "(1 | ptID)", "estimate"]

  testthat::expect_equal(
    object = virus.est,
    expected = -0.8614,
    tolerance = 0.001
  )

  testthat::expect_equal(
    object = asthma.est,
    expected = -0.1164,
    tolerance = 0.001
  )

  testthat::expect_equal(
    object = virus_asthma.est,
    expected = -0.4536,
    tolerance = 0.001
  )

  testthat::expect_equal(
    object = rand.eff.est,
    expected = 3.1392,
    tolerance = 0.001
  )
})
