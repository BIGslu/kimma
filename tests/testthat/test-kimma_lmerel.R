# TODO:
# ?? add tests for kimma_lmerel if:
# if(nrow(p.kin)==nrow(est.kin)){


test_tomodel_path <- testthat::test_path("test_data", "to_model_ls.Rds")

to.model.ls <- readRDS(test_tomodel_path)
to.model.df <- to.model.ls[["to.model"]]
tst.df <- to.model.df[to.model.df$rowname == "ENSG00000002587", ]


testthat::test_that("kimma_lmerel produces correct results with gene weights", {

  lmerel.res <- kimma_lmerel(
    model_lme = "expression~virus*asthma+(1|ptID)",
    to_model_gene = tst.df,
    gene = "ENSG00000002587",
    kin_subset = kimma::example.kin,
    use_weights = TRUE,
    patientID = "ptID",
    metrics = TRUE
  )

  res.tab <- lmerel.res$results
  virus.est <- res.tab[res.tab$variable == "virus", "estimate"]
  asthma.est <- res.tab[res.tab$variable == "asthma", "estimate"]
  virus_asthma.est <- res.tab[res.tab$variable == "virus:asthma", "estimate"]
  rand.eff.est <- res.tab[res.tab$variable == "(1 | ptID)", "estimate"]

  testthat::expect_equal(virus.est, -0.8614, tolerance = 0.001)
  testthat::expect_equal(asthma.est, -0.1164, tolerance = 0.001)
  testthat::expect_equal(virus_asthma.est, -0.4536, tolerance = 0.001)
  testthat::expect_equal(rand.eff.est, 3.1392, tolerance = 0.001)
})


testthat::test_that("kimma_lmerel produces correct results without gene weights", {

  lmerel.res <- kimma_lmerel(
    model_lme = "expression~virus*asthma+(1|ptID)",
    to_model_gene = tst.df,
    gene = "ENSG00000002587",
    kin_subset = kimma::example.kin,
    use_weights = FALSE,
    patientID = "ptID",
    metrics = TRUE
  )

  res.tab <- lmerel.res$results
  virus.est <- res.tab[res.tab$variable == "virus", "estimate"]
  asthma.est <- res.tab[res.tab$variable == "asthma", "estimate"]
  virus_asthma.est <- res.tab[res.tab$variable == "virus:asthma", "estimate"]
  rand.eff.est <- res.tab[res.tab$variable == "(1 | ptID)", "estimate"]

  testthat::expect_equal(virus.est, -0.9075, tolerance = 0.001)
  testthat::expect_equal(asthma.est, -0.1324, tolerance = 0.001)
  testthat::expect_equal(virus_asthma.est, -0.4129, tolerance = 0.001)
  testthat::expect_equal(rand.eff.est, 2.5444, tolerance = 0.001)
})
