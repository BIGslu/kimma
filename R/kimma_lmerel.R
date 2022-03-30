#' Run kimma linear mixed effects model with kinship
#'
#' @param model.lme Character model created in kmFit
#' @param to.model.gene Data frame formatted in kmFit, subset to gene of interest
#' @param gene Character of gene to model
#' @param kin.subset Pairwise kinship matrix
#' @param use.weights Logical if gene specific weights should be used in model. Default is FALSE
#' @param patientID Character of variable name to match dat$targets to kinship row and column names.

kimma_lmerel <- function(model.lme, to.model.gene, gene, kin.subset, use.weights,
                         patientID){
  rowname <- NULL
  #Place holder LME results
  p.kin <- NaN
  sigma.kin <- 0
  results.kin <- NULL
  .GlobalEnv$to.model.gene <- to.model.gene

  #Format kinship matrix to list
  kin.ls <- list()
  kin.ls[[patientID]] <- as.matrix(kin.subset)
  #Fit kin model
    if(use.weights){
    fit.kin <- lme4qtl::relmatLmer(model.lme, data=to.model.gene,
                                   relmat = kin.ls,
                                   weights=to.model.gene$gene_weight)
  } else{
    fit.kin <- lme4qtl::relmatLmer(model.lme, data=to.model.gene,
                                   relmat = kin.ls)
  }

  #Estimate P-value
  p.kin <- broom::tidy(car::Anova(fit.kin))
  p.rand <- as.data.frame(lmerTest::rand(fit.kin))
  #Get estimates
  est.kin <- as.data.frame(stats::coef(summary(fit.kin))) %>%
    tibble::rownames_to_column() %>%
    dplyr::filter(rowname != "(Intercept)")
  #Calculate sigma
  sigma.kin <- stats::sigma(fit.kin)

  #Extract results
  #If no 3+ level variables
  if(nrow(p.kin)==nrow(est.kin)){
    results.kin <- data.frame(
      model = "lmerel",                      #Label model as lmerel
      gene = gene,                        #gene name
      variable = c(p.kin$term, rownames(p.rand)[-1]),      #variables in model
      estimate = c(est.kin$Estimate, p.rand$LRT[-1]),      #estimate in model
      pval = c(p.kin$p.value, p.rand$`Pr(>Chisq)`[-1]))    #P-value
  } else {
    #If 3+ variable
    results.kin <- data.frame(
      model = "lmerel",                      #Label model as lmerel
      gene = gene,                        #gene name
      variable = c(p.kin$term, p.rand$t[-1]),      #variables in model
      estimate = "seeContrasts",                  #estimate in model
      pval = c(p.kin$p.value, p.rand$p.value[-1])) #P-value
  }

  #Calculate R-squared
  if(use.weights){
    null <- stats::glm(formula = expression ~ 1, data = to.model.gene,
                       weights = to.model.gene$gene_weight)
  } else{
    null <- stats::glm(formula = expression ~ 1, data = to.model.gene)
  }

  L0 <- as.vector(stats::logLik(null))
  L1 <- as.vector(stats::logLik(fit.kin))
  n <- stats::nobs(fit.kin)
  ret <- 1 - exp(-2 / n * (L1 - L0))
  max.r2 <- 1 - exp(2 / n * L0)

  #Model fit metrics
  fit.metrics <- data.frame(
    model="lmerel.fit",
    gene=gene,
    sigma = stats::sigma(fit.kin),
    AIC = stats::AIC(fit.kin),
    BIC = stats::BIC(fit.kin),
    Rsq = ret,
    adj_Rsq = ret / max.r2
  )

  results.kin.ls <- list()
  results.kin.ls[["fit"]] <- fit.kin
  results.kin.ls[["metrics"]] <- fit.metrics
  results.kin.ls[["results"]] <- results.kin
  return(results.kin.ls)

}
