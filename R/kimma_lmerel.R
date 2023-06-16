#' Run kimma linear mixed effects model with kinship
#'
#' @param model_lme Character model created in kmFit
#' @param to_model_gene Data frame formatted in kmFit, subset to gene of interest
#' @param gene Character of gene to model
#' @param kin_subset Pairwise kinship matrix
#' @param use_weights Logical if gene specific weights should be used in model. Default is FALSE
#' @param patientID Character of variable name to match dat$targets to kinship row and column names.
#' @param metrics Logical if should calculate model fit metrics such as AIC, BIC, R-squared. Default is FALSE
#'
#' @return Linear mixed effect with kinship results data frame for 1 gene
#' @keywords internal

kimma_lmerel <- function(model_lme, to_model_gene, gene, kin_subset, use_weights,
                         patientID, metrics){
  rowname <- NULL
  #Place holder LME results
  p.kin <- NaN
  sigma.kin <- 0
  results.kin <- NULL
  .GlobalEnv$to_model_gene <- to_model_gene

  #Format kinship matrix to list
  kin.ls <- list()
  kin.ls[[patientID]] <- as.matrix(kin_subset)
  #Fit kin model
    if(use.weights){
    fit.kin <- lme4qtl::relmatLmer(model.lme, data=to.model.gene,
                                   relmat = kin.ls,
                                   weights=to_model_gene$gene_weight)
  } else{
    fit.kin <- lme4qtl::relmatLmer(model.lme, data=to.model.gene,
                                   relmat = kin.ls)
  }

  #Estimate P-value
  p.kin <- broom::tidy(car::Anova(fit.kin))
  p.rand <- as.data.frame(lmerTest::rand(fit.kin)) %>%
    tibble::rownames_to_column() %>%
    dplyr::filter(rowname != "<none>")
  #Get estimates
  est.kin <- as.data.frame(stats::coef(summary(fit.kin))) %>%
    tibble::rownames_to_column() %>%
    dplyr::filter(rowname != "(Intercept)")

  #Extract results
  #If no 3+ level variables
  if(nrow(p.kin)==nrow(est.kin)){
    results.kin <- data.frame(
      model = "lmerel",                               #Label model as lme
      gene = gene,                                    #gene name
      variable = c(p.kin$term, p.rand$rowname),       #variables in model
      estimate = c(est.kin$Estimate, p.rand$LRT),     #estimate in model
      statistic = c(p.kin$statistic, rep(NA, nrow(p.rand))),  #test statistic
      #df = NA,                                        #degrees of freedom
      pval = c(p.kin$p.value, p.rand$`Pr(>Chisq)`))   #P-value
  } else{
    #If 3+ variable
    results.kin <- data.frame(
      model = "lmerel",                             #Label model as lme
      gene = gene,                                  #gene name
      variable = c(p.kin$term, p.rand$rowname),     #variables in model
      statistic = "seeContrasts",                   #test statistic
      #df = "seeContrasts",                          #degrees of freedom
      estimate = "seeContrasts",                    #estimate in model
      pval = c(p.kin$p.value, p.rand$`Pr(>Chisq)`)) #P-value
  }

  if(metrics){
    #Calculate R-squared
    if(use_weights){
      null <- stats::glm(formula = expression ~ 1, data = to_model_gene,
                         weights = to_model_gene$gene_weight)
    } else{
      null <- stats::glm(formula = expression ~ 1, data = to_model_gene)
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
      adj_Rsq = ret / max.r2)
  } else{
    fit.metrics <- NULL
  }

  results.kin.ls <- list()
  results.kin.ls[["fit"]] <- fit.kin
  results.kin.ls[["metrics"]] <- fit.metrics
  results.kin.ls[["results"]] <- results.kin
  return(results.kin.ls)

}
