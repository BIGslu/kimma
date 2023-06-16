#' Run kimma linear mixed effects model
#'
#' @param model_lme Character model created in kmFit
#' @param to_model_gene Data frame formatted in kmFit, subset to gene of interest
#' @param gene Character of gene to model
#' @param use_weights Logical if gene specific weights should be used in model. Default is FALSE
#' @param metrics Logical if should calculate model fit metrics such as AIC, BIC, R-squared. Default is FALSE
#'
#' @return Linear mixed effect results data frame for 1 gene
#' @keywords internal

kimma_lme <- function(model_lme, to_model_gene, gene, use_weights, metrics){
    rowname <- NULL
    #Place holder LME results
    p.lme <- NaN
    sigma.lme <- 0
    results.lme <- NULL
    .GlobalEnv$to_model_gene <- to_model_gene

    #Fit LME model
    if(use.weights){
      fit.lme <- lme4::lmer(model.lme, data=to.model.gene,
                            weights=to.model.gene$gene_weight)
    } else{
      fit.lme <- lme4::lmer(model.lme, data=to.model.gene, weights=NULL)
    }

    #Estimate P-value
    p.lme <- broom::tidy(car::Anova(fit.lme))
    p.rand <- as.data.frame(lmerTest::rand(fit.lme)) %>%
      tibble::rownames_to_column() %>%
      dplyr::filter(rowname != "<none>")
    #Get estimates
    est.lme <- as.data.frame(stats::coef(summary(fit.lme))) %>%
      tibble::rownames_to_column() %>%
      dplyr::filter(rowname != "(Intercept)")

    #Extract results
    #If no 3+ level variables
    if(nrow(p.lme)==nrow(est.lme)){
      results.lme <- data.frame(
        model = "lme",                                  #Label model as lme
        gene = gene,                                    #gene name
        variable = c(p.lme$term, p.rand$rowname),       #variables in model
        statistic = c(p.lme$statistic, rep(NA, nrow(p.rand))),  #test statistic
        #df = NA,                                        #degrees of freedom
        estimate = c(est.lme$Estimate, p.rand$LRT),     #estimate in model
        pval = c(p.lme$p.value, p.rand$`Pr(>Chisq)`))   #P-value
    } else{
      #If 3+ variable
      results.lme <- data.frame(
        model = "lme",                                #Label model as lme
        gene = gene,                                  #gene name
        variable = c(p.lme$term, p.rand$rowname),     #variables in model
        statistic = "seeContrasts",                   #test statistic
        #df = "seeContrasts",                          #degrees of freedom
        estimate = "seeContrasts",                    #estimate in model
        pval = c(p.lme$p.value, p.rand$`Pr(>Chisq)`)) #P-value
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
      L1 <- as.vector(stats::logLik(fit.lme))
      n <- stats::nobs(fit.lme)
      ret <- 1 - exp(-2 / n * (L1 - L0))
      max.r2 <- 1 - exp(2 / n * L0)

      #Model fit metrics
      fit.metrics <- data.frame(
        model="lme.fit",
        gene=gene,
        sigma = stats::sigma(fit.lme),
        AIC = stats::AIC(fit.lme),
        BIC = stats::BIC(fit.lme),
        Rsq = ret,
        adj_Rsq = ret / max.r2
      )
    } else{
      fit.metrics <- NULL
    }

    results.lme.ls <- list()
    results.lme.ls[["fit"]] <- fit.lme
    results.lme.ls[["metrics"]] <- fit.metrics
    results.lme.ls[["results"]] <- results.lme
    return(results.lme.ls)

}
