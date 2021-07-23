#' Run kimma linear mixed effects model
#'
#' @param model.lme Character model created in kmFit
#' @param to.model.gene Data frame formatted in kmFit
#' @param gene Character of gene to model
#'
#' @return Linear model effect results data frame for 1 gene
#' @keywords internal

kimma_lme <- function(model.lme, to.model.gene, gene, contrast.mat, model){
    rowname <- NULL
    #Place holder LME results
    p.lme <- NaN
    sigma.lme <- 0
    results.lme <- NULL

    #Fit LME model
    fit.lme <- lme4::lmer(model.lme, data=to.model.gene)
    #Estimate P-value
    p.lme <- broom::tidy(car::Anova(fit.lme))
    #Get estimates
    est.lme <- as.data.frame(stats::coef(summary(fit.lme))) %>%
      tibble::rownames_to_column() %>%
      dplyr::filter(rowname != "(Intercept)")
    #Calculate sigma
    sigma.lme <- stats::sigma(fit.lme)

    #Extract results
    results.lme <- data.frame(
      model = rep("lme", nrow(p.lme)),    #Label model as lme
      gene = rep(gene, nrow(p.lme)),      #gene name
      variable = p.lme$term,              #variables in model
      estimate = est.lme$Estimate,        #estimate in model
      pval = p.lme$p.value,               #P-value
      sigma = rep(sigma.lme, nrow(p.lme)))#sigma

    results.lme.ls <- list()
    results.lme.ls[["fit"]] <- fit.lme
    results.lme.ls[["results"]] <- results.lme
    return(results.lme.ls)

}