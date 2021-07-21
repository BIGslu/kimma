#' Run kimma linear model
#'
#' @param model.lm Character model created in kmFit
#' @param to.model.gene Data frame formatted in kmFit
#' @param gene Character of gene to model
#'
#' @return Linear model results data frame for 1 gene
#' @keywords internal

kimma_lm <- function(model.lm, to.model.gene, gene, contrast.mat, model){
    #Place holder LM results
    p.lm <- NaN
    sigma.lm <- 0
    results.lm <- NULL

    #Fit model
    fit.lm <- stats::lm(model.lm, data=to.model.gene)
    p.lm <- broom::tidy(fit.lm)
    sigma.lm <- stats::sigma(fit.lm)

    #Extract results
    results.lm <- data.frame(
      model = rep("lm", nrow(p.lm)),    #Label model as lm
      gene = rep(gene, nrow(p.lm)),     #gene name
      variable = p.lm$term,             #variables in model
      estimate = p.lm$estimate,         #estimates in model
      pval = p.lm$p.value,              #P-value
      sigma = rep(sigma.lm, nrow(p.lm)))#sigma

    results.lm.ls <- list()
    results.lm.ls[["fit"]] <- fit.lm
    results.lm.ls[["results"]] <- results.lm
    return(results.lm.ls)
}
