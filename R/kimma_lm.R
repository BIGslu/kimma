#' Run kimma linear model
#'
#' @param model.lm Character model created in kmFit
#' @param to.model.gene Data frame formatted in kmFit, subset to gene of interest
#' @param gene Character of gene to model
#' @param use.weights Logical if gene specific weights should be used in model. Default is FALSE
#' @param metrics Logical if should calculate model fit metrics such as AIC, BIC, R-squared. Default is FALSE
#'
#' @return Linear model results data frame for 1 gene
#' @keywords internal

kimma_lm <- function(model.lm, to.model.gene, gene, use.weights, metrics){
    #Place holder LM results
    p.lm <- NaN
    sigma.lm <- 0
    results.lm <- NULL
    .GlobalEnv$to.model.gene <- to.model.gene

    #Fit model
    if(use.weights){
      fit.lm <- stats::lm(model.lm, data=to.model.gene, weights=to.model.gene$gene_weight)
    } else{
      fit.lm <- stats::lm(model.lm, data=to.model.gene, weights=NULL)
    }

    p.lm <- broom::tidy(fit.lm)

    #Extract results
    results.lm <- data.frame(
      model = rep("lm", nrow(p.lm)),    #Label model as lm
      gene = rep(gene, nrow(p.lm)),     #gene name
      variable = p.lm$term,             #variables in model
      estimate = p.lm$estimate,         #estimates in model
      pval = p.lm$p.value)              #P-value

    #Model fit metrics
    if(metrics){
      fit.metrics <- data.frame(
        model="lm.fit",
        gene=gene,
        sigma = stats::sigma(fit.lm),
        AIC = stats::AIC(fit.lm),
        BIC = stats::BIC(fit.lm),
        Rsq = summary(fit.lm)$r.squared,
        adj_Rsq = summary(fit.lm)$adj.r.squared)
    } else{
      fit.metrics <- NULL
    }

    results.lm.ls <- list()
    results.lm.ls[["fit"]] <- fit.lm
    results.lm.ls[["metrics"]] <- fit.metrics
    results.lm.ls[["results"]] <- results.lm
    return(results.lm.ls)
}
