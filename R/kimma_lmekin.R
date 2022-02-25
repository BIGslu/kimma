#' Run kimma linear mixed effects model with kinship
#'
#' @param model.lme Character model created in kmFit
#' @param to.model.gene Data frame formatted in kmFit, subset to gene of interest
#' @param gene Character of gene to model
#' @param kin.subset Pairwise kinship matrix
#' @param use.weights Logical if gene specific weights should be used in model. Default is FALSE
#'
#' @return Linear model effect with kinship results data frame for 1 gene
#' @keywords internal

kimma_lmekin <- function(model.lme, to.model.gene, gene, kin.subset, use.weights){
    #Place holder LMEKIN results
    p.kin <- NaN
    sigma.kin <- 0
    results.kin <- NULL
    .GlobalEnv$to.model.gene <- to.model.gene

    #Fit LMEKIN model
    if(use.weights){
      fit.kin <- coxme::lmekin(stats::as.formula(model.lme),
                               data=to.model.gene, varlist=as.matrix(kin.subset),
                               weights=to.model.gene$gene_weight,
                               method="REML")
    } else{
      fit.kin <- coxme::lmekin(stats::as.formula(model.lme),
                               data=to.model.gene, varlist=as.matrix(kin.subset),
                               weights=NULL, method="REML")
    }

    #Calculate stats
    beta <- fit.kin$coefficients$fixed
    nvar <- length(beta)
    nfrail <- nrow(fit.kin$var) - nvar
    se <- sqrt(diag(fit.kin$var)[nfrail + 1:nvar])
    t <- beta/se
    p.kin <- stats::pchisq((t)^2, 1, lower.tail=FALSE)

    #Extract results
    results.kin <- data.frame(
      model = rep("lmekin", length(p.kin)),  #Label model as lmekin
      gene = rep(gene, length(p.kin)),       #gene name
      variable = names(p.kin),               #variables in model
      estimate = beta,                       #estimates in model
      pval = p.kin)                          #P-value

    #Calculate R-squared
    if(use.weights){
      null <- stats::glm(formula = expression ~ 1, data = to.model.gene,
                  weights = to.model.gene$gene_weight)
    } else{
      null <- stats::glm(formula = expression ~ 1, data = to.model.gene)
    }

    L0 <- as.vector(stats::logLik(null))
    L1 <- as.vector(fit.kin$loglik)
    n <- nrow(to.model.gene)
    ret <- 1 - exp(-2 / n * (L1 - L0))
    max.r2 <- 1 - exp(2 / n * L0)

    #Model fit metrics
    fit.metrics <- data.frame(
      model="lmekin.fit",
      gene=gene,
      sigma = fit.kin$sigma,
      AIC = AICcmodavg::AICc(fit.kin, second.ord=FALSE),
      BIC = AICcmodavg::useBIC(fit.kin),
      Rsq = ret,
      adj_Rsq = ret / max.r2
    )

    results.kin.ls <- list()
    results.kin.ls[["fit"]] <- fit.kin
    results.kin.ls[["metrics"]] <- fit.metrics
    results.kin.ls[["results"]] <- results.kin
    return(results.kin.ls)
  }
