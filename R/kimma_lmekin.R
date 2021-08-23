#' Run kimma linear mixed effects model with kinship
#'
#' @param model.lme Character model created in kmFit
#' @param to.model.gene Data frame formatted in kmFit, subset to gene of interest
#' @param gene Character of gene to model
#' @param kin.subset Pairwise kinship matrix
#'
#' @return Linear model effect with kinship results data frame for 1 gene
#' @keywords internal

kimma_lmekin <- function(model.lme, to.model.gene, gene, kin.subset){
    #Place holder LMEKIN results
    p.kin <- NaN
    sigma.kin <- 0
    results.kin <- NULL

    #Fit LMEKIN model
    fit.kin <- coxme::lmekin(stats::as.formula(model.lme),
                             data=to.model.gene, varlist=as.matrix(kin.subset))
    #Calulate stats
    beta <- fit.kin$coefficients$fixed
    nvar <- length(beta)
    nfrail <- nrow(fit.kin$var) - nvar
    se <- sqrt(diag(fit.kin$var)[nfrail + 1:nvar])
    t <- beta/se
    p.kin <- stats::pchisq((t)^2, 1, lower.tail=FALSE)
    sigma.kin <- fit.kin$sigma

    #Extract results
    results.kin <- data.frame(
      model = rep("lmekin", length(p.kin)),  #Label model as lmekin
      gene = rep(gene, length(p.kin)),       #gene name
      variable = names(p.kin),               #variables in model
      estimate = beta,                       #estimates in model
      pval = p.kin,                          #P-value
      sigma = rep(sigma.kin, length(p.kin))) #sigma

    results.kin.ls <- list()
    results.kin.ls[["fit"]] <- fit.kin
    results.kin.ls[["results"]] <- results.kin
    return(results.kin.ls)
  }
