#' Run pairwise model comparisons with emmeans
#'
#' @param fit model fit from lm( ) or lmer( )
#' @param contrast.vars character vector of variables to test
#' @param contrast.list contrast matrix formatted as list
#'
#' @return data frame with contrast model results
#' @keywords internal

kmFit_contrast <- function(fit, contrast.vars, contrast.list){
  contrast.lme <- data.frame()

  if(is.null(contrast.list)){
    ##All contrasts in variables
    for(var.i in contrast.vars){
      contrast.lme <- emmeans::emmeans(fit, var.i) %>%
        emmeans::contrast(method="pairwise") %>%
        broom::tidy() %>%
        dplyr::bind_rows(contrast.lme)
      }
    } else if(!is.null(contrast.list)){
      ##From model matrix
      for(var.i in contrast.vars){
        contrast.lme <- emmeans::lsmeans(fit, var.i) %>%
          emmeans::contrast(contrast.list) %>%
          broom::tidy() %>%
          dplyr::bind_rows(contrast.lme)
      }}

  return(contrast.lme)
}

