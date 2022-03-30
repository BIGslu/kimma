#' Run pairwise model comparisons with emmeans
#'
#' @param fit model fit from lm( ) or lmer( )
#' @param contrast.var Character vector of variable in model to run contrasts of
#' @param to.model.gene Formatted data from kimma_cleaning( ), subset to gene of interest
#' @param metrics Logical if should calculate model fit metrics such as AIC, BIC, R-squared. Default is FALSE
#'
#' @return data frame with contrast model results
#' @keywords internal

kmFit_contrast <- function(fit, contrast.var, to.model.gene, metrics){
  contrast.i <- term <- p.value <- contrast_ref <- contrast_lvl <- contrast <- null.value <- NULL
  contrast.result <- data.frame()

 #Fit variables in contrast.var
  for(contrast.i in contrast.var){
    contrast.result.temp <- NULL
    #Class
    i.split <- strsplit(contrast.i, split=":")[[1]]
    contrast.is.numeric <- unlist(lapply(to.model.gene[,i.split], is.numeric))

    #If any numeric
    if(any(contrast.is.numeric)){
      contrast.result.temp <- tryCatch({
        emmeans::emtrends(fit, adjust="none",var=i.split[contrast.is.numeric],
                          stats::as.formula(paste("pairwise~", i.split[!contrast.is.numeric],
                                           sep="")))$contrasts %>%
          broom::tidy() %>%
          dplyr::mutate(term = gsub(":","*", contrast.i)) %>%
          tidyr::separate(contrast, into=c("contrast_ref","contrast_lvl"),
                          sep=" - ")
        }, error=function(e){ return(NULL) })

      #if model ran, add to results
      if(is.data.frame(contrast.result.temp)){
        contrast.result <- contrast.result.temp %>%
          dplyr::bind_rows(contrast.result)
        }
    } else {
      #If not numeric
      contrast.result.temp <- tryCatch({
        emmeans::emmeans(fit, adjust="none",
                         stats::as.formula(paste("pairwise~", contrast.i,sep="")))$contrasts %>%
          broom::tidy() %>%
          dplyr::mutate(term = gsub(":","*", contrast.i)) %>%
          tidyr::separate(contrast, into=c("contrast_ref","contrast_lvl"),
                          sep=" - ")
      }, error=function(e){ return(NULL) })

      #if model ran, add to results
      if(is.data.frame(contrast.result.temp)){
        contrast.result <- contrast.result.temp %>%
          dplyr::bind_rows(contrast.result)
      }
        }}
  contrast.result.format <- contrast.result %>%
    dplyr::rename(variable=term, pval=p.value) %>%
    dplyr::select(-null.value)
  return(contrast.result.format)
}

