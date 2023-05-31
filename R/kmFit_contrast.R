#' Run pairwise model comparisons with emmeans
#'
#' @param fit model fit from lm( ) or lmer( )
#' @param contrast_var Character vector of variable in model to run contrasts of
#' @param to_model_gene Formatted data from kimma_cleaning( ), subset to gene of interest
#' @param genotype_name Character string. Used internally for kmFit_eQTL
#'
#' @return data frame with contrast model results
#' @keywords internal

kmFit_contrast <- function(fit, contrast_var, to_model_gene, genotype_name){
  contrast.i <- term <- p.value <- contrast_ref <- contrast_lvl <- contrast <- null.value <- estimate <- NULL
  contrast.result <- data.frame()

  #Fit variables in contrast_var
  for(contrast.i in contrast_var){
    contrast.result.temp <- NULL
    #Class
    i.split <- strsplit(contrast.i, split=":")[[1]]
    contrast.is.numeric <- unlist(lapply(to_model_gene[,i.split], is.numeric))

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
        #fix genotype names
        if(!is.null(genotype_name)){
          if(grepl(genotype_name, contrast.i)){
          contrast.result.temp <- contrast.result.temp %>%
            dplyr::mutate(dplyr::across(c(contrast_ref, contrast_lvl),
                                 ~gsub(genotype_name, paste0(genotype_name,"_"), .)))
        }}

        contrast.result <- contrast.result.temp %>%
          dplyr::bind_rows(contrast.result)
      }
    }


  }

  contrast.result.format <- contrast.result %>%
    dplyr::rename(variable=term, pval=p.value) %>%
    dplyr::select(-null.value) %>%
    #Switch estimate sign to match lvl minus ref calculation
    #THE REF AND LVL VALUES ARE INCORRECT UNTIL YOU DO THIS
    dplyr::mutate(estimate = -estimate)

  return(contrast.result.format)
}

