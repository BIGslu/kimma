#' Summarise kmFit FDR results
#'
#' Summarise number of significant genes at various FDR cutoffs. Can split by up/down fold change as well.
#'
#' @param fdr data.frame output by kimma::kmFit( ). Main model or contrasts accepted
#' @param fdr.cutoff numeric vector of FDR cutoffs to summarise at
#' @param p.cutoff numeric vector of P-value cutoffs to summarise at. No FDR summary given if p.cutoff is provided
#' @param FCgroup logical if should separate summary by up/down fold change groups
#' @param intercept logical if should include intercept variable in summary
#'
#' @return Data frame with total significant genes for each variable at various FDR cutoffs
#' @export
#'
#' @examples
#' # Run kimma model
#' model_results <- kmFit(dat = example.voom,
#'       kin = example.kin,
#'       run.lme = TRUE, run.lmerel=TRUE, run.contrast=TRUE,
#'       subset.genes = c("ENSG00000250479","ENSG00000250510","ENSG00000255823"),
#'       model = "~ virus + (1|ptID)")
#'
#' # Summarise results
#' summarise_kmFit(fdr = model_results$lmerel, fdr.cutoff = c(0.01, 0.5),
#'                 FCgroup = TRUE)
#' summarise_kmFit(fdr = model_results$lme.contrast, fdr.cutoff = c(0.01, 0.5),
#'                 FCgroup = FALSE)
#'
#' #No significant genes. No run
#' ## summarise_kmFit(fdr = model_results$lmerel, fdr.cutoff = c(0.001))
#'

summarise_kmFit <- function(fdr, fdr.cutoff = c(0.05,0.1,0.2,0.3,0.4,0.5),
                            p.cutoff = NULL, FCgroup = FALSE, intercept = FALSE){
  group <- n <- variable <- fdr.var <- estimate <- gene <- contrast_ref <- contrast_lvl <- NULL
  if(intercept | "contrast_ref" %in% colnames(fdr)){
    fdr.filter <- fdr
  } else{
    fdr.filter <- dplyr::filter(fdr, variable != '(Intercept)') %>%
      droplevels()
  }

  # Use p-values if specified
  if(!is.null(p.cutoff)){
    fdr.cutoff <- p.cutoff
    fdr.var <- "pval"
  } else {
    fdr.var <- "FDR"
  }

  if(FCgroup){
    fdr.filter.FC <- fdr.filter %>%
      dplyr::mutate(FCgroup = ifelse(estimate < 0, "down", "up"))
    #Blank df for results
    result <- data.frame()

    for(FDR.i in fdr.cutoff){
      name.fdr <- paste("fdr",FDR.i, sep="_")
      #Calculate total, nonredundant signif genes at different levels
      total.temp <- fdr.filter.FC %>%
        dplyr::filter(get(fdr.var) <= FDR.i) %>%
        dplyr::distinct(gene, FCgroup) %>%
        dplyr::count(FCgroup, .drop = FALSE) %>%
        dplyr::mutate(variable = "total (nonredundant)")

      #Summarize signif genes per variable at various levels
      if("contrast_ref" %in% colnames(fdr.filter.FC)){
        group.temp <- fdr.filter.FC %>%
          dplyr::filter(get(fdr.var) <= FDR.i) %>%
          dplyr::count(variable, contrast_ref, contrast_lvl, FCgroup, .drop = FALSE)

        result.temp <- total.temp %>%
          dplyr::bind_rows(group.temp) %>%
          dplyr::mutate(group = name.fdr)

        result <- dplyr::bind_rows(result, result.temp)
      } else{
        group.temp <- fdr.filter.FC %>%
          dplyr::filter(get(fdr.var) <= FDR.i) %>%
          dplyr::count(variable, FCgroup, .drop = FALSE)

        result.temp <- dplyr::bind_rows(total.temp, group.temp) %>%
          dplyr::mutate(group = name.fdr)

        result <- dplyr::bind_rows(result, result.temp)
      }
    }
  } else {
    #Blank df for results
    result <- data.frame()

    for(FDR.i in fdr.cutoff){
      name.fdr <- paste("fdr",FDR.i, sep="_")
      #Calculate total, nonredundant signif genes at different levels
      total.temp <- fdr.filter %>%
        dplyr::filter(get(fdr.var) <= FDR.i) %>%
        dplyr::distinct(gene) %>%
        dplyr::mutate(variable = "total (nonredundant)") %>%
        dplyr::count(variable, .drop = FALSE)

      #Summarize signif genes per variable at various levels
      if("contrast_ref" %in% colnames(fdr.filter)){
        group.temp <- fdr.filter %>%
          dplyr::filter(get(fdr.var) <= FDR.i) %>%
          dplyr::count(variable, contrast_ref, contrast_lvl, .drop = FALSE)

        result.temp <- total.temp %>%
          dplyr::bind_rows(group.temp) %>%
          dplyr::mutate(group = name.fdr)

        result <- dplyr::bind_rows(result, result.temp)
      } else{
        group.temp <- fdr.filter %>%
          dplyr::filter(get(fdr.var) <= FDR.i) %>%
          dplyr::count(variable, .drop = FALSE)

        result.temp <- dplyr::bind_rows(total.temp, group.temp) %>%
          dplyr::mutate(group = name.fdr)

        result <- dplyr::bind_rows(result, result.temp)
      }
    }
  }

  #Stop if not signif
  if(nrow(result) == 0){
    stop("No significant genes. Please increase FDR or P-value cutoffs.")
  }

  #Format to wide output
  result.format <- tidyr::pivot_wider(result, names_from = group,
                                        values_from = n) %>%
      dplyr::mutate(variable = forcats::fct_relevel(factor(variable),
                                                    "total (nonredundant)",
                                                    after=Inf)) %>%
      dplyr::arrange(variable)

  #rename for P-value if specified
  if(!is.null(p.cutoff)) {
    result.format <- result.format %>%
      dplyr::rename_all(~gsub("fdr_","p_",.))
  }

  return(result.format)
}

#' @rdname summarise_kmFit
#' @export
summarize_kmFit <- summarise_kmFit
