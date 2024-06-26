#' Summarise kmFit FDR results
#'
#' Summarise number of significant genes at various FDR cutoffs. Can split by up/down fold change as well.
#'
#' @param fdr data.frame output by kimma::kmFit( ). Main model or contrasts accepted
#' @param fdr_cutoff numeric vector of FDR cutoffs to summarise at
#' @param p_cutoff numeric vector of P-value cutoffs to summarise at. No FDR summary given if p.cutoff is provided
#' @param FCgroup logical if should separate summary by up/down fold change groups
#' @param intercept logical if should include intercept variable in summary

#' @param fdr.cutoff Deprecated form of fdr_cutoff
#' @param p.cutoff Deprecated form of p_cutoff
#'
#' @return Data frame with total significant genes for each variable at various FDR cutoffs
#' @export
#'
#' @examples
#' # Run kimma model
#' model_results <- kmFit(dat = example.voom,
#'       kin = example.kin,
#'       run_lme = TRUE, run_lmerel=TRUE, run_contrast=TRUE,
#'       subset_genes = c("ENSG00000250479","ENSG00000250510","ENSG00000255823"),
#'       model = "~ virus + asthma + (1|ptID)")
#'
#' # Or extract limma results
#' # design <- model.matrix(~ virus + asthma, data = example.voom$targets)
#' # fit <- limma::eBayes(limma::lmFit(example.voom$E, design))
#' # model_results <- extract_lmFit(design = design, fit = fit)
#'
#' # Summarise results
#' summarise_kmFit(fdr = model_results$lmerel, fdr_cutoff = c(0.01, 0.5),
#'                 FCgroup = TRUE)
#' summarise_kmFit(fdr = model_results$lme.contrast, fdr_cutoff = c(0.01, 0.5),
#'                 FCgroup = FALSE)
#'
#' #No significant genes. No run
#' ## summarise_kmFit(fdr = model_results$lmerel, fdr_cutoff = c(0.001))
#'

summarise_kmFit <- function(fdr, fdr_cutoff = c(0.05,0.1,0.2,0.3,0.4,0.5),
                            p_cutoff = NULL, FCgroup = FALSE, intercept = FALSE,
                            fdr.cutoff = NULL, p.cutoff = NULL){
  group <- n <- variable <- fdr.var <- estimate <- gene <- contrast_ref <- contrast_lvl <- NULL

  #Back compatibility
  if(!is.null(fdr.cutoff)){fdr_cutoff <- fdr.cutoff}
  if(!is.null(p.cutoff)){p_cutoff <- p.cutoff}

  if(intercept | "contrast_ref" %in% colnames(fdr)){
    fdr.filter <- fdr %>%
      dplyr::mutate(variable = as.factor(variable))
  } else{
    fdr.filter <- dplyr::filter(fdr, variable != '(Intercept)') %>%
      dplyr::mutate(variable = as.factor(variable)) %>%
      droplevels()
  }

  # Use p-values if specified
  if(!is.null(p_cutoff)){
    fdr_cutoff <- p_cutoff
    fdr.var <- "pval"
  } else {
    fdr.var <- "FDR"
  }

  if(FCgroup){
    fdr.filter.FC <- fdr.filter %>%
      dplyr::mutate(FCgroup = dplyr::case_when(estimate < 0 ~ "down",
                                               estimate > 0 ~ "up",
                                               TRUE ~ "0"),
                    FCgroup = as.factor(FCgroup))
    #Blank df for results
    result <- data.frame()

    for(FDR.i in fdr_cutoff){
      name.fdr <- paste("fdr",FDR.i, sep="_")
      #Calculate total, nonredundant signif genes at different levels
      total.temp <- fdr.filter.FC %>%
        dplyr::filter(get(fdr.var) <= FDR.i) %>%
        droplevels() %>%
        dplyr::distinct(gene, FCgroup) %>%
        dplyr::count(FCgroup, .drop = FALSE) %>%
        dplyr::mutate(variable = "total (nonredundant)")

      #Summarize signif genes per variable at various levels
      if("contrast_ref" %in% colnames(fdr.filter.FC)){
        group.temp <- fdr.filter.FC %>%
          dplyr::mutate(dplyr::across(c(contrast_ref, contrast_lvl),
                                      ~as.factor(.))) %>%
          dplyr::filter(get(fdr.var) <= FDR.i) %>%
          dplyr::count(variable, contrast_ref, contrast_lvl, FCgroup, .drop = FALSE) %>%
          dplyr::inner_join(dplyr::distinct(fdr.filter.FC,
                                            variable, contrast_ref, contrast_lvl),
                            by = c("variable", "contrast_ref", "contrast_lvl"))

        result.temp <- total.temp %>%
          dplyr::bind_rows(group.temp) %>%
          dplyr::mutate(group = name.fdr)

        result <- dplyr::bind_rows(result, result.temp)
      } else{
        group.temp <- fdr.filter.FC %>%
          dplyr::filter(get(fdr.var) <= FDR.i) %>%
          dplyr::count(variable, FCgroup, .drop = FALSE) %>%
          dplyr::filter(FCgroup != "0" | n > 0)

        result.temp <- dplyr::bind_rows(total.temp, group.temp) %>%
          dplyr::mutate(group = name.fdr)

        result <- dplyr::bind_rows(result, result.temp)
      }
    }
  } else {
    #Blank df for results
    result <- data.frame()

    for(FDR.i in fdr_cutoff){
      name.fdr <- paste("fdr",FDR.i, sep="_")
      #Calculate total, nonredundant signif genes at different levels
      total.temp <- fdr.filter %>%
        dplyr::filter(get(fdr.var) <= FDR.i) %>%
        droplevels() %>%
        dplyr::distinct(gene) %>%
        dplyr::mutate(variable = "total (nonredundant)") %>%
        dplyr::count(variable, .drop = FALSE)

      #Summarize signif genes per variable at various levels
      if("contrast_ref" %in% colnames(fdr.filter)){
        group.temp <- fdr.filter %>%
          dplyr::mutate(dplyr::across(c(contrast_ref, contrast_lvl),
                                      ~as.factor(.))) %>%
          dplyr::filter(get(fdr.var) <= FDR.i) %>%
          dplyr::count(variable, contrast_ref, contrast_lvl, .drop = FALSE) %>%
          dplyr::inner_join(dplyr::distinct(fdr.filter,
                                            variable, contrast_ref, contrast_lvl),
                            by = c("variable", "contrast_ref", "contrast_lvl"))

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
  if(sum(result$n) == 0){
    stop("No significant genes. Please increase FDR or P-value cutoffs.")
  }

  #Format to wide output
  result.format <- tidyr::pivot_wider(result, names_from = group,
                                      values_from = n) %>%
    dplyr::mutate(variable = forcats::fct_relevel(factor(variable),
                                                  "total (nonredundant)",
                                                  after=Inf)) %>%
    dplyr::arrange(variable) %>%
    dplyr::mutate(dplyr::across(dplyr::contains("fdr"),
                                ~tidyr::replace_na(., 0))) %>%
    dplyr::rename_all(~gsub("fdr_","fdr<",.))

    #rename for P-value if specified
    if(!is.null(p_cutoff)) {
      result.format <- result.format %>%
        dplyr::rename_all(~gsub("fdr_","p_",.))
    }

  return(result.format)
}

#' @rdname summarise_kmFit
#' @export
summarize_kmFit <- summarise_kmFit
#' @rdname summarise_kmFit
#' @export
summarise_lmFit <- summarise_kmFit
#' @rdname summarise_kmFit
#' @export
summarize_lmFit <- summarise_kmFit
