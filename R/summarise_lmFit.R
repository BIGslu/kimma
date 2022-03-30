#' Summarise lmFit FDR results
#'
#' Summarise number of significant genes at various FDR cutoffs. Can split by up/down fold change as well.
#'
#' @param fdr data.frame output by kimma::extract_lmFit( )
#' @param fdr.cutoff numeric vector of FDR cutoffs to summarise at
#' @param p.cutoff numeric vector of P-value cutoffs to summarise at. No FDR summary given if p.cutoff is provided
#' @param FCgroup logical if should separate summary by up/down fold change groups
#' @param intercept logical if should include intercept variable in summary
#'
#' @return Data frame with total significant genes for each variable at various FDR cutoffs
#' @export
#' @examples
#'# Run limma model
#' design <- model.matrix(~ virus, data = example.voom$targets)
#' fit <- limma::eBayes(limma::lmFit(example.voom$E, design))
#'
#' ## Get results
#' model_results <- extract_lmFit(design = design, fit = fit)
#'
#' # Summarise results
#' summarise_lmFit(fdr = model_results, fdr.cutoff = c(0.05, 0.5), FCgroup = TRUE)
#'
#' # Contrasts model
#' design <- model.matrix(~ 0 + virus, data = example.voom$targets)
#' fit <- limma::lmFit(example.voom$E, design)
#' contrast.mat <- limma::makeContrasts(virusHRV-virusnone, levels = design)
#' fit <- limma::eBayes(limma::contrasts.fit(fit, contrast.mat))
#'
#' ## Get results
#' model_results <- extract_lmFit(design = design, fit = fit,
#'                                contrast.mat = contrast.mat)
#'
#' # Summarise results
#' summarise_lmFit(fdr = model_results, fdr.cutoff = c(0.05, 0.5), FCgroup = FALSE)
#' ## No signif results (not run)
#' # summarise_lmFit(fdr = model_results, fdr.cutoff = 0.0001, FCgroup = FALSE)
#'

summarise_lmFit <- function(fdr, fdr.cutoff = c(0.05,0.1,0.2,0.3,0.4,0.5),
                            p.cutoff = NULL,
                            FCgroup = FALSE, intercept = FALSE){
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
          dplyr::count(contrast_ref, contrast_lvl, FCgroup, .drop = FALSE)

        result.temp <- total.temp %>%
          dplyr::rename(contrast_ref=variable) %>%
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
          dplyr::count(contrast_ref, contrast_lvl, .drop = FALSE)

        result.temp <- total.temp %>%
          dplyr::rename(contrast_ref=variable) %>%
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
  if("contrast_ref" %in% colnames(fdr.filter)){
    result.format <- tidyr::pivot_wider(result, names_from = group,
                                        values_from = n) %>%
      dplyr::mutate(contrast_ref = forcats::fct_relevel(factor(contrast_ref),
                                                        "total (nonredundant)",
                                                    after=Inf)) %>%
      dplyr::arrange(contrast_ref, contrast_lvl)
  } else{
    result.format <- tidyr::pivot_wider(result, names_from = group,
                                        values_from = n) %>%
      dplyr::mutate(variable = forcats::fct_relevel(factor(variable),
                                                    "total (nonredundant)",
                                                    after=Inf)) %>%
      dplyr::arrange(variable)
  }
  #rename for P-value if specified
  if(!is.null(p.cutoff)) {
    result.format <- result.format %>%
      dplyr::rename_all(~gsub("fdr_","p_",.))
  }

  return(result.format)
}

#' @rdname summarise_lmFit
#' @export
summarize_lmFit <- summarise_lmFit
