#' Summarise lmFit FDR results
#'
#' Summarise number of significant genes at various FDR cutoffs. Can split by up/down fold change as well.
#'
#' @param fdr data.frame output by kimma::extract_lmFit( )
#' @param fdr.cutoff numeric vector of FDR cutoffs to summarise at
#' @param FCgroup logical if should separate summary by up/down fold change groups
#' @param intercept logical if should include intercept variable in summary
#'
#' @return Data frame with total significant genes for each variable at various FDR cutoffs
#' @export
#'
#' @examples
#'# Run limma model
#' design <- model.matrix(~ virus, data = dat.voom.example$targets)
#' fit <- limma::eBayes(limma::lmFit(dat.voom.example$E, design))
#'
#' ## Get results
#' fdr <- extract_lmFit(design = design, fit = fit)
#'
#' # Summarise results
#' fdr.summary <- summarise_lmFit(fdr = fdr, fdr.cutoff = c(0.05, 0.5), FCgroup = TRUE)
#'

summarise_lmFit <- function(fdr, fdr.cutoff = c(0.05,0.1,0.2,0.3,0.4,0.5),
                            FCgroup = FALSE, intercept = FALSE){
  adj.P.Val <- geneName <- group <- n <- variable <- NULL
  if(intercept){
    fdr.filter <- fdr
  } else{
    fdr.filter <- dplyr::filter(fdr, variable != '(Intercept)')
  }

  if(FCgroup){
    #Blank df for results
    result <- data.frame()

    for(FDR.i in fdr.cutoff){
      name.fdr <- paste("fdr",FDR.i, sep="_")
      #Calculate total, nonredundant signif genes at different levels
      total.temp <- fdr.filter %>%
        dplyr::filter(adj.P.Val <= FDR.i) %>%
        dplyr::distinct(geneName, FCgroup) %>%
        dplyr::count(FCgroup, .drop = FALSE) %>%
        dplyr::mutate(variable = "total (nonredundant)")

      #Summarize signif genes per variable at various levels
      group.temp <- fdr.filter %>%
        dplyr::filter(adj.P.Val <= FDR.i) %>%
        dplyr::count(variable, FCgroup, .drop = FALSE)

      result.temp <- dplyr::bind_rows(total.temp, group.temp) %>%
        dplyr::mutate(group = name.fdr)

      result <- dplyr::bind_rows(result, result.temp)
    }
    } else {
      #Blank df for results
      result <- data.frame()

      for(FDR.i in fdr.cutoff){
        name.fdr <- paste("fdr",FDR.i, sep="_")
        #Calculate total, nonredundant signif genes at different levels
        total.temp <- fdr.filter %>%
          dplyr::filter(adj.P.Val <= FDR.i) %>%
          dplyr::distinct(geneName) %>%
          dplyr::mutate(variable = "total (nonredundant)") %>%
          dplyr::count(variable, .drop = FALSE)

        #Summarize signif genes per variable at various levels
        group.temp <- fdr.filter %>%
          dplyr::filter(adj.P.Val <= FDR.i) %>%
          dplyr::count(variable, .drop = FALSE)

        result.temp <- dplyr::bind_rows(total.temp, group.temp) %>%
          dplyr::mutate(group = name.fdr)

        result <- dplyr::bind_rows(result, result.temp)
      }
    }
  #Format to wide output
  result.format <- tidyr::pivot_wider(result, names_from = group, values_from = n) %>%
    dplyr::mutate(variable = forcats::fct_relevel(factor(variable), "total (nonredundant)",
                                                  after=Inf)) %>%
    dplyr::arrange(variable)

  return(result.format)
}
