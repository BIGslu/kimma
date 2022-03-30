#' Summarise kmFit FDR results
#'
#' Summarise number of significant genes at various FDR cutoffs. Can split by up/down fold change as well.
#'
#' @param fdr data.frame output by kimma::kmFit( )
#' @param fdr.cutoff numeric vector of FDR cutoffs to summarise at
#' @param p.cutoff numeric vector of P-value cutodds to summarise at. No FDR summary given if p.cutoff is provided
#' @param contrast logical if should separate summary by pairwise contrasts within variables
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
#'       run.lme = TRUE, run.lmerel=TRUE,
#'       subset.genes = c("ENSG00000250479","ENSG00000250510","ENSG00000255823"),
#'       model = "~ virus + (1|ptID)")
#'
#' # Summarise results
#' summarise_kmFit(fdr = model_results$lmerel, fdr.cutoff = c(0.01, 0.5),
#'                 FCgroup = TRUE)
#' summarise_kmFit(fdr = model_results$lmerel, fdr.cutoff = c(0.001))
#'

summarise_kmFit <- function(fdr, fdr.cutoff = c(0.05,0.1,0.2,0.3,0.4,0.5),
                            p.cutoff = NULL,
                            contrast = FALSE, FCgroup = FALSE, intercept = FALSE){
  estimate <- variable <- FDR <- gene <- group <- n <- model <- fdr.var <- NULL

  if(intercept){
    fdr.filter <- fdr %>%
      dplyr::mutate(variable=factor(variable)) %>%
      dplyr::mutate(FCgroup = ifelse(estimate<0,"down","up"),
                    FCgroup = factor(FCgroup),
                    model = factor(model))

  } else{
    fdr.filter <- dplyr::filter(fdr, variable != '(Intercept)') %>%
      dplyr::mutate(variable=factor(variable)) %>%
      dplyr::mutate(FCgroup = ifelse(estimate<0,"down","up"),
                    FCgroup = factor(FCgroup),
                    model = factor(model))
  }

  # Use p-values if specified
  if(!is.null(p.cutoff)){
    fdr.cutoff <- p.cutoff
    fdr.var <- "pval"
  } else {
    fdr.var <- "FDR"
  }

  if(FCgroup & contrast){
    #Blank df for results
    result <- data.frame()

    for(FDR.i in fdr.cutoff){
      name.fdr <- paste("fdr",FDR.i, sep="_")
      #Calculate total, nonredundant signif genes at different levels
      total.temp <- fdr.filter %>%
        dplyr::filter(get(fdr.var) <= FDR.i) %>%
        dplyr::distinct(model, gene, contrast, FCgroup) %>%
        dplyr::count(model, FCgroup, .drop = FALSE) %>%
        dplyr::mutate(variable = "total (nonredundant)")

      #Summarize signif genes per variable at various levels
      group.temp <- fdr.filter %>%
        dplyr::filter(get(fdr.var) <= FDR.i) %>%
        dplyr::count(model, variable, contrast, FCgroup, .drop = FALSE)

      result.temp <- dplyr::bind_rows(total.temp, group.temp) %>%
        dplyr::mutate(group = name.fdr,
                      contrast = gsub("contrast","",contrast))

      result <- dplyr::bind_rows(result, result.temp)
    }
  } else if(FCgroup){
    #Blank df for results
    result <- data.frame()

    for(FDR.i in fdr.cutoff){
      name.fdr <- paste("fdr",FDR.i, sep="_")
      #Calculate total, nonredundant signif genes at different levels
      total.temp <- fdr.filter %>%
        dplyr::filter(get(fdr.var) <= FDR.i) %>%
        dplyr::distinct(model, gene, FCgroup) %>%
        dplyr::count(model, FCgroup, .drop = FALSE) %>%
        dplyr::mutate(variable = "total (nonredundant)")

      #Summarize signif genes per variable at various levels
      group.temp <- fdr.filter %>%
        dplyr::filter(get(fdr.var) <= FDR.i) %>%
        dplyr::count(model, variable, FCgroup, .drop = FALSE)

      result.temp <- dplyr::bind_rows(total.temp, group.temp) %>%
        dplyr::mutate(group = name.fdr)

      result <- dplyr::bind_rows(result, result.temp)
    }
  } else if(contrast){
    #Blank df for results
    result <- data.frame()

    for(FDR.i in fdr.cutoff){
      name.fdr <- paste("fdr",FDR.i, sep="_")
      #Calculate total, nonredundant signif genes at different levels
      total.temp <- fdr.filter %>%
        dplyr::filter(get(fdr.var) <= FDR.i) %>%
        dplyr::distinct(model, gene) %>%
        dplyr::mutate(variable = "total (nonredundant)") %>%
        dplyr::count(model, variable, .drop = FALSE)


      #Summarize signif genes per variable at various levels
      group.temp <- fdr.filter %>%
        dplyr::filter(get(fdr.var) <= FDR.i) %>%
        dplyr::count(model, variable, contrast, .drop = FALSE)

      result.temp <- dplyr::bind_rows(total.temp, group.temp) %>%
        dplyr::mutate(group = name.fdr,
                      contrast = gsub("contrast","",contrast))

      result <- dplyr::bind_rows(result, result.temp)
    }} else {
    #Blank df for results
    result <- data.frame()

    for(FDR.i in fdr.cutoff){
      name.fdr <- paste("fdr",FDR.i, sep="_")
      #Calculate total, nonredundant signif genes at different levels
      total.temp <- fdr.filter %>%
        dplyr::filter(get(fdr.var) <= FDR.i) %>%
        dplyr::distinct(model, gene) %>%
        dplyr::mutate(variable = "total (nonredundant)") %>%
        dplyr::count(model, variable, .drop = FALSE)%>%
        dplyr::mutate(variable = "total (nonredundant)")

      #Summarize signif genes per variable at various levels
      group.temp <- fdr.filter %>%
        dplyr::filter(get(fdr.var) <= FDR.i) %>%
        dplyr::count(model, variable, .drop = FALSE)

      #If 0 genes, make zero df
      if(nrow(total.temp)==0){
        result.temp <- data.frame(variable=c(unique(fdr.filter$variable),
                                             "total (nonredundant)"),
                                  n=0,
                                  model=unique(fdr.filter$model),
                                  group=name.fdr)
      } else {
        result.temp <- dplyr::bind_rows(total.temp, group.temp) %>%
          dplyr::mutate(group = name.fdr)
      }
      result <- dplyr::bind_rows(result, result.temp)
    }
  }

  #Format to wide output
  if("total (nonredundant)" %in% result$variable){
  result.format <- tidyr::pivot_wider(result, names_from = group, values_from = n) %>%
    dplyr::mutate(variable = forcats::fct_relevel(variable,
                                                  "total (nonredundant)",
                                                  after=Inf)) %>%
    dplyr::arrange(model, variable)
  } else{
    result.format <- tidyr::pivot_wider(result, names_from = group, values_from = n) %>%
      dplyr::arrange(model, variable)
  }

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
