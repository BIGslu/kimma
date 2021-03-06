#' Extract lmFit model results
#'
#' Extract model fit and significance for all individual variables and/or contrasts in a limma model
#'
#' @param design model matrix output by model.matrix( )
#' @param fit MArrayLM model fit output by limma::eBayes( )
#' @param contrast.mat contrast matrix output by limma::makeContrasts( ). NOTE: When using constrasts, the result will not exactly match extract_kmFit due to limma's naming of contrast levels as variableLEVEL
#' @param dat.genes data frame with additional gene annotations. Optional.
#' @param name.genes character for variable name in dat.genes that matches gene names in fit
#'
#' @return Data frame with model fit and significance for all variable and genes. Format as in limma::topTable( )
#' @export
#'
#' @examples
#' # Run limma model
#' design <- model.matrix(~ virus, data = example.voom$targets)
#' fit <- limma::eBayes(limma::lmFit(example.voom$E, design))
#'
#' ## Get results
#' result <- extract_lmFit(design = design, fit = fit)
#' ## Get results and add gene annotations
#' fdr <- extract_lmFit(design = design, fit = fit,
#'                         dat.genes = example.voom$genes, name.genes = "geneName")
#'
#' # Run limma contrasts model
#' design <- model.matrix(~ 0 + virus, data = example.voom$targets)
#' fit <- limma::lmFit(example.voom$E, design)
#' contrast.mat <- limma::makeContrasts(virusHRV-virusnone, levels = design)
#' fit <- eBayes(contrasts.fit(fit, contrast.mat))
#'
#' ## Get contrast results
#' fdr <- extract_lmFit(design = design, fit = fit, contrast.mat = contrast.mat)

extract_lmFit <- function(design, fit, contrast.mat=NULL,
                          dat.genes=NULL, name.genes="geneName"){
  #Empty df to hold results
  pval.result <- data.frame()
  logFC <- FCgroup <- variable <- geneName <- P.Value <- adj.P.Val <- model <- gene <- contrast_ref <- contrast_lvl <- estimate <- pval <- FDR <- B <- NULL

  #List variables of interest
  if(!is.null(contrast.mat)){
    vars <- colnames(contrast.mat)
  } else{
    vars <- colnames(design)
  }

  for(var in 1:length(vars)){

    pval.temp <- limma::topTable(fit, coef=var,
                                 number=nrow(fit),
                                 adjust.method = "BH")

    #Move gene names from rownames if exist
    if (is.numeric(pval.temp[,1])){
      pval.temp <- pval.temp %>%
        tibble::rownames_to_column("geneName") %>%
        dplyr::mutate(variable = vars[var])
    } else{
      pval.temp <- dplyr::mutate(pval.temp, variable = vars[var])
    }

    pval.result <- dplyr::bind_rows(pval.result, pval.temp)

  }

  #Add hgnc symbol
  if(!is.null(dat.genes)){
    pval.result <- dplyr::left_join(pval.result, dat.genes, by=name.genes)
  }

  pval.result.format <- pval.result %>%
    #Convert groups to ordered factors
    dplyr::mutate(variable = factor(variable, levels=vars)) %>%
    #rename to match kmFit
    dplyr::rename(gene=geneName, estimate=logFC, pval=P.Value, FDR=adj.P.Val) %>%
    dplyr::mutate(model="limma")

  #Split contrast if present
  if(!is.null(contrast.mat)){
    pval.result.format2 <- pval.result.format %>%
      tidyr::separate(variable, into=c("contrast_lvl","contrast_ref"), sep=" - ") %>%
      dplyr::select(model, gene, contrast_ref, contrast_lvl,
                    estimate, pval, FDR, t, B, dplyr::everything())

      paste(Reduce(intersect, strsplit(c(pval.result.format2$contrast_lvl,
                                   pval.result.format2$contrast_ref),"")),
            collapse="")
  } else{
    pval.result.format2 <- pval.result.format %>%
      dplyr::select(model, gene, variable, estimate, pval, FDR, t, B,
                    dplyr::everything())
  }

  return(pval.result.format2)
}
