#' Extract lmFit model results
#'
#' Extract model fit and significance for all individual variables and/or contrasts in a limma model
#'
#' @param design model matrix output by model.matrix( )
#' @param fit MArrayLM model fit output by limma::eBayes( )
#' @param contrast.mat contrast matrix output by limma::makeContrasts( ). NOTE: When using constrasts, the result will not exactly match extract_kmFit due to limma's naming of contrast levels as variableLEVEL
#' @param dat.genes data frame with additional gene annotations. Optional. If not provided, the fit object is also checked for gene annotation information.
#' @param name.genes character for variable name in dat.genes that matches gene names in fit
#'
#' @return List with data frames. One for model fit (sigma) and one for significance for all variable and genes. Variables names as in limma::topTable( )
#' @export
#'
#' @examples
#' # Run limma model
#' design <- model.matrix(~ virus, data = example.voom$targets)
#' fit <- limma::eBayes(limma::lmFit(example.voom$E, design))
#'
#' ## Get results
#' fdr <- extract_lmFit(design = design, fit = fit)
#' ## Get results and add gene annotations
#' fdr <- extract_lmFit(design = design, fit = fit,
#'                         dat.genes = example.voom$genes, name.genes = "geneName")
#'
#' # Run limma contrasts model
#' design <- model.matrix(~ 0 + virus, data = example.voom$targets)
#' fit <- limma::lmFit(example.voom$E, design)
#' contrast.mat <- limma::makeContrasts(virusHRV-virusnone, levels = design)
#' fit <- limma::eBayes(limma::contrasts.fit(fit, contrast.mat))
#'
#' ## Get contrast results
#' fdr <- extract_lmFit(design = design, fit = fit, contrast.mat = contrast.mat)

extract_lmFit <- function(design, fit, contrast.mat=NULL,
                          dat.genes=NULL, name.genes="geneName"){
  #Empty df to hold results
  pval.result <- data.frame()
  logFC <- FCgroup <- variable <- geneName <- P.Value <- adj.P.Val <- model <- gene <- contrast_ref <- contrast_lvl <- estimate <- pval <- FDR <- B <- sigma <- NULL

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
    if (!is.numeric(rownames(pval.temp)[1])){
      if(name.genes %in% colnames(pval.temp)){
        pval.temp <- pval.temp %>% dplyr::select(-dplyr::all_of(name.genes))
      }
      pval.temp <- pval.temp %>%
        tibble::rownames_to_column(name.genes) %>%
        dplyr::mutate(variable = vars[var])
    } else{
      pval.temp <- dplyr::mutate(pval.temp, variable = vars[var])
    }

    pval.result <- pval.temp %>%
      dplyr::select(dplyr::all_of(name.genes), logFC:variable) %>%
      dplyr::bind_rows(pval.result)

  }

  #Add hgnc symbol
  if(!is.null(dat.genes)){
    pval.result <- suppressMessages(dplyr::left_join(pval.result, dat.genes))
  } else if(!is.null(fit$genes)){
    pval.result <- suppressMessages(dplyr::left_join(pval.result, fit$genes))
  }

  pval.result.format <- pval.result %>%
    #Convert groups to ordered factors
    dplyr::mutate(variable = factor(variable, levels=vars)) %>%
    #rename to match kmFit
    dplyr::rename(gene=geneName, estimate=logFC, pval=P.Value, FDR=adj.P.Val)

  #Split contrast if present
  if(!is.null(contrast.mat)){
    pval.result.format2 <- pval.result.format %>%
      tidyr::separate(variable, into=c("contrast_lvl","contrast_ref"), sep=" - ") %>%
      dplyr::mutate(variable = paste(Reduce(intersect, strsplit(c(contrast_lvl,
                                                                  contrast_ref),"")),
                                     collapse=""),
                    model="limma.contrast") %>%
      dplyr::select(model, gene, variable, contrast_ref, contrast_lvl,
                    estimate, pval, FDR, t, B, dplyr::everything())

    paste(Reduce(intersect, strsplit(c(pval.result.format2$contrast_lvl,
                                       pval.result.format2$contrast_ref),"")),
          collapse="")
  } else{
    pval.result.format2 <- pval.result.format %>%
      dplyr::mutate(model="limma") %>%
      dplyr::select(model, gene, variable, estimate, pval, FDR, t, B,
                    dplyr::everything())
  }

  #Goodness of fit metric sigma
  lm.fit <- data.frame(
    model = "limma.fit",
    gene = names(fit$Amean),
    sigma = fit$sigma
  )

  ## add gene annotation if exists
  rownames(lm.fit) <- NULL
  if(!is.null(dat.genes)){
    lm.fit <- dplyr::left_join(lm.fit, dat.genes, by=c("gene"=name.genes))
  } else if(!is.null(fit$genes)){
    lm.fit <- dplyr::left_join(lm.fit, fit$genes, by=c("gene"=name.genes))
  }

  result.ls <- list()
  if(is.null(contrast.mat)){
    result.ls[["lm"]] <- pval.result.format2
  } else {
    result.ls[["lm.contrast"]] <- pval.result.format2
  }
  result.ls[["lm.fit"]] <- lm.fit
  return(result.ls)
}
