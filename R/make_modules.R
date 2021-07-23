#' Construct WGCNA modules and associated data
#'
#' Make WGCNA modules from gene expression data with dynamic soft threshold selection. Also outputs mean module expression and DAVID formatted gene lists
#'
#' @param dat limma EList output by voom( )
#' @param genes Character vector of genes to used in module building. Must match rownames in dat. If not set, all genes in dat are used
#' @param Rsq.min Numeric minimum R-squared for soft threshold selection. If set, sft.value is not used
#' @param sft.value Numeric soft threshold. Set when minimum R-squared is no used
#' @param minModuleSize Numeric minimum module size
#' @param deepSplit Integer value between 0 and 4. Provides a simplified control over how sensitive module detection should be to module splitting, with 0 least and 4 most sensitive
#' @param nThread Integer for number of threads to use
#'
#' @return List including:
#' \itemize{
#'   \item{genes} Character vector of genes used in module building
#'   \item{sft} Data frame with soft thresholding selected for module building. Includes power, minimum R-squared, and connectivity
#'   \item{sft.plot} ggplot object of soft thresholding topology and connectivity
#'   \item{mods} Data frame of genes in modules
#'   \item{mods.voom} Data frame of mean module expression in each library
#'   \item{david} DAVID formatted data frame of genes in modules
#' }
#' @export
#'
#' @examples
#' dat.mods <- make_modules(dat = dat.voom.example, sft.value = 1)

make_modules <- function(dat, genes=NULL,
                         Rsq.min = NULL, sft.value = NULL, minModuleSize = 20,
                         deepSplit = 3, nThread=2){

  ##### Check data format #####
  if(!is.numeric(dat$E) & is.numeric(dat$E[,1])){
    stop("Gene names must be rownames or the first column of dat.")
  }
  rowname <- SFT.R.sq <- Power <- slope <- fit <- mean.k. <- value <- cutoff <- module <- module.char <- geneName <- NULL

  ##### Subset significant genes #####
  dat.signif <- dat
  #If gene names as rownames
  if(!is.null(genes) & is.numeric(dat$E)){
    dat.signif$E <- as.data.frame(dat.signif$E) %>%
      tibble::rownames_to_column() %>%
      dplyr::filter(rowname %in% genes) %>%
      tibble::column_to_rownames()
    #If gene names are within df
  } else if(!is.null(genes) & !is.numeric(dat$E[,1])){
    gene.column <- colnames(dat$E)[1]

    dat.signif$E <- as.data.frame(dat.signif$E) %>%
      dplyr::filter(get(gene.column) %in% genes) %>%
      tibble::column_to_rownames(gene.column)
  }

  ##### Soft-thresholding #####
  WGCNA::allowWGCNAThreads(nThreads=nThread)
  #Calculate
  sft <- WGCNA::pickSoftThreshold(t(dat.signif$E),
                           powerVector=c(1:30), verbose=0,
                           networkType = "signed")
  sft.format <- as.data.frame(sft$fitIndices)

  #Select threshold
  if(!is.null(Rsq.min)){
    sft.select <- dplyr::filter(sft.format, SFT.R.sq  >= Rsq.min)
    if(nrow(sft.select)>0){
      power.t <- min(sft.select$Power)
    } else {
      stop("R-squared minimum not reached. Please input lower Rsq.min or set sft.value instead.")
    }
  } else if(!is.null(sft.value)){
    sft.select <- dplyr::filter(sft.format, Power == sft.value)
    power.t <- unique(sft.select$Power)
  } else{
    stop("Please Rsq.min or sft.value.")
  }

  #Plot
  sft.plot <- sft.format %>%
    dplyr::mutate(fit = -sign(slope)*SFT.R.sq) %>%
    tidyr::pivot_longer(c(fit,mean.k.), names_to = "name", values_to = "value") %>%

    ggplot2::ggplot(ggplot2::aes(x=Power, y=value, label=Power)) +
    ggplot2::geom_text(size=3) +
    ggplot2::facet_wrap(~name, scales = "free_y",
               labeller=ggplot2::labeller(name=
                                   c(fit="Scale free topology model fit,\nsigned R^2",
                                     mean.k.="Mean connectivity"))) +
    ggplot2::theme_classic() +
    ggplot2::labs(y="", x="Soft threshold (power)") +
    ggplot2::geom_hline(data=data.frame(name=c("fit","mean.k."),
                               cutoff=c(sft.select$SFT.R.sq[1],
                                        sft.select$mean.k.[1])),
               ggplot2::aes(yintercept = cutoff), color="red") +
    ggplot2::geom_vline(xintercept = power.t, color="red")

  ##### Build modules #####
  cor <- WGCNA::cor #Reassign cor fxn to prevent namespace error with stats
  mod.net <- WGCNA::blockwiseModules(t(dat.signif$E),
                              power=power.t,
                              networkType="signed",
                              TOMType="signed",
                              maxBlockSize=500,
                              minModuleSize=minModuleSize,
                              deepSplit=deepSplit,
                              numericLabels=TRUE,
                              saveTOMFileBase="TOM-blockwise",
                              nthreads=nThread)
  cor <- stats::cor #Revert cor fxn assignment

  mods <- as.data.frame(mod.net$colors) %>%
    tibble::rownames_to_column("geneName") %>%
    dplyr::rename(module = "mod.net$colors") %>%
    dplyr::left_join(dat.signif$genes, by="geneName") %>%
    #add leading 0 to module names for correct sorting of factor
    dplyr::mutate(module.char = ifelse(module <= 9,
                                paste("0", module, sep=""),
                                module)) %>%
    #Add color var
    dplyr::mutate(mod.color = WGCNA::labels2colors(mod.net$colors))

  ##### Mean module expression #####
  mods.voom <- mods %>%
    #Combine count and module data
    dplyr::select(geneName, module.char) %>%
    dplyr::left_join(tibble::rownames_to_column(as.data.frame(dat.signif$E), "geneName"),
              by="geneName") %>%
    #Calculate mean by module
    dplyr::group_by(module.char) %>%
    dplyr::summarise_if(is.numeric, mean, na.rm = TRUE) %>%
    tibble::rownames_to_column() %>%
    #Make module names
    dplyr::mutate(rowname=paste("module", module.char, sep="_")) %>%
    tibble::column_to_rownames()

  ##### DAVID format #####
  #list modules
  mod.names <- sort(rownames(mods.voom))
  #calculate maximum module gene list length
  max.mod <- max(table(mod.net$colors))
  #Create data frame with that number of rows
  david <- data.frame(rowname=1:max.mod)

  for (i in 1:length(mod.names)){
    #Create module name as "module#"
    mod.name <- mod.names[i]
    #Filter gene list to module of interest
    gene.list <- mods %>%
      dplyr::filter(module == i-1) %>%
      dplyr::select(geneName)
    #Calculate the total number of rows that need to be added for nrows to match
    add.genes <- max.mod - nrow(gene.list)
    #Add NAs to fill out gene.list and convert to data frame for merging
    gene.list <- c(gene.list$geneName, rep(NA, times=add.genes))
    gene.list <- as.data.frame(gene.list)
    colnames(gene.list) <- mod.name

    #Combine module lists and rename to mod.name
    david <- david %>%
      dplyr::bind_cols(gene.list)
  }

  ##### Combine results into list #####
  dat.mods <- list()
  if(!is.null(genes)){
    dat.mods[["genes"]] <- genes
  } else {
    dat.mods[["genes"]] <- dat$genes$geneName
  }
  dat.mods[["sft"]] <- sft.select[1,]
  dat.mods[["sft.plot"]] <- sft.plot
  dat.mods[["mods"]] <- mods
  dat.mods[["mods.voom"]] <- mods.voom
  dat.mods[["david"]] <- david
  return(dat.mods)
}
