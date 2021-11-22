#' kmFit data cleaning
#'
#' Data
#' @param dat EList object output by voom( ). Contains counts (dat$E), meta (dat$targets), and genes (dat$genes).
#' @param kin Matrix with pairwise kinship values between individuals. Must be numeric with rownames.
#' @param patientID Character of variable name to match dat$targets to kinship row and column names.
#' @param libraryID Character of variable name to match dat$targets to dat$E colnames
#'
#' Alternate data if not using EList object
#' @param counts Matrix of normalized expression. Rows are genes, columns are libraries.
#' @param meta Matrix or data frame of sample and individual metadata.
#' @param genes Matrix or data frame of gene metadata.
#'
#' Subset data (optional)
#' @param subset.var Character list of variable name(s) to filter data by.
#' @param subset.lvl Character list of variable value(s) or level(s) to filter data to. Must match order of subset.var
#' @param subset.genes Character vector of genes to include in models.
#' @return Data frame formatted for use in kmFit
#' @keywords internal

kimma_cleaning <- function(dat=NULL, kin=NULL, patientID="ptID", libraryID="libID",
                           counts=NULL, meta=NULL, genes=NULL,
                           subset.var = NULL, subset.lvl = NULL, subset.genes = NULL){
  i <- rowname <- NULL
  #If data are NOT a voom EList, create a mock version
  if(is.null(dat)) {
    dat.format <- list()

    #Expression data
    ##Move rownames to column if exist
    ##Order columns as in metadata and genes as in genes
    if(rownames(counts)[1]!=1){
      counts.format <- as.data.frame(counts) %>%
        tibble::rownames_to_column() %>%
        dplyr::select(rowname, tidyselect::all_of(meta[,libraryID])) %>%
        dplyr::arrange(match(rowname, genes$geneName)) %>%
        tibble::column_to_rownames()
    } else {
      counts.format <- as.data.frame(counts) %>%
        dplyr::rename_if(is.character, ~"rowname")%>%
        dplyr::select(rowname, tidyselect::all_of(meta[,libraryID])) %>%
        dplyr::arrange(match(rowname, genes$geneName)) %>%
        tibble::column_to_rownames()
    }

    #Metadata
    ##Remove samples not in expression data
    meta.format <- meta %>%
      dplyr::filter(get(libraryID) %in% colnames(counts.format))

    #Put in list
    dat.format$E <- counts.format
    dat.format$targets <- meta.format
    dat.format$genes <- genes
  } else {
    dat.format <- dat
  }

  #Format data
  #If has rownames, move into df
  if(is.numeric(dat.format$E[,1])){
    dat.format$E <- tibble::rownames_to_column(as.data.frame(dat.format$E))
  } else {
    #Rename 1st column
    colnames(dat.format$E)[1] <- "rowname"
  }

  ###### Subset to variable of interest if selected ######
  dat.subset <- dat.format

  #Subset samples
  if(!is.null(subset.var)){
    for(i in 1:length(subset.var)) {
      dat.subset$targets <- dplyr::filter(dat.subset$targets,
                                          get(subset.var[[i]]) %in% subset.lvl[[i]])

      dat.subset$E <- dplyr::select(as.data.frame(dat.subset$E),rowname,
                                    tidyselect::all_of(dat.subset$targets$libID))
    }
  }

  #Subset genes
  if(!is.null(subset.genes)){
    dat.subset$E <- dplyr::filter(as.data.frame(dat.subset$E), rowname %in% subset.genes)
  }

  ###### Format data for modeling ####
  if(!is.null(kin)){
    #Combine expression data (E) and sample metadata (targets)
    to.model <- dat.subset$E %>%
      tidyr::pivot_longer(-rowname, names_to = "libID", values_to = "expression") %>%
      dplyr::inner_join(dat.subset$targets, by=c("libID"=libraryID)) %>%
      #Remove samples missing kinship
      dplyr::filter(get(patientID) %in% colnames(kin)) %>%
      dplyr::arrange(get(patientID))

    #Remove samples from kinship missing expression data
    #Order kinship as in to.model
    to.keep <- unique(unlist(to.model[,patientID]))
    kin.subset <- as.data.frame(kin) %>%
      tibble::rownames_to_column() %>%
      dplyr::filter(rowname %in% to.keep) %>%
      dplyr::select(rowname, tidyselect::all_of(to.keep)) %>%
      dplyr::arrange(rowname) %>%
      tibble::column_to_rownames()

    #Compute number of samples to run in models
    rna.no <- dat.subset$targets %>%
      dplyr::distinct(get(patientID)) %>% nrow()
    kin.no <- to.model %>%
      dplyr::distinct(get(patientID)) %>% nrow()

    message(paste("Running models on", kin.no, "individuals.",
                  rna.no-kin.no, "individuals missing kinship data."))
  }else{
    #Combine expression data (E) and sample metadata (targets)
    to.model <- dat.subset$E %>%
      tidyr::pivot_longer(-rowname, names_to = "libID", values_to = "expression") %>%
      dplyr::inner_join(dat.subset$targets, by=c("libID"=libraryID))

    kin.subset <- NULL
    #Compute number of samples to run in models
    if(patientID == libraryID){
      rna.no <- to.model %>%
        dplyr::distinct(libID) %>% nrow()
    } else {
      rna.no <- to.model %>%
        dplyr::distinct(get(patientID)) %>% nrow()
    }

    message(paste("Running models on", rna.no, "individuals. No kinship provided."))
  }

  #Combine for saving
  to.model.ls <- list()
  to.model.ls[["to.model"]] <- to.model
  to.model.ls[["kin.subset"]] <- kin.subset
  to.model.ls[["dat.subset"]] <- dat.subset

  return(to.model.ls)
}
