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
#' @param weights Matrix or data frame of gene specific weights
#'
#' Subset data (optional)
#' @param subset.var Character list of variable name(s) to filter data by.
#' @param subset.lvl Character list of variable value(s) or level(s) to filter data to. Must match order of subset.var
#' @param subset.genes Character vector of genes to include in models.
#' @return Data frame formatted for use in kmFit
#' @keywords internal

kimma_cleaning <- function(dat=NULL, kin=NULL, patientID="ptID", libraryID="libID",
                           counts=NULL, meta=NULL, genes=NULL, weights=NULL,
                           subset.var = NULL, subset.lvl = NULL, subset.genes = NULL){
  i <- rowname <- libID <- NULL
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
    meta.format <- as.data.frame(meta) %>%
      dplyr::filter(get(libraryID) %in% colnames(counts.format))

    #Weights if provided
    if(!is.null(weights)){
      weights.format <- as.data.frame(weights)
      if(rownames(counts)[1]!=1){
        rownames(weights.format) <- rownames(counts)
        colnames(weights.format) <- colnames(counts)
      } else {
        rownames(weights.format) <- as.data.frame(counts) %>%
          dplyr::select_if(is.character) %>% unlist(use.names = FALSE)
        colnames(weights.format) <- colnames(counts)[-1]
      }

      weights.format <- weights.format[order(rownames(weights.format)),
                                       colnames(counts.format)]
    } else{
      weights.format <- NULL
    }

    #Put in list
    dat.format$E <- counts.format
    dat.format$targets <- meta.format
    dat.format$genes <- genes
    dat.format$weights <- weights.format
  } else {
    dat.format <- dat

    #Format weights from voom object
    if(!is.null(dat$weights)){
      weights.format <- as.data.frame(dat$weights)
      if(rownames(dat$E)[1]!=1){
        rownames(weights.format) <- rownames(dat$E)
        colnames(weights.format) <- colnames(dat$E)
      } else {
        rownames(weights.format) <- as.data.frame(dat$E) %>%
          dplyr::select_if(is.character) %>% unlist(use.names = FALSE)
        colnames(weights.format) <- colnames(dat$E)[-1]
      }
      dat.format$weights <- weights.format[order(rownames(weights.format)),
                                            colnames(dat$E)]
    } else {
      dat.format$weights <- NULL
    }


  }

  #Format data
  #If has rownames, move into df
  if(is.numeric(dat.format$E[,1])){
    dat.format$E <- tibble::rownames_to_column(as.data.frame(dat.format$E))
  } else {
    #Rename 1st column
    colnames(dat.format$E)[1] <- "rowname"
  }

  if(!is.null(dat.format$weights) & is.numeric(dat.format$weights[,1])){
    dat.format$weights <- tibble::rownames_to_column(dat.format$weights)
  } else if(!is.null(dat.format$weights)){
    #Rename 1st column
    colnames(dat.format$weights)[1] <- "rowname"
  }

  ###### Subset to variable of interest if selected ######
  dat.subset <- dat.format

  #Subset samples
  if(!is.null(subset.var)){
    for(i in 1:length(subset.var)) {
      dat.subset$targets <- dplyr::filter(dat.subset$targets,
                                          get(subset.var[[i]]) %in% subset.lvl[[i]])

      dat.subset$E <- dplyr::select(as.data.frame(dat.subset$E),
                                    rowname,
                                    tidyselect::all_of(dat.subset$targets[,libraryID]))

      if(!is.null(dat.subset$weights)){
        dat.subset$weights <- dplyr::select(as.data.frame(dat.subset$weights),
                                            rowname,
                                            tidyselect::all_of(dat.subset$targets[,libraryID]))
      }
    }
  }

  #Subset genes
  if(!is.null(subset.genes)){
    dat.subset$E <- dplyr::filter(as.data.frame(dat.subset$E),
                                  rowname %in% subset.genes)
    if(!is.null(dat.subset$weights)){
      dat.subset$weights <- dplyr::filter(as.data.frame(dat.subset$weights),
                                          rowname %in% subset.genes)
    }
  }

  ###### Format data for modeling ####
  # Convert IDs if they are the same variable
  if(patientID == libraryID){ to.modelID <- "libID" } else{ to.modelID <- patientID }

  if(!is.null(kin)){
    #Combine expression data (E) and sample metadata (targets)
    to.model <- dat.subset$E %>%
      tidyr::pivot_longer(-rowname, names_to = "libID", values_to = "expression") %>%
      dplyr::inner_join(dat.subset$targets, by=c("libID"=libraryID)) %>%
      #Remove samples missing kinship
      dplyr::filter(get(to.modelID) %in% colnames(kin)) %>%
      dplyr::arrange(get(to.modelID))

    #Add weights if available
    if(!is.null(dat.subset$weights)){
      to.model <- to.model %>%
        dplyr::left_join(tidyr::pivot_longer(dat.subset$weights, -rowname,
                                             names_to = "libID", values_to = "weight"),
                         by=c("rowname", "libID"))
    } else{
      to.model <- to.model %>%
        dplyr::mutate(weight = NA)
    }

    #Remove samples from kinship missing expression data
    #Order kinship as in to.model
    to.keep <- unique(unlist(to.model[,to.modelID]))
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
      dplyr::distinct(get(to.modelID)) %>% nrow()

    message(paste("Running models on", kin.no, "individuals.",
                  rna.no-kin.no, "individuals missing kinship data."))
  }else{
    #Combine expression data (E) and sample metadata (targets)
    to.model <- dat.subset$E %>%
      tidyr::pivot_longer(-rowname, names_to = "libID", values_to = "expression") %>%
      dplyr::inner_join(dat.subset$targets, by=c("libID"=libraryID))

    #Add weights if available
    if(!is.null(dat.subset$weights)){
      to.model <- to.model %>%
        dplyr::left_join(tidyr::pivot_longer(dat.subset$weights, -rowname,
                                             names_to = "libID", values_to = "weight"),
                         by=c("rowname", "libID"))
    } else{
      to.model <- to.model %>%
        dplyr::mutate(weight = NA)
    }

    kin.subset <- NULL

    rna.no <- to.model %>%
      dplyr::distinct(get(patientID)) %>% nrow()

    message(paste("Running models on", rna.no, "individuals. No kinship provided."))
  }

  #Put back ptID if was renamed as libID
  if(!(patientID %in% colnames(to.model))){ to.model[[patientID]] <- to.model$libID }

  #Combine for saving
  to.model.ls <- list()
  to.model.ls[["to.model"]] <- to.model
  to.model.ls[["kin.subset"]] <- kin.subset
  to.model.ls[["dat.subset"]] <- dat.subset

  return(to.model.ls)
}
