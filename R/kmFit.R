#' Linear mixed effects models with kinship for RNA-seq
#'
#' Run lmekin and corresponding lm or lme without kinship of gene expression in RNA-seq data
#'
#' @param dat EList object output by voom( ). Contains counts (dat$E), meta (dat$targets), and genes (dat$genes).
#' @param kin Matrix with pairwise kinship values between individuals. Must be numeric with rownames.
#' @param patientID Character of variable name to match dat$targets to kinship row and column names.
#' @param libraryID Character of variable name to match dat$targets to dat$E colnames
#' @param counts Matrix of normalized expression. Rows are genes, columns are libraries.
#' @param meta Matrix or data frame of sample and individual metadata.
#' @param genes Matrix or data frame of gene metadata.
#' @param subset.var Character list of variable name(s) to filter data by.
#' @param subset.lvl Character list of variable value(s) or level(s) to filter data to. Must match order of subset.var
#' @param subset.genes Character vector of genes to include in models.
#' @param model Character vector of model starting with ~ Should include (1|patientID) if mixed effects will be run
#' @param run.lm Logical if should run lm model without kinship
#' @param run.lme Logical if should run lme model without kinship
#' @param run.lmekin Logical if should run lmekin model with kinship
#' @param run.contrast Logical if should run pairwise contrasts. If no matrix provided, all possible pairwise comparisons are completed.
#' @param contrast.var Character vector of variable in model to run contrasts of. Interaction terms must be specified as "var1:var2". If NULL (default), all contrasts for all variables in the model are run
#' @param processors Numeric processors to run in parallel. Default is 2 less than the total available
#' @param p.method Character of FDR adjustment method. Values as in p.adjust( )
#'
#' @return Dataframe with model fit and significance for each gene
#' @importFrom foreach %dopar%
#' @export
#'
#' @examples
#' # All samples and all genes
#' # Not run
#' # kmFit(dat = example.voom,
#' #       patientID = "donorID", libraryID = "libID",
#' #       kin = example.kin, run.lmekin = TRUE,
#' #       model = "~ virus + (1|donorID)")
#'
#' # Subset samples and genes
#' kmFit(dat = example.voom,
#'       patientID = "donorID", libraryID = "libID",
#'       run.lme = TRUE,
#'       subset.var = list("asthma"), subset.lvl = list(c("asthma")),
#'       subset.genes = c("ENSG00000250479","ENSG00000250510","ENSG00000255823"),
#'       model = "~ virus + (1|donorID)")
#' # Pairwise contrasts
#' kmFit(dat = example.voom,
#'       patientID = "donorID", libraryID = "libID",
#'       run.lme = TRUE, run.contrast = TRUE,
#'       subset.genes = c("ENSG00000250479","ENSG00000250510","ENSG00000255823"),
#'       model = "~ virus+asthma * median_cv_coverage + (1|donorID)",
#'       contrast.var=c("virus","asthma:median_cv_coverage"))
#'
#' kmFit(dat = example.voom, kin = example.kin,
#'       patientID = "donorID", libraryID = "libID",
#'       run.lmekin = TRUE, run.contrast = TRUE,
#'       subset.genes = c("ENSG00000250479","ENSG00000250510","ENSG00000255823"),
#'       model = "~ virus*asthma + (1|donorID)",
#'       contrast.var=c("virus","virus:asthma"))
#'
#' # Model with failed genes
#' kmFit(dat = example.voom,
#'       patientID = "donorID", libraryID = "libID",
#'       kin = example.kin, run.lmekin = TRUE, run.lm = TRUE,
#'       subset.genes = c("ENSG00000250479","ENSG00000250510","ENSG00000255823"),
#'       model = "~ virus*asthma + lib.size + norm.factors + median_cv_coverage + donorID+(1|donorID)")

kmFit <- function(dat=NULL, kin=NULL, patientID="ptID", libraryID="libID",
                  counts=NULL, meta=NULL, genes=NULL,
                  subset.var = NULL, subset.lvl = NULL, subset.genes = NULL,
                  model, run.lm = FALSE, run.lme = FALSE, run.lmekin = FALSE,
                  run.contrast = FALSE, contrast.var = NULL,
                  processors = NULL, p.method = "BH"){

  rowname <- libID <- variable <- pval <- group <- gene <- V1 <- V2 <- combo <- term <- p.value <- estimate <- contrast <- contrast.i <- NULL

  ###### Parallel ######
  #setup parallel processors
  chk <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")

  if (nzchar(chk) && chk == "TRUE") {
    #Use 2 in CRAN/Travis/AppVeyor
    processors.to.use <- 2
  } else if (is.null(processors)){
    #Use 2 less than total if not user defined
    processors.to.use <- parallel::detectCores()-2
    if(processors.to.use == 0){
      stop("Error processors: Default resulted in 0. Please correct.")}
  } else {
    #Use user defined number
    processors.to.use <- processors
  }

  cl <- parallel::makeCluster(processors.to.use)

  ###### Check common input parameter errors #####
  if(!is.null(dat) & !(libraryID %in% colnames(dat$targets))){
      stop("LibraryID column not found in dat$targets.")}
  if(!is.null(meta) & !(libraryID %in% colnames(meta))){
      stop("LibraryID column not found in meta.")}
  if(is.null(subset.var) & !is.null(subset.lvl)){
    stop("Sample subsetting has been selected. Please also provide subset.var")}
  if(!is.null(subset.var) & is.null(subset.lvl)){
    stop("Sample subsetting has been selected. Please also provide subset.lvl")}
  if(run.lmekin & !grepl("\\|", model)){
    stop("Kinship models require a random effect in the model as in (1 | ptID)")}
  if(is.null(kin) & run.lmekin){
    stop("Kinship matrix is required to run lmekin")}
  if(!run.lm & !run.lme & !run.lmekin & !run.contrast){
    stop("At least 1 model type must be selected. Please set one run parameter to TRUE.")}
  if(!run.lm & !run.lme & !run.lmekin & run.contrast){
      stop("Contrast models must be run with an accompanying linear model.")}

  ###### Data #####
  #print("Format data")
  to.model.ls <- kimma_cleaning(dat, kin, patientID, libraryID,
                             counts, meta, genes,
                             subset.var, subset.lvl, subset.genes)

  #Make formulae. as.formula does not work
  if(grepl("\\|", model)){
    model.temp <- strsplit(gsub(" ", "", model), split = "\\+\\(1")[[1]][1]
    model.lm <- paste("expression", model.temp, sep="")
  } else {
    model.lm <- paste("expression", model, sep="")
  }
  model.lme <- paste("expression", gsub(" ", "", model), sep="")

  #Model message
  if(run.lm){ message(paste("lm model:",model.lm))}
  if(run.lme | run.lmekin){ message(paste("lme/lmekin model:",model.lme))}

  #If no contrast variable set, us all
  if(run.contrast & is.null(contrast.var)){
    contrast.temp <- strsplit(gsub(" |~", "", model), split = "\\+\\(1")[[1]][1]
    #main term
    contrast.main <- strsplit(contrast.temp, split = "\\+|\\*")[[1]]
    #interaction term
    if(grepl("\\*",model)){
      contrast.interact <- strsplit(contrast.temp, split = "\\+")[[1]]
      contrast.interact <- contrast.interact[grepl("\\*",contrast.interact)]
      #Check for triple interactions
      if(any(stringr::str_count(contrast.interact, "\\*") > 1)){
        stop("Contrasts of triple interactions are not supported.")}
      contrast.interact <- gsub("\\*",":", contrast.interact)
    } else {
      contrast.interact <- NULL
    }

    contrast.var <- c(contrast.main, contrast.interact)
  }

  #If contrast variables given, force run contrast model
  if(!is.null(contrast.var)){run.contrast <- TRUE}

  ###### Run models ######
  #create blank df to hold results
  fit.results <- data.frame()
  doParallel::registerDoParallel(cl)

  #Loop through each gene
  fit.results <- data.table::rbindlist(fill=TRUE,
                    foreach::foreach(gene=unique(to.model.ls[["to.model"]]$rowname),
                                     .packages = c("dplyr","magrittr","stats","broom","lme4",
                                                   "car","tibble","coxme","utils","emmeans",
                                                   "data.table","foreach","doParallel"),
                                     .export = c("kimma_lm","kimma_lme","kimma_lmekin",
                                                 "kmFit_contrast","kmFit_contrast_kin")) %dopar% {
    #### Prepare data ####
    #Filter data to gene
    to.model.gene <- to.model.ls[["to.model"]] %>%
      dplyr::filter(rowname == gene) %>%
      dplyr::arrange(patientID)

    #### LM model #####
    #Run linear model without kinship
    results.lm.ls <- NULL

    if(run.lm){
    #Wrap model run in error catch to allow loop to continue even if a single model fails
     results.lm.ls <- tryCatch({
       kimma_lm(model.lm, to.model.gene, gene)
     }, error=function(e){
       results.lm.ls[["error"]] <- data.frame(model="lm",
                                               gene=gene,
                                               message=conditionMessage(e))
       return(results.lm.ls)
     })
    }

    #### LME model #####
    results.lme.ls <- NULL
    if(run.lme){
      #Wrap model run in error catch to allow loop to continue even if a single model fails
      results.lme.ls <- tryCatch({
        kimma_lme(model.lme, to.model.gene, gene)
        }, error=function(e){
          results.lme.ls[["error"]] <- data.frame(model="lme",
                                                  gene=gene,
                                                  message=conditionMessage(e))
          return(results.lme.ls)
        })
    }

    ##### Kinship model ######
    results.kin.ls <- NULL
    if(run.lmekin){
      #Wrap model run in error catch to allow loop to continue even if a single model fails
      results.kin.ls <- tryCatch({
        kimma_lmekin(model.lme, to.model.gene, gene, to.model.ls[["kin.subset"]])
        }, error=function(e){
          results.kin.ls[["error"]] <- data.frame(model="lmekin",
                                                  gene=gene,
                                                  message=conditionMessage(e))
          return(results.kin.ls)
        })
    }

    #### Contrasts ####
    contrast.lm <- NULL
    contrast.lme <- NULL
    contrast.kin <- NULL
    contrast.results <- NULL

    if(run.contrast){
      if(!is.null(results.lm.ls)){
        contrast.lm <- kmFit_contrast(results.lm.ls[["fit"]], contrast.var, to.model.gene)%>%
          dplyr::mutate(model="lm.contrast")
      }

      if(!is.null(results.lme.ls)){
        contrast.lme <- kmFit_contrast(results.lme.ls[["fit"]], contrast.var, to.model.gene) %>%
          dplyr::mutate(model="lme.contrast")
      }

      if(!is.null(results.kin.ls)){
        contrast.kin <- kmFit_contrast_kin(contrast.var, to.model.gene,
                                           patientID, to.model.ls, gene) %>%
          dplyr::mutate(model="lmekin.contrast")
      }
      #Combine contrast results
      contrast.results <- dplyr::bind_rows(contrast.lm, contrast.lme, contrast.kin) %>%
        dplyr::mutate(gene=gene) %>%
        dplyr::select(model, gene, variable, contrast, estimate, pval, estimate)
      }

    #### Combine results #####
    #All models for this gene
    #If any estimate are character (seeContrasts), force for merging
    if(is.character(results.lm.ls[["results"]]$estimate) |
       is.character(results.lme.ls[["results"]]$estimate) |
       is.character(results.kin.ls[["results"]]$estimate) |
       is.character(contrast.results$estimate)){
      results.lm.ls[["results"]]$estimate <- as.character(results.lm.ls[["results"]]$estimate)
      results.lme.ls[["results"]]$estimate <- as.character(results.lme.ls[["results"]]$estimate)
      results.kin.ls[["results"]]$estimate <- as.character(results.kin.ls[["results"]]$estimate)
      contrast.results$estimate <- as.character(contrast.results$estimate)

    }

    fit.results <- results.lm.ls[["results"]] %>%
      dplyr::bind_rows(results.lme.ls[["results"]]) %>%
      dplyr::bind_rows(results.kin.ls[["results"]]) %>%
      dplyr::bind_rows(contrast.results) %>%
      dplyr::bind_rows(results.lm.ls[["error"]]) %>%
      dplyr::bind_rows(results.lme.ls[["error"]]) %>%
      dplyr::bind_rows(results.kin.ls[["error"]])
  })
  parallel::stopCluster(cl)

  #### Completion messages ####
  all <- length(unique(to.model.ls[["to.model"]]$rowname))
  message(paste(all, "genes complete."))

  if("message" %in% colnames(fit.results)){
    fail <- fit.results %>%
      tidyr::drop_na(message) %>%
      dplyr::distinct(gene, message) %>% nrow()

    message(paste(fail, "genes failed one or more models. See results[['model_error']]"))

    #move message to separate df
    error.results <- fit.results %>%
      dplyr::filter(!is.na(message)) %>%
      dplyr::distinct(model, gene, message)
    fit.results <- fit.results %>%
      dplyr::filter(is.na(message)) %>%
      dplyr::select(-message)

  } else {
    message("No genes failed.")
    error.results <- NULL
  }

  #### Final formatting ####
  kmFit.ls <- list()

  if(nrow(fit.results) > 0 ){
    # Calculate FDR
    if(run.contrast){
      kmFit.results <- fit.results %>%
        #Within model and variable
        dplyr::group_by(model, variable, contrast) %>%
        dplyr::mutate(FDR=stats::p.adjust(pval, method=p.method)) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(contrast = gsub("contrast","",contrast))
    }else{
      kmFit.results <- fit.results %>%
        #Within model and variable
        dplyr::group_by(model, variable) %>%
        dplyr::mutate(FDR=stats::p.adjust(pval, method=p.method)) %>%
        dplyr::ungroup()
    }
    # Split into list
    for(result.i in unique(kmFit.results$model)){
      kmFit.ls[[result.i]] <- dplyr::filter(kmFit.results, model==result.i)
      #Turn estimate numeric if needed
      estimates <- unique(kmFit.ls[[result.i]]$estimate)
      estimates <- estimates[!is.na(estimates)]
      if(all(estimates != "seeContrasts")){
        kmFit.ls[[result.i]] <- dplyr::filter(kmFit.results, model==result.i) %>%
          dplyr::mutate(estimate=as.numeric(estimate))
      }
    }
  }

 # Split error messages into list object
  if(!is.null(error.results)){
    for(result.i in unique(error.results$model)){
      kmFit.ls[[paste(result.i,"error",sep="_")]] <- dplyr::filter(error.results, model==result.i)
    }}

  #### Save ####
  return(kmFit.ls)
}
