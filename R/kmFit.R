#' Linear mixed effects models with kinship for RNA-seq
#'
#' Run lmerel and corresponding lm or lme without kinship of gene expression in RNA-seq data
#'
#' @param dat EList object output by voom( ). Must contain counts (dat$E) and meta (dat$targets). Optionally also contains gene metadata (dat$genes) and weights (dat$weights)
#' @param kin Matrix with pairwise kinship values between individuals. Must be numeric with rownames.
#' @param patientID Character of variable name to match dat$targets to kinship row and column names.
#' @param libraryID Character of variable name to match dat$targets to dat$E colnames
#' @param counts Matrix of normalized expression. Rows are genes, columns are libraries.
#' @param meta Matrix or data frame of sample and individual metadata.
#' @param genes Optional matrix or data frame of gene metadata.
#' @param weights Optional matrix of data frame of gene specific weights. Usually calculated with limma::voomWithQualityWeights().
#' @param subset.var Character list of variable name(s) to filter data by.
#' @param subset.lvl Character list of variable value(s) or level(s) to filter data to. Must match order of subset.var
#' @param subset.genes Character vector of genes to include in models.
#' @param model Character vector of model starting with ~ Should include (1|patientID) if mixed effects will be run
#' @param use.weights Logical if gene specific weights should be used in model. Default is FALSE
#' @param run.lm Logical if should run lm model without kinship
#' @param run.lme Logical if should run lme model without kinship
#' @param run.lmerel Logical if should run lmerel model with kinship
#' @param run.contrast Logical if should run pairwise contrasts. If no matrix provided, all possible pairwise comparisons are completed.
#' @param contrast.var Character vector of variable in model to run contrasts of. Interaction terms must be specified as "var1:var2". If NULL (default), all contrasts for all variables in the model are run
#' @param metrics Logical if should calculate model fit metrics such as AIC, BIC, R-squared. Default is FALSE
#' @param processors Numeric processors to run in parallel. Default is 2 less than the total available
#' @param p.method Character of FDR adjustment method. Values as in p.adjust( )
#' @param genotype.name Character string. Used internally for kmFit_eQTL
#' @param run.lmekin Depreciated. Please use run.lmerel
#'
#' @return List of data frames including
#'    - lm/lme/lmerel: model estimates and significance
#'    - *.contrast: model estimates and significance for pairwise contrasts with variables in the original model
#'    - *.fit: model fit metrics such as sigma, AIC, BIC, R-squared (optional with metrics paramater)
#'    - *.error: error messages for genes that failed model fitting
#'
#' @importFrom foreach %dopar%
#' @export
#'
#' @examples
#' # All samples and all genes
#' ## Not run
#' # kmFit(dat = example.voom,
#' #     kin = example.kin, run.lmerel = TRUE,
#' #     model = "~ virus + (1|ptID)")
#'
#' # Subset samples and genes
#' ## Also with weights
#' kmFit(dat = example.voom,
#'       run.lm = TRUE, use.weights = FALSE,
#'       subset.var = list("asthma"), subset.lvl = list(c("asthma")),
#'       subset.genes = c("ENSG00000250479","ENSG00000250510","ENSG00000255823"),
#'       model = "~ virus + (1|ptID)")
#'
#' # Pairwise contrasts
#' ## Continuous interaction
#' kmFit(dat = example.voom,
#'       run.lme = TRUE, run.contrast = TRUE,
#'       subset.genes = c("ENSG00000250479","ENSG00000250510","ENSG00000255823"),
#'       model = "~ virus + asthma * median_cv_coverage + (1|ptID)",
#'       contrast.var=c("asthma:median_cv_coverage"))
#'
#' ## Categorical interaction
#' kmFit(dat = example.voom, kin = example.kin,
#'       run.lmerel = TRUE, run.contrast = TRUE, metrics=TRUE,
#'       subset.genes = c("ENSG00000250479","ENSG00000250510","ENSG00000255823"),
#'       model = "~ virus*asthma + (1|ptID)",
#'       contrast.var=c("virus:asthma"))
#'
#' # Model with failed genes
#' kmFit(dat = example.voom,
#'       kin = example.kin, run.lmerel = TRUE, run.lm = TRUE,
#'       subset.genes = c("ENSG00000250479","ENSG00000250510","ENSG00000255823"),
#'       model = "~ virus*asthma + lib.size + norm.factors + median_cv_coverage + ptID + (1|ptID)")
#'
#' # Non-dat data
#' kmFit(counts = example.voom$E, meta = example.voom$targets,
#'       run.lm = TRUE, use.weights = FALSE,
#'       subset.genes = c("ENSG00000250479","ENSG00000250510","ENSG00000255823"),
#'       model = "~ virus + (1|ptID)")
#'
#' # Three level variable
#' example.voom$targets$lvl <- rep(c("A","B","C"), length(example.voom$targets$libID)/3)
#' kmFit(dat = example.voom,
#'       run.lme= TRUE, run.contrast = TRUE,
#'       subset.genes = c("ENSG00000250479","ENSG00000250510","ENSG00000255823"),
#'       model = "~ lvl + (1|ptID)")
#'

kmFit <- function(dat=NULL, kin=NULL, patientID="ptID", libraryID="libID",
                  counts=NULL, meta=NULL, genes=NULL, weights=NULL,
                  subset.var = NULL, subset.lvl = NULL, subset.genes = NULL,
                  model, use.weights=FALSE,
                  run.lm = FALSE, run.lme = FALSE, run.lmerel = FALSE,
                  metrics = FALSE, run.lmekin = NULL,
                  run.contrast = FALSE, contrast.var = NULL,
                  processors = NULL, p.method = "BH",
                  genotype.name=NULL){

  rowname <- libID <- variable <- statistic <- df <- pval <- group <- gene <- V1 <- V2 <- combo <- term <- p.value <- estimate <- contrast <- contrast.i <- weights.gene <- FDR <- contrast_ref <- contrast_lvl <- NULL

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
      stop("libraryID column not found in dat$targets.")}
  if(!is.null(meta) & !(libraryID %in% colnames(meta))){
      stop("libraryID column not found in meta.")}
  if(!is.null(dat) & !(patientID %in% colnames(dat$targets))){
    stop("patientID column not found in dat$targets.")}
  if(!is.null(meta) & !(patientID %in% colnames(meta))){
    stop("patientID column not found in meta.")}
  if(is.null(subset.var) & !is.null(subset.lvl)){
    stop("Sample subsetting has been selected. Please also provide subset.var")}
  if(!is.null(subset.var) & is.null(subset.lvl)){
    stop("Sample subsetting has been selected. Please also provide subset.lvl")}
  if(run.lmerel & !grepl("\\|", model)){
    stop("Kinship models require a random effect in the model as in (1 | ptID)")}
  if(is.null(kin) & run.lmerel){
    stop("Kinship matrix is required to run lmerel")}
  if(!run.lm & !run.lme & !run.lmerel & !run.contrast){
    stop("At least 1 model type must be selected. Please set one run parameter to TRUE.")}
  if(!run.lm & !run.lme & !run.lmerel & run.contrast){
      stop("Contrast models must be run with an accompanying linear model.")}
  if(use.weights & is.null(weights) & is.null(dat$weights)){
    stop("When use.weights is TRUE, must provide gene weights is dat object or separate data frame.")
  }
  if("gene_weight" %in% c(colnames(meta), colnames(dat$targets))){
    stop("Variable gene_weight is present in dat$targets or meta. This name is used for model weights. Please change variable name in your data.")
  }
  if(!is.null(run.lmekin)){
    stop("run.lmekin no longer supported. Please use run.lmerel")
  }
  if(grepl("\n", model)){
    stop("Model cannot contain hard returns \n. Please correct.")
  }
  if(!use.weights){
    if(!is.null(weights) | !is.null(dat$weights)){
      message("WARNING: To use weights provided in dat$weights or weights, set use.weights = TRUE\n")
    }
  }

  ###### Formulae #####
  #Make formulae. as.formula does not work
  if(grepl("\\|", model)){
    model.temp <- strsplit(gsub(" ", "", model), split = "\\+\\(1")[[1]][1]
    model.lm <- paste("expression", model.temp, sep="")
  } else {
    model.lm <- paste("expression", model, sep="")
    model.lm <- gsub(" ","",model.lm)
  }
  model.lme <- paste("expression", gsub(" ", "", model), sep="")

  #Model message
  if(run.lm){ message(paste("lm model:",model.lm))}
  if(run.lme | run.lmerel){ message(paste("lme/lmerel model:",model.lme))}

  #If no contrast variable set, as all
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

    #Set eQTL variable
    if(!is.null(genotype.name)){
      contrast.var <- gsub("genotype", genotype.name, contrast.var)
    }
  }

  #If contrast variables given, force run contrast model
  if(!is.null(contrast.var)){run.contrast <- TRUE}

  #Check contrast variable is character/factor
  if(run.contrast){
    if(!is.null(dat)){
      meta.temp <- dat$targets
    } else {
      meta.temp <- meta
    }
    for(v in contrast.var){
      #Deal with interaction terms. Only 1 need be non-numeric
      if(grepl("[*]|:", v)){
        v.sep <- strsplit(v, split="[*]|:")[[1]]
        v.class1 <- class(meta.temp[,v.sep[1]])
        v.class2 <- class(meta.temp[,v.sep[2]])
        if(all(c(v.class1,v.class2) %in% c("numeric","integer"))){
          stop(paste("Contrast variable", v, "is numeric. Please specify only character/factor contrasts."))
        }
      } else {
        v.class <- class(meta.temp[,v])
        if(v.class %in% c("numeric","integer")){
          stop(paste("Contrast variable", v, "is numeric. Please specify only character/factor contrasts."))
        }
      }
    }
  }

  ###### Data #####
  to.model.ls <- kimma_cleaning(dat, kin, patientID, libraryID,
                                counts, meta, genes, weights,
                                subset.var, subset.lvl, subset.genes,
                                model.lm, genotype.name)

  ###### Run models ######
  #create blank df to hold results
  fit.results <- data.frame()
  doParallel::registerDoParallel(cl)

  #Loop through each gene
  fit.results <- data.table::rbindlist(fill=TRUE,
                    foreach::foreach(gene=unique(to.model.ls[["to.model"]]$rowname),
                                     .packages = c("dplyr","magrittr","stats","broom","lme4",
                                                   "car","tibble","lme4qtl","utils","emmeans",
                                                   "data.table","foreach","doParallel"),
                                     .export = c("kimma_lm","kimma_lme","kimma_lmerel",
                                                 "kmFit_contrast")) %dopar% {

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
       kimma_lm(model.lm, to.model.gene, gene, use.weights, metrics)
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
        kimma_lme(model.lme, to.model.gene, gene, use.weights, metrics)
        }, error=function(e){
          results.lme.ls[["error"]] <- data.frame(model="lme",
                                                  gene=gene,
                                                  message=conditionMessage(e))
          return(results.lme.ls)
        })
    }

    ##### Kinship model ######
    results.kin.ls <- NULL
    if(run.lmerel){
      #Wrap model run in error catch to allow loop to continue even if a single model fails
      results.kin.ls <- tryCatch({
        kimma_lmerel(model.lme, to.model.gene, gene, to.model.ls[["kin.subset"]],
                     use.weights, patientID, metrics)
        }, error=function(e){
          results.kin.ls[["error"]] <- data.frame(model="lmerel",
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
      if(!is.null(results.lm.ls[["results"]])){
        contrast.lm <- tryCatch({
          kmFit_contrast(results.lm.ls[["fit"]], contrast.var, to.model.gene,
                         genotype.name) %>%
            dplyr::mutate(model="lm.contrast")
        }, error=function(e){
          contrast.lm.error <- data.frame(model="lm.contrast",
                                     gene=gene,
                                     message=conditionMessage(e))
          return(contrast.lm.error)
        })
      }

      if(!is.null(results.lme.ls[["results"]])){
        contrast.lme <- tryCatch({
          kmFit_contrast(results.lme.ls[["fit"]], contrast.var, to.model.gene,
                         genotype.name) %>%
            dplyr::mutate(model="lme.contrast")
        }, error=function(e){
          contrast.lme.error <- data.frame(model="lme.contrast",
                                     gene=gene,
                                     message=conditionMessage(e))
          return(contrast.lme.error)
        })
      }

      if(!is.null(results.kin.ls[["results"]])){
        contrast.kin <- tryCatch({
          kmFit_contrast(results.kin.ls[["fit"]], contrast.var, to.model.gene,
                         genotype.name) %>%
            dplyr::mutate(model="lmerel.contrast")
        }, error=function(e){
          contrast.kin.error <- data.frame(model="lmerel.contrast",
                                           gene=gene,
                                           message=conditionMessage(e))
          return(contrast.kin.error)
        })
      }

      #Combine contrast results
      contrast.results <- dplyr::bind_rows(contrast.lm, contrast.lme, contrast.kin) %>%
        dplyr::mutate(gene=gene)
      }

    #### Combine results #####
    #All models for this gene
    #If any estimate are character (seeContrasts), force for merging
    if(any(is.character(results.lm.ls[["results"]]$estimate),
       is.character(results.lme.ls[["results"]]$estimate),
       is.character(results.kin.ls[["results"]]$estimate),
       is.character(contrast.results$estimate))){
      #Estimate
      results.lm.ls[["results"]]$estimate <- as.character(results.lm.ls[["results"]]$estimate)
      results.lme.ls[["results"]]$estimate <- as.character(results.lme.ls[["results"]]$estimate)
      results.kin.ls[["results"]]$estimate <- as.character(results.kin.ls[["results"]]$estimate)
      contrast.results$estimate <- as.character(contrast.results$estimate)
      #Statistic
      results.lm.ls[["results"]]$statistic <- as.character(results.lm.ls[["results"]]$statistic)
      results.lme.ls[["results"]]$statistic <- as.character(results.lme.ls[["results"]]$statistic)
      results.kin.ls[["results"]]$statistic <- as.character(results.kin.ls[["results"]]$statistic)
      contrast.results$statistic <- as.character(contrast.results$statistic)
      #df
      # results.lm.ls[["results"]]$df <- as.character(results.lm.ls[["results"]]$df)
      # results.lme.ls[["results"]]$df <- as.character(results.lme.ls[["results"]]$df)
      # results.kin.ls[["results"]]$df <- as.character(results.kin.ls[["results"]]$df)
      # contrast.results$df <- as.character(contrast.results$df)
    }


    fit.results <- results.lm.ls[["results"]] %>%
      dplyr::bind_rows(results.lme.ls[["results"]]) %>%
      dplyr::bind_rows(results.kin.ls[["results"]]) %>%
      dplyr::bind_rows(contrast.results) %>%
      dplyr::bind_rows(results.lm.ls[["error"]]) %>%
      dplyr::bind_rows(results.lme.ls[["error"]]) %>%
      dplyr::bind_rows(results.kin.ls[["error"]]) %>%
      dplyr::bind_rows(results.lm.ls[["metrics"]]) %>%
      dplyr::bind_rows(results.lme.ls[["metrics"]]) %>%
      dplyr::bind_rows(results.kin.ls[["metrics"]])
  })
  parallel::stopCluster(cl)

  #### Completion messages ####
  all <- length(unique(to.model.ls[["to.model"]]$rowname))
  message("Complete: ", all, " genes")

  if("message" %in% colnames(fit.results)){
    fail <- fit.results %>%
      tidyr::drop_na(message) %>%
      dplyr::distinct(gene, message) %>% nrow()

    message("Failed: ", fail, " genes. See results[['model_error']]")

    #move message to separate df
    error.results <- fit.results %>%
      dplyr::filter(!is.na(message)) %>%
      dplyr::distinct(model, gene, message)
    fit.results <- fit.results %>%
      dplyr::filter(is.na(message)) %>%
      dplyr::select(-message)

  } else {
    message("Failed: 0 genes")
    error.results <- NULL
  }

  #### Final formatting ####
  kmFit.ls <- list()

  if(nrow(fit.results) > 0 ){
    # Calculate FDR
    if(run.contrast){
      kmFit.results <- fit.results %>%
        #Within model and variable
        dplyr::group_by(model, variable, contrast_ref, contrast_lvl) %>%
        dplyr::mutate(FDR=stats::p.adjust(pval, method=p.method)) %>%
        dplyr::ungroup() %>%
        dplyr::select(model:statistic, contrast_ref, contrast_lvl,
                      estimate, pval, FDR, dplyr::everything())
    }else{
      kmFit.results <- fit.results %>%
        #Within model and variable
        dplyr::group_by(model, variable) %>%
        dplyr::mutate(FDR=stats::p.adjust(pval, method=p.method)) %>%
        dplyr::ungroup() %>%
        dplyr::select(model:statistic, estimate, pval, FDR,
                      dplyr::everything())
    }

    # Add gene info if available
    # REMOVE because uses too much RAM with large models / contrasts
    # if(!is.null(dat$genes)){
    #   genes <- as.data.frame(dat$genes)
    # }
    #
    # if(!is.null(genes)){
    #   #Find matching column
    #   nameID <- which(apply(genes, 2, function(x)
    #     any(grepl(kmFit.results$gene[1], x))))
    #   name <- colnames(genes)[nameID]
    #
    #   kmFit.results.anno <- kmFit.results %>%
    #     dplyr::left_join(genes, by=c("gene"=name))
    # } else{
    #   kmFit.results.anno <- kmFit.results
    # }
    kmFit.results.anno <- kmFit.results

    # Split into list
    for(result.i in unique(kmFit.results.anno$model)){
      result.temp <- dplyr::filter(kmFit.results.anno, model==result.i)
      #Turn estimate numeric if needed
      estimates <- unique(result.temp$estimate)
      estimates <- estimates[!is.na(estimates)]
      if(all(estimates != "seeContrasts") & !is.null(estimates)){
        result.temp <- dplyr::filter(kmFit.results.anno, model==result.i) %>%
          dplyr::mutate(estimate=as.numeric(estimate),
                        statistic=as.numeric(statistic))
      }


      kmFit.ls[[result.i]] <- result.temp %>%
        dplyr::select_if(function(x) any(!is.na(x)))
    }
  }

 # Split error messages into list object
  if(!is.null(error.results)){
    for(result.i in unique(error.results$model)){
      kmFit.ls[[paste(result.i,"error",sep=".")]] <- dplyr::filter(error.results, model==result.i)
    }}

  #### Save ####
  return(kmFit.ls)
}
