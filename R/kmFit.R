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
#' @param compare.lm Logical if should run corresponding lm model without kinship
#' @param compare.lme Logical if should run corresponding lme model without kinship
#' @param contrast Logical if should run pairwise contrasts. If no matrix provided, all possible pairwise comparisons are completed.
#' @param contrast.mat Numeric contrast matrix created limma::makeContrasts( )
#' @param processors Numeric processor to run in parallel
#' @param p.method Character of FDR adjustment method. Values as in p.adjust( )
#'
#' @return Dataframe with model fit and significance for each gene
#' @importFrom foreach %dopar%
#' @export
#'
#' @examples
#' # All samples and all genes
#' # Not run
#' # kmFit(dat = dat.voom.example,
#' #       patientID = "donorID", libraryID = "libID",
#' #       kin = kin.example, compare.lme = TRUE,
#' #       model = "~ virus + (1|donorID)")
#'
#' # Subset samples and genes
#' kmFit(dat = dat.voom.example,
#'       patientID = "donorID", libraryID = "libID",
#'       kin = kin.example,
#'       compare.lme = TRUE,
#'       subset.var = list("donorID"), subset.lvl = list(c("donor1","donor2")),
#'       subset.genes = c("ENSG00000250479","ENSG00000250510","ENSG00000255823"),
#'       model = "~ virus + (1|donorID)")
#' # Pairwise contrasts
#' kmFit(dat = dat.voom.example,
#'       patientID = "donorID", libraryID = "libID",
#'       kin = kin.example,
#'       compare.lm = TRUE, contrast = TRUE,
#'       subset.genes = c("ENSG00000250479","ENSG00000250510","ENSG00000255823"),
#'       model = "~ donorID + (1|donorID)")

kmFit <- function(dat=NULL, kin=NULL, patientID="ptID", libraryID="libID",
                  counts=NULL, meta=NULL, genes=NULL,
                  subset.var = NULL, subset.lvl = NULL, subset.genes = NULL,
                  model, compare.lm = FALSE, compare.lme = FALSE,
                  contrast = FALSE, contrast.mat = NULL,
                  processors = 1, p.method = "BH"){

  rowname <- libID <- libraryID <- variable <- pval <- group <- i <- V1 <- V2 <- combo <- term <- p.value <- estimate <- NULL

  #Log start time
  old <- Sys.time()
  `%notin%` <- Negate(`%in%`)

  ###### Parallel ######
  #setup parallel processors
  doParallel::registerDoParallel(processors)

  ###### Check common input parameter errors #####
  if(is.null(subset.var) & !is.null(subset.lvl)){
    stop("Sample subsetting has been selected. Please also provide subset.var")}
  if(!is.null(subset.var) & is.null(subset.lvl)){
    stop("Sample subsetting has been selected. Please also provide subset.lvl")}
  if(!is.null(kin) & !grepl("\\|", model)){
    stop("Kinship models require a random effect in the model as in (1 | ptID)")
  }

  ###### Data #####
  print("Format data")
  to.model.ls <- kimma_cleaning(dat, kin, patientID, libraryID,
                             counts, meta, genes,
                             subset.var, subset.lvl, subset.genes)
  #Make formulae. as.formula does not work
  if(grepl("\\|", model)){
    model.temp <- gsub(" ", "", model)
    model.temp <- strsplit(model.temp, split = "\\+\\(1")[[1]][1]
    model.lm <- paste("expression", model.temp, sep="")
  } else {
    model.lm <- paste("expression", model, sep="")
  }
  model.lme <- paste("expression", model, sep="")

  #If pairwise contrasts requested
  if(contrast | !is.null(contrast.mat)){
    model.temp <- strsplit(model, split = "\\(1")[[1]][1]
    model.temp <- gsub("~| ", "", model.temp)
    contrast.vars <- strsplit(model.temp, split = "\\+")[[1]]
  }

  if(!is.null(contrast.mat)){
    contrast.list <- as.list(as.data.frame(contrast.mat))
  } else {
    contrast.list <- NULL
  }

  ###### Run models ######
  print("Run models")

  #create blank df to hold results
  fit.results <- data.frame()

  #Loop through each gene
  fit.results <- data.table::rbindlist(fill=TRUE,
                    foreach::foreach(i=1:length(unique(to.model.ls[["to.model"]]$rowname))) %dopar% {
    #### Prepare data ####
    #Get gene name
    gene <- unique(to.model.ls[["to.model"]]$rowname)[i]

    #Filter data to gene
    to.model.gene <- to.model.ls[["to.model"]] %>%
      dplyr::filter(rowname == gene) %>%
      dplyr::arrange(patientID)

    #### Simple LM models, if selected #####
    #Run linear model without kinship
    results.lm.ls <- NULL
    if(compare.lm){
    #Wrap model run in error catch to allow loop to continue even if a single model fails
     results.lm.ls <- tryCatch({
       kimma_lm(model.lm, to.model.gene, gene)
     }, error=function(e){})
    }

    #### Simple LME models, if selected #####
    results.lme.ls <- NULL
    if(compare.lme){
      #Wrap model run in error catch to allow loop to continue even if a single model fails
      results.lme.ls <- tryCatch({
        kimma_lme(model.lme, to.model.gene, gene)
        }, error=function(e){})
    }

    ##### Kinship model ######
    results.kin.ls <- NULL
    if(!is.null(kin)){
      #Wrap model run in error catch to allow loop to continue even if a single model fails
      results.kin.ls <- tryCatch({
        kimma_lmekin(model.lme, to.model.gene, gene, to.model.ls[["kin.subset"]])
        }, error=function(e){})
    }

    #### Contrasts ####
    contrast.lm <- NULL
    contrast.lme <- NULL
    contrast.kin <- NULL
    contrast.results <- NULL

    if(contrast | !is.null(contrast.mat)){
      if(!is.null(results.lm.ls)){
        contrast.lm <- kmFit_contrast(results.lm.ls[["fit"]], contrast.vars, contrast.list)%>%
          dplyr::mutate(model="lm.contrast")
      }

      if(!is.null(results.lme.ls)){
        contrast.lme <- kmFit_contrast(results.lme.ls[["fit"]], contrast.vars, contrast.list) %>%
          dplyr::mutate(model="lme.contrast")
      }

      if(!is.null(results.kin.ls)){
        #emmeans does not work for lmekin object. Instead, run pairwise contrasts as subsets
        for(var.i in contrast.vars){
          var.class <- unlist(to.model.gene[,var.i])
          #only run on categorical variables
          if(is.character(var.class) | is.factor(var.class)){
          model.contrast <- paste("expression~",var.i, "+(1|", patientID, ")", sep="")

          #list all possible combos
          contrast.combo <- as.data.frame(t(utils::combn(unlist(
            unique(to.model.gene[,var.i])), m=2))) %>%
            dplyr::mutate(combo = paste(var.i, V2, " - ",var.i, V1, sep=""))

          #if contrast matrix provided, remove unused combos
          if(!is.null(contrast.mat)){
            contrast.combo <- contrast.combo %>%
              dplyr::filter(combo %in% colnames(contrast.mat))
            }


          if(nrow(contrast.combo)>0){
          for(row.i in 1:nrow(contrast.combo)){
            #filter data to each pairwise comparison of interest
            to.model.gene.contrast <- dplyr::filter(to.model.gene,
                                                    get(var.i) %in% unlist(contrast.combo[row.i,]))

            to.keep <- unique(unlist(to.model.gene.contrast[,patientID]))
            kin.contrast <- to.model.ls[["kin.subset"]][rownames(to.model.ls[["kin.subset"]]) %in%
                                                          to.keep, to.keep]

            contrast.kin.temp <- tryCatch({
              kimma_lmekin(model.contrast, to.model.gene.contrast, gene, kin.contrast)
              }, error=function(e){})

            contrast.kin <- contrast.kin.temp[["results"]] %>%
              dplyr::filter(variable != "(Intercept)") %>%
              dplyr::mutate(variable = var.i,
                            contrast = contrast.combo[row.i, "combo"],
                            model = "lmekin.contrast") %>%
              dplyr::bind_rows(contrast.kin)
          }
            }
          }
        }
      }

      #Combine contrast results
      contrast.results <- contrast.lm %>%
        dplyr::bind_rows(contrast.lme) %>%
        dplyr::rename(variable=term, pval=p.value) %>%
        dplyr::mutate(gene=gene) %>%
        dplyr::bind_rows(contrast.kin)  %>%
        dplyr::select(model, gene, variable, contrast, estimate, pval, estimate)
    }

    #### Combine results #####
    #All models for this gene
    results <- results.lm.ls[["results"]] %>%
      dplyr::bind_rows(results.lme.ls[["results"]]) %>%
      dplyr::bind_rows(results.kin.ls[["results"]]) %>%
      dplyr::bind_rows(contrast.results)

    #This gene to all previous gene results
    fit.results <- dplyr::bind_rows(results, fit.results)
  })

  #### Calculate FDR ####
  if(contrast | !is.null(contrast.mat)){
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

  #### Save ####
  return(kmFit.results)

  ###### Fin ######
  print("All models complete")
  #Print total time to run
  new <- Sys.time() - old
  print(new)
}
