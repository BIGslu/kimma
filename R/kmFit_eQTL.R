#' eQTL linear mixed effects models with kinship
#'
#' Run lmerel and corresponding lm or lme without kinship of gene expression in RNA-seq data against associated genotypes
#'
#' @param dat.snp Data frame containing 0,1,2 genotypes. Can be numeric or character/factor depending on desired analysis. Rows are donors and columns are SNP.
#' @param dat.map Data frame mapping geneID to genotypeID. If provided, genotypes will be modeled against only the genes they are paired with in the map. If not provided, genotypes will be modeled against all genes in dat$E / counts
#' @param geneID Character variable name to match dat.snp dat with gene identifiers. Values should match rownames in express data (dat$E or counts) Default it "gene"
#' @param genotypeID Character variable name to match dat.snp dat with SNP identifiers. Default it "snpID"
#' @param model Character vector of model starting with ~ Should include (1|patientID) if mixed effects will be run. NOTE THAT FOR EQTL, the model MUST include genotype as a place holder for the genotype data. For example "~ genotype + (1|ptID)"
#'
#' Same as kmFit( )
#' @param dat EList object output by voom( ). Contains counts (dat$E), meta (dat$targets), and genes (dat$genes).
#' @param kin Matrix with pairwise kinship values between individuals. Must be numeric with rownames.
#' @param patientID Character of variable name to match dat$targets to kinship row and column names.
#' @param libraryID Character of variable name to match dat$targets to dat$E colnames
#' @param counts Matrix of normalized expression. Rows are genes, columns are libraries.
#' @param meta Matrix or data frame of sample and individual metadata.
#' @param genes Matrix or data frame of gene metadata.
#' @param weights Matrix of data frame of gene specific weights. Usually calculated with limma::voomWithQualityWeights().
#' @param subset.var Character list of variable name(s) to filter data by.
#' @param subset.lvl Character list of variable value(s) or level(s) to filter data to. Must match order of subset.var
#' @param subset.genes Character vector of genes to include in models.
#' @param use.weights Logical if gene specific weights should be used in model. Default is FALSE
#' @param run.lm Logical if should run lm model without kinship
#' @param run.lme Logical if should run lme model without kinship
#' @param run.lmerel Logical if should run lmerel model with kinship
#' @param run.contrast Logical if should run pairwise contrasts. If no matrix provided, all possible pairwise comparisons are completed.
#' @param contrast.var Character vector of variable in model to run contrasts of. Interaction terms must be specified as "var1:var2". If NULL (default), all contrasts for all variables in the model are run
#' @param metrics Logical if should calculate model fit metrics such as AIC, BIC, R-squared. Default is FALSE
#' @param processors Numeric processors to run in parallel. Default is 2 less than the total available
#' @param p.method Character of FDR adjustment method. Values as in p.adjust( )
#' @param run.lmekin Depreciated. Please use run.lmerel
#'
#' @return List of data frames including
#'    - lm/lme/lmerel: model estimates and significance
#'    - *.contrast: model estimates and significance for pairwise contrasts with variables in the original model
#'    - *.fit: model fit metrics such as sigma, AIC, BIC, R-squared (optional with metrics paramater)
#'    - *.error: error messages for genes that failed model fitting
#'
#' @export
#'
#' @examples
#' #Determine data size
#' nSNP <- 4
#' nSamp <- length(unique(example.voom$targets$ptID))
#' #Create random genotype data
#' snp <- data.frame(matrix(sample(c(0,1,2), nSNP*nSamp, replace=TRUE),
#'                          nrow=nSamp)) %>%
#'   #Add SNP names
#'   dplyr::mutate(ptID = unique(example.voom$targets$ptID)) %>%
#'   dplyr::select(ptID, dplyr::everything())
#' colnames(snp) <- c("ptID", paste0("rs",1:nSNP))
#'
#' #Create gene to SNP map
#'  map <- data.frame(snpID = paste0("rs",1:nSNP),
#'                  gene = rep(example.voom$genes$geneName[1:2], nSNP/2))
#'
#' #Run SNP interaction model
#' eQTL_result <- kmFit_eQTL(dat.snp = snp, dat.map = map,
#'                           dat = example.voom, kin = example.kin,
#'                           run.lmerel = TRUE, run.contrast = TRUE, metrics=TRUE,
#'                           model = "~ virus*genotype + (1|ptID)")
#'
#' #Run SNP as factors with contrasts
#' snp_fct <- snp %>%
#'  dplyr::mutate_if(is.numeric, ~as.factor(.))
#'
#' eQTL_result2 <- kmFit_eQTL(dat.snp = snp_fct, dat.map = map,
#'                            dat = example.voom, kin = example.kin,
#'                            run.lmerel = TRUE, run.contrast = TRUE, metrics=TRUE,
#'                            model = "~ genotype + (1|ptID)")

kmFit_eQTL <- function(dat.snp = NULL, dat.map = NULL,
                       geneID = "gene", genotypeID = "snpID",
                       #Same as kmFit
                       dat=NULL, kin=NULL, patientID="ptID", libraryID="libID",
                       counts=NULL, meta=NULL, genes=NULL, weights=NULL,
                       subset.var = NULL, subset.lvl = NULL, subset.genes = NULL,
                       model, use.weights=FALSE,
                       run.lm = FALSE, run.lme = FALSE, run.lmerel = FALSE,
                       metrics = FALSE, run.lmekin = NULL,
                       run.contrast = FALSE, contrast.var = NULL,
                       processors = NULL, p.method = "BH"){

  snp <- genotype <- NULL

  #### Checks ####
  if(!is.null(dat.map) & !is.null(subset.genes)){
    stop("Only one of dat.map and subset.genes can be provided.")}
  if(!(genotypeID %in% colnames(dat.map))){
    stop("genotypeID variable not in dat.map. Please set genotypeID to the correct variable name for your data.")}
  if(!(geneID %in% colnames(dat.map))){
    stop("geneID variable not in dat.map Please set geneID to the correct variable name for your data.")}

  num.check <- as.matrix(dat.snp[,colnames(dat.snp) != patientID])

  if((run.contrast | !is.null(contrast.var)) & is.numeric(num.check)){
    message("WARNING: Genotype data are numeric. If you would like to treat 0,1,2 as distinct levels in the contrast model, please change dat.snp to character/factor.")
  }

  #### Data ####
  #Move rownames if present in dat.snp
  if(!patientID %in% colnames(dat.snp)){
    snp.format <- dat.snp %>%
      tibble::rownames_to_column(patientID)
  } else {
    snp.format <- dat.snp
  }


  if(!is.null(dat)){
    dat$targets <- dat$targets %>%
      dplyr::left_join(snp.format, by = patientID)
  }

  if(!is.null(meta)){
    meta <- meta %>%
      dplyr::left_join(snp.format, by = patientID)
  }

  #### Model ####
  all.snp <- colnames(snp.format)[colnames(snp.format) != patientID]
  result.ls <- list()

  for(s in all.snp){
    message(paste0("\nRunning ", s))

    #Make model for SNP of interest
    m <- gsub("genotype", s, model)

    #Match SNP to subset of genes if within gene is selected
    if(!is.null(dat.map)){
      #List genes of interest
      subset.genes <- dat.map %>%
        dplyr::filter(get(genotypeID) == s) %>%
        dplyr::pull(geneID)
    }

    temp <- kmFit(dat, kin, patientID, libraryID,
                  counts, meta, genes, weights,
                  subset.var, subset.lvl, subset.genes,
                  model=m, use.weights,
                  run.lm, run.lme, run.lmerel,
                  metrics, run.lmekin,
                  run.contrast, contrast.var,
                  processors, p.method, genotype.name=s)

    #combine with previous snp outputs
    for(df in unique(c(names(temp), names(result.ls)))){
      #Reconcile "seeContrasts" is a character for estimate
      if(!is.null(result.ls[[df]]) & "estimate" %in% colnames(temp[[df]])){
        if(any(is.character(result.ls[[df]]$estimate), is.character(temp[[df]]$estimate))){
          temp[[df]]$estimate <- as.character(temp[[df]]$estimate)
          result.ls[[df]]$estimate <- as.character(result.ls[[df]]$estimate)
        }}

      result.ls[[df]] <- temp[[df]] %>%
        dplyr::mutate(genotype = s) %>%
        dplyr::select(genotype, dplyr::everything()) %>%
        dplyr::bind_rows(result.ls[[df]])
    }
  }
  return(result.ls)
}
