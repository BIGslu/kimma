#' eQTL linear mixed effects models with kinship
#'
#' Run lmerel and corresponding lm or lme without kinship of gene expression in RNA-seq data against associated genotypes
#'
#' @param dat_snp Data frame containing 0,1,2 genotypes. Can be numeric or character/factor depending on desired analysis. If running contrasdts, must be character/factor for pairwise comparisons to run. Rows are donors and columns are SNP.
#' @param dat_map Data frame mapping geneID to genotypeID. If provided, genotypes will be modeled against only the genes they are paired with in the map. If not provided, genotypes will be modeled against all genes in dat$E / counts
#' @param geneID Character variable name to match dat_snp dat with gene identifiers. Values should match rownames in express data (dat$E or counts) Default it "gene"
#' @param genotypeID Character variable name to match dat_snp dat with SNP identifiers. Default it "snpID"
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
#' @param subset_var Character list of variable name(s) to filter data by.
#' @param subset_lvl Character list of variable value(s) or level(s) to filter data to. Must match order of subset.var
#' @param subset_genes Character vector of genes to include in models.
#' @param use_weights Logical if gene specific weights should be used in model. Default is FALSE
#' @param run_lm Logical if should run lm model without kinship
#' @param run_lme Logical if should run lme model without kinship
#' @param run_lmerel Logical if should run lmerel model with kinship
#' @param run_contrast Logical if should run pairwise contrasts. If no matrix provided, all possible pairwise comparisons are completed.
#' @param contrast_var Character vector of variable in model to run contrasts of. Interaction terms must be specified as "var1:var2". If NULL (default), all contrasts for all variables in the model are run
#' @param metrics Logical if should calculate model fit metrics such as AIC, BIC, R-squared. Default is FALSE
#' @param processors Numeric processors to run in parallel. Default is 2 less than the total available
#' @param p_method Character of FDR adjustment method. Values as in p.adjust( )
#'
#' @param dat.snp Deprecated form of dat_snp
#' @param dat.map Deprecated form of dat.map
#' @param run.lmekin Deprecated. Please use run_lmerel
#' @param subset.var Deprecated form of subset_var
#' @param subset.lvl Deprecated form of subset_lvl
#' @param subset.genes Deprecated form of subset_genes
#' @param use.weights Deprecated form of use_weights
#' @param run.lm Deprecated form of run_lm
#' @param run.lme Deprecated form of run_lme
#' @param run.lmerel Deprecated form of run_lmerel
#' @param run.contrast Deprecated form of run_contrast
#' @param contrast.var Deprecated form of contrast_var
#' @param p.method Deprecated form of p_method
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
#' eQTL_result <- kmFit_eQTL(dat_snp = snp, dat_map = map,
#'                           dat = example.voom, kin = example.kin,
#'                           run_lmerel = TRUE, metrics=TRUE,
#'                           run_contrast = TRUE, contrast_var="virus:genotype",
#'                           model = "~ virus*genotype + (1|ptID)")
#'
#' #Run SNP as factors with contrasts
#' snp_fct <- snp %>%
#'  dplyr::mutate_if(is.numeric, ~as.factor(.))
#'
#' eQTL_result2 <- kmFit_eQTL(dat_snp = snp_fct, dat_map = map,
#'                            dat = example.voom, kin = example.kin,
#'                            run_lmerel = TRUE, run_contrast = TRUE, metrics=TRUE,
#'                            model = "~ genotype + (1|ptID)")

kmFit_eQTL <- function(dat_snp = NULL, dat_map = NULL,
                       geneID = "gene", genotypeID = "snpID",
                       #Same as kmFit
                       dat=NULL, kin=NULL, patientID="ptID", libraryID="libID",
                       counts=NULL, meta=NULL, genes=NULL, weights=NULL,
                       subset_var = NULL, subset_lvl = NULL, subset_genes = NULL,
                       model, use_weights=FALSE,
                       run_lm = FALSE, run_lme = FALSE, run_lmerel = FALSE,
                       metrics = FALSE,
                       run_contrast = FALSE, contrast_var = NULL,
                       processors = NULL, p_method = "BH",
                       #deprecated
                       dat.snp = NULL, dat.map = NULL,
                       run.lmekin = NULL,
                       subset.var = NULL, subset.lvl = NULL, subset.genes = NULL,
                       use.weights=FALSE,
                       run.lm = FALSE, run.lme = FALSE, run.lmerel = FALSE,
                       run.contrast = FALSE, contrast.var = NULL,
                       p.method = NULL){

  snp <- genotype <- NULL
  if(any(!is.null(subset.var),!is.null(subset.lvl),
         !is.null(subset.genes),use.weights,run.lm,run.lme,
         run.lmerel,run.contrast,!is.null(contrast.var),
         !is.null(p.method),!is.null(dat.snp), !is.null(dat.map))){
    message("WARNING: Arguments with '.' have been deprecated. Please use '_' versions.")
  }
  if(!is.null(dat.snp)){dat_snp <- dat.snp}
  if(!is.null(dat.map)){dat_map <- dat.map}
  if(!is.null(subset.var)){subset_var <- subset.var}
  if(!is.null(subset.lvl)){subset_lvl <- subset.lvl}
  if(!is.null(subset.genes)){subset_genes <- subset.genes}
  if(use.weights){use_weights <- use.weights}
  if(run.lm){run_lm <- run.lm}
  if(run.lme){run_lme <- run.lme }
  if(run.lmerel){run_lmerel <- run.lmerel}
  if(run.contrast){run_contrast <- run.contrast}
  if(!is.null(contrast.var)){contrast_var <- contrast.var}
  if(!is.null(p.method)){p_method <- p.method}

  #### Checks ####
  if(!is.null(dat_map) & !is.null(subset_genes)){
    stop("Only one of dat_map and subset.genes can be provided.")}
  if(!is.null(dat_map)){
    if(!(genotypeID %in% colnames(dat_map))){
      stop("genotypeID variable not in dat_map. Please set genotypeID to the correct variable name for your data.")}
    if(!(geneID %in% colnames(dat_map))){
      stop("geneID variable not in dat_map Please set geneID to the correct variable name for your data.")}}

  num.check <- as.matrix(dat_snp[,colnames(dat_snp) != patientID])

  if((run_contrast | !is.null(contrast_var)) & is.numeric(num.check)){
    message("WARNING: Genotype data are numeric. If you would like to treat 0,1,2 as distinct levels in the contrast model, please change dat_snp to character/factor.")
  }

  #### Data ####
  #Move rownames if present in dat_snp
  if(!patientID %in% colnames(dat_snp)){
    snp.format <- dat_snp %>%
      tibble::rownames_to_column(patientID)
  } else {
    snp.format <- dat_snp
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
    if(!is.null(dat_map)){
      #List genes of interest
      subset_genes <- dat_map %>%
        dplyr::filter(get(genotypeID) == s) %>%
        dplyr::pull(geneID)
    }

    #Fix contrast to specify SNP
    if(!is.null(contrast_var)){
      contrast_var2 <- gsub("genotype", s, contrast_var)
    } else{
      contrast_var2 <- contrast_var
    }

    temp <- kmFit(dat, kin, patientID, libraryID,
                  counts, meta, genes, weights,
                  subset_var, subset_lvl, subset_genes,
                  model=m, use_weights,
                  run_lm, run_lme, run_lmerel,
                  metrics,
                  run_contrast, contrast_var=contrast_var2,
                  processors, p_method, genotype_name=s,
                  #deprecated
                  run.lmekin,
                  subset.var, subset.lvl, subset.genes,
                  use.weights,
                  run.lm, run.lme, run.lmerel,
                  run.contrast, contrast.var,
                  p.method)


    #combine with previous snp outputs
    for(datf in unique(c(names(temp), names(result.ls)))){
      #Reconcile "seeContrasts" is a character for estimate
      if(!is.null(result.ls[[datf]]) & "estimate" %in% colnames(temp[[datf]])){
        if(any(is.character(result.ls[[datf]]$estimate),
               is.character(temp[[datf]]$estimate))){
          temp[[datf]]$estimate <- as.character(temp[[datf]]$estimate)
          result.ls[[datf]]$estimate <- as.character(result.ls[[datf]]$estimate)

          temp[[datf]]$statistic <- as.character(temp[[datf]]$statistic)
          result.ls[[datf]]$statistic <- as.character(result.ls[[datf]]$statistic)
        }}

      result.ls[[datf]] <- temp[[datf]] %>%
        dplyr::mutate(genotype = s) %>%
        dplyr::select(genotype, dplyr::everything()) %>%
        dplyr::bind_rows(result.ls[[datf]])
    }
  }
  return(result.ls)
}
