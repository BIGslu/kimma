#' eQTL linear mixed effects models with kinship
#'
#' Run lmerel and corresponding lm or lme without kinship of gene expression in RNA-seq data against associated genotypes
#'
#' @param dat.snp Data frame containing numeric (0,1,2) genotypes. Rows are SNP and columns are donors. Option column matching SNP to genes of interest (see within.gene).
#' @param within.gene Logical if should only run paired SNP and genes in rows of dat.snp. Default is FALSE. Does NOT work in combination with subset.genes
#' @param geneID Character variable name to match dat.snp dat with gene identifiers. Values should match rownames in express data (dat$E or counts) Default it "gene"
#' @param genotypeID Character variable name to match dat.snp dat with SNP identifiers. Default it "snpID"
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
#'                          nrow=nSNP)) %>%
#'   #Add SNP and gene names
#'   dplyr::mutate(snpID = paste0("rs",1:nSNP),
#'                 gene = rep(example.voom$genes$geneName[1:2], nSNP/2)) %>%
#'  dplyr::select(snpID, gene, dplyr::everything())
#'  colnames(snp) <- c("snpID", "gene", unique(example.voom$targets$ptID))
#'
#' eQTL_result <- kmFit_eQTL(dat.snp = snp,
#'                           dat = example.voom, kin = example.kin,
#'                           run.lmerel = TRUE, run.contrast = TRUE, metrics=TRUE,
#'                           subset.genes = c("ENSG00000250479","CENSG00000250510"),
#'                           model = "~ virus*genotype + (1|ptID)")
#'
#' eQTL_result2 <- kmFit_eQTL(dat.snp = snp, within.gene = TRUE,
#'                            dat = example.voom, kin = example.kin,
#'                            run.lmerel = TRUE, run.contrast = TRUE, metrics=TRUE,
#'                            model = "~ virus*genotype + (1|ptID)")


kmFit_eQTL <- function(dat.snp = NULL, within.gene = FALSE,
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
  if(!is.null(subset.genes) & within.gene){
    stop("Only one of within.gene and subset.genes can be provided.")}
  if(within.gene & !(geneID %in% colnames(dat.snp))){
    stop("geneID variable must be in dat.snp in order for within.gene to run properly.")}
  if(!(genotypeID %in% colnames(dat.snp))){
    stop("genotypeID variable not in dat.snp. Please set genotypeID to the correct variable name for your data.")}

  num.check <- as.matrix(dat.snp[,!colnames(dat.snp) %in% c(genotypeID, geneID)])
  if(!is.numeric(num.check)){
    stop("Non-numeric genotypes in dat.snp. Please provide only genotypeID, geneID (optional), and numeric 0,1,2 SNP data.")}

  #### Data ####
  #combine genotype and metadata
  snp.format <- dat.snp %>%
    tibble::column_to_rownames(genotypeID) %>%
    dplyr::select_if(is.numeric) %>%
    t() %>% as.data.frame() %>%
    tibble::rownames_to_column(patientID)

  if(!is.null(dat)){
    dat$targets <- dat$targets %>%
      dplyr::left_join(snp.format, by = patientID)
  }

  if(!is.null(meta)){
    meta <- meta %>%
      dplyr::left_join(snp.format, by = patientID)
  }

  #### Model ####
  result.ls <- list()
  for(s in unique(dat.snp[,genotypeID])){
    message(paste0("\nRunning ", s))

    #Make model for SNP of interest
    m <- gsub("genotype", s, model)

    #Match SNP to subset of genes if within gene is selected
    if(within.gene){
      #List genes of interest
      subset.genes <- dat.snp %>%
        dplyr::filter(get(genotypeID) == s) %>%
        dplyr::pull(geneID)
    }

    temp <- kmFit(dat, kin, patientID, libraryID,
                  counts, meta, genes, weights,
                  subset.var, subset.lvl, subset.genes,
                  m, use.weights,
                  run.lm, run.lme, run.lmerel,
                  metrics, run.lmekin,
                  run.contrast, contrast.var,
                  processors, p.method)

    #combine with previous snp outputs
    for(df in unique(c(names(temp), names(result.ls)))){
      result.ls[[df]] <- temp[[df]] %>%
        dplyr::mutate(genotype = s) %>%
        dplyr::select(genotype, dplyr::everything()) %>%
        dplyr::bind_rows(result.ls[[df]])
    }
  }

  return(result.ls)
}
