filter_rare <- function(dat, min.CPM, gene.var="geneName",
                        min.sample=NULL, min.pct=NULL,
                        name="dat.filter"){
  ##### Packages #####
  `%notin%` <- Negate(`%in%`)

  ##### Check parameters #####
  #Correct input object type?
  if(class(dat) != "DGEList"){ stop("dat object must be a DGEList object") }
  #Set at least one min sample param
  if(is.null(min.sample) & is.null(min.pct)){ stop("Please provide one of min.sample or min.pct") }
  #min samples or percent only?
  if(!is.null(min.sample) & !is.null(min.pct)){ stop("Please provide only one of min.sample or min.pct") }
  #percent as non-decimal?
  if(!is.null(min.pct) & min.pct<1){ warning("min.pct should be percentage between 0 and 100. Your value appears to be a proporion less than 1. Please verify.") }

  ##### Define min number of samples #####
  #Calculate min samples based on percent if provided
  if(!is.null(min.pct)){
    min.sample <- round(nrow(dat$samples)*min.pct/100, digits=0)
  }

  ##### List not rare genes #####
  # Convert counts to counts per million
  dat.cpm <- cpm(dat$counts)
  # Calculate number of samples that meet the cutoff per gene
  not.rare.samples <- rowSums(dat.cpm >= min.CPM)
  # List not rare genes to be RETAINED
  not.rare.genes <- names(not.rare.samples[not.rare.samples >= min.sample])

  ##### Filter data to remove rare genes #####
  dat.filter <- dat
  # Filter counts
  dat.filter$counts <- dat.filter$counts[rownames(dat.filter$counts) %in% not.rare.genes,]

  #If gene info exists, filter as well
  if(!is.null(dat.filter$genes)){
    #Gene name column in gene info table?
    if(gene.var %notin% colnames(dat.filter$genes)){
      stop("Gene name varible not present in gene info (dat$genes)") }

    # Filter gene key
    dat.filter$genes <- dat.filter$genes[dat.filter$genes[,gene.var] %in% not.rare.genes,]
  }

  ##### Save to environment #####
  assign(name, dat.filter, envir = .GlobalEnv)
}
