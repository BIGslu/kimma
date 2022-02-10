#' Run pairwise model comparisons within lmekin model
#'
#' @param contrast.var Character vector of variable in model to run contrasts of
#' @param to.model.gene Formatted data from kimma_cleaning( ), subset to gene of interest
#' @param patientID Character of variable name to match dat$targets to kinship row and column names.
#' @param to.model.ls Formatted data from kimma_cleaning( )
#' @param gene Character name of gene of interest
#' @param use.weights Logical if gene specific weights should be used in model
#'
#' @return data frame with contrast model results
#' @keywords internal

kmFit_contrast_kin <- function(contrast.var, to.model.gene, patientID,
                               to.model.ls, gene, use.weights){
  contrast.kin <- contrast.i <- contrast.i.format <- V1 <- V2 <- variable <- contrast <- NULL
  to.model.contrast <- to.model.gene

  #emmeans does not work for lmekin object. Instead, run pairwise contrasts as subsets
  for(contrast.i in contrast.var){
    #Class
    i.split <- strsplit(contrast.i, split=":")[[1]]
    contrast.isnt.numeric <- !(unlist(lapply(to.model.contrast[,i.split], is.numeric)))

    #only run on categorical variable
    if(any(contrast.isnt.numeric)){
      #non-interaction variables
      if(!grepl(":", contrast.i)){
        model.contrast <- paste("expression~",contrast.i, "+(1|", patientID, ")", sep="")
        to.model.contrast.format <- to.model.contrast
        contrast.i.format <- contrast.i
      } else{
        model.contrast <- paste("expression~", i.split[1], "_", i.split[2],
                                "+(1|", patientID, ")", sep="")
        #Create variable in data
        to.model.contrast.format <- to.model.contrast
        to.model.contrast.format[,paste0(i.split[1], "_", i.split[2])] <-
          paste(unlist(to.model.contrast[,i.split[1]]),
                unlist(to.model.contrast[,i.split[2]]), sep="_")
        contrast.i.format <- gsub(":","_",contrast.i)
      }

        #list all possible combos
        levels.OI <- as.character(unlist(unique(to.model.contrast.format[,contrast.i.format])))
        contrast.combo <- as.data.frame(t(utils::combn(levels.OI, m=2))) %>%
          dplyr::mutate(combo = paste(contrast.i.format, V2, " - ",contrast.i.format, V1, sep=""))

        if(nrow(contrast.combo)>0){
          for(row.i in 1:nrow(contrast.combo)){
            to.model.current <- to.model.contrast.format
            contrast.kin.temp <- NULL
            #filter data to each pairwise comparison of interest
            to.model.contrast.gene <- dplyr::filter(to.model.current,
                                      get(contrast.i.format) %in% unlist(contrast.combo[row.i,]))

            to.keep <- unique(unlist(to.model.current[,patientID]))
            kin.contrast <- to.model.ls[["kin.subset"]][rownames(to.model.ls[["kin.subset"]]) %in%
                                                          to.keep, to.keep]

            contrast.kin.temp <- tryCatch({
              kimma_lmekin(model.contrast, to.model.gene.contrast, gene, kin.contrast, use.weights)
            }, error=function(e){})

            if(!is.null(contrast.kin.temp)){
            contrast.kin <- contrast.kin.temp[["results"]] %>%
              dplyr::filter(variable != "(Intercept)") %>%
              dplyr::mutate(variable = contrast.i.format,
                            contrast = contrast.combo[row.i, "combo"],
                            model = "lmekin.contrast") %>%
              dplyr::mutate(contrast = gsub(contrast.i.format,"",contrast)) %>%
              dplyr::bind_rows(contrast.kin)
            }}
          }
    }
  }
  rownames(contrast.kin) <- NULL
  return(contrast.kin)
}
