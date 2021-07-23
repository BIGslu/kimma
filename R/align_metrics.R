#' Extract and format cleaning and alignment metrics
#'
#' Extract data from FastQC trim settings, Picard, samtools flagstat, and featureCounts output by the RNA-seq fastq pipeline
#'
#' @param data.dir Character string of directory containing all associated files
#' @param trim Logical if should include FastQC trim settings
#' @param bam Logical if should include samtools flagstat for raw alignments
#' @param picard Logical if should include Picard for raw alignments
#' @param bam.filter Logical if should include samtools flagstat for filtered alignments
#' @param count Logical if should include featureCounts total reads in genes
#'
#' @return Data frame cleaning and alignment metrics for all libraries
#' @export

align_metrics <- function(data.dir=NULL,
                          trim=TRUE, bam=TRUE, picard=TRUE, bam.filter=TRUE, count=TRUE){
  #Blank df for results
  summ.all <- data.frame(libID=NA)
  CODING_BASES <- PCT_PF_ALIGNED <- PF_ALIGNED_BASES <- PF_BASES <- Status <- X1 <- i <- libID <- name <- value <- . <- NULL

  #### Raw and adapter trimmed ####
  if(trim){
    # Raw and adapter trimmed seqs
    trim.files <- list.files(data.dir, pattern="*settings",
                             all.files=FALSE, full.names=TRUE, recursive=TRUE) %>%
      gsub("//", "/", .)

    trim.summ <- data.frame()

    for (file in trim.files){
      libID.temp <- basename(file) %>%
        gsub(".settings", "", .)

      trim.summ.temp <- readr::read_table(file, col_names = FALSE,
                                   col_types = readr::cols()) %>%
        dplyr::filter(startsWith(X1, 'Total number of read pairs') |
                 startsWith(X1, 'Number of retained reads')) %>%
        tidyr::separate(X1, into=c("metric","value"), sep=": ") %>%
        dplyr::mutate(libID = libID.temp) %>%
        tidyr::pivot_wider(names_from = "metric", values_from = "value") %>%
        dplyr::rename(raw="Total number of read pairs",
               trim="Number of retained reads") %>%
        dplyr::mutate_at(dplyr::vars(raw, trim), as.numeric) %>%
        dplyr::mutate(raw=raw*2)

      trim.summ <- dplyr::bind_rows(trim.summ, trim.summ.temp)
    }

    # Combine into results df
    summ.all <- dplyr::full_join(summ.all, trim.summ, by = "libID")
  }

  #### Alignment ####
  if(bam){
    align.file <- list.files(data.dir, pattern="summary.alignment|bam.summary",
                             all.files=FALSE, full.names=TRUE, recursive=TRUE) %>%
      gsub("//", "/", .)

    align.lib <- readr::read_tsv(align.file, col_names = FALSE,
                          col_types = readr::cols()) %>%
      dplyr::filter(grepl("bam",X1)) %>%
      dplyr::distinct(X1)

    # Align summary
    align.summ <- readr::read_tsv(align.file, col_names = FALSE,
                           col_types = readr::cols()) %>%
      #Get sample name from filename
      dplyr::mutate(libID = factor(X1, levels=align.lib$X1),
             libID = basename(as.character(libID)),
             libID = gsub("_Aligned.sortedByCoord.out.bam", "", libID)) %>%
      tidyr::fill(libID) %>%
      #Separate data to column
      tidyr::separate(X1, into=c("h","i"), sep=" \\+ 0 ", fill="right") %>%
      tidyr::drop_na(i) %>%
      #Recode data types (f)
      tidyr::separate(i, into=c("i"), sep="[(]", extra="drop") %>%
      dplyr::mutate(i = forcats::fct_recode(factor(i),
                            to.be.aligned="in total ",
                            secondary.align="secondary",
                            chimeric.align="supplementary",
                            PCR.dups="duplicates",
                            align="mapped ",
                            paired="paired in sequencing",
                            R1.paired="read1", R2.paired="read2",
                            align.paired="properly paired ",
                            both.align.paired= "with itself and mate mapped",
                            one.align.paired="singletons " ,
                            both.align.paired.diffCHR="with mate mapped to a different chr",
                            both.align.paired.diffCHR.mapq="with mate mapped to a different chr ")) %>%
      tidyr::pivot_wider(names_from = "i", values_from = "h") %>%
      dplyr::mutate_at(dplyr::vars(-libID), as.numeric)

    # Combine into results df
    summ.all <- dplyr::full_join(summ.all, align.summ, by = "libID")
  }

  #### BAM summary ####
  if(picard){
    bam1 <- list.files(data.dir, pattern="bam.metrics",
                       all.files=FALSE, full.names=TRUE, recursive=TRUE) %>%
      gsub("//", "/", .)

    bam.colnames <- readr::read_delim(bam1, delim = " ", col_names=FALSE)[7,1] %>%
      tidyr::separate(X1, into=as.character(c(1:30)), sep="\t") %>%
      unlist(use.names = FALSE)

    bam.summ <- readr::read_delim(bam1, delim = " ", col_names=FALSE, comment = "#",
                           col_types = readr::cols()) %>%
      #Get sample name from file name
      dplyr::mutate(libID = factor(X1, levels=align.lib$X1),
             libID = basename(as.character(libID)),
             libID = gsub("_Aligned.sortedByCoord.out.bam", "", libID)) %>%
      tidyr::fill(libID) %>%
      #Sep data columns
      tidyr::separate(X1, into = bam.colnames, sep="\t", fill="right") %>%
      #Remove histogram columns
      tidyr::drop_na(CODING_BASES) %>%
      #Keep only data rows
      dplyr::filter(!startsWith(as.character(PF_BASES), 'PF')) %>%
      #Make numeric
      dplyr::mutate_at(dplyr::vars(-libID), as.numeric) %>%
      #Calculate perc align
      dplyr::mutate(PCT_PF_ALIGNED = PF_ALIGNED_BASES/PF_BASES) %>%
      dplyr::select(libID, PCT_PF_ALIGNED, dplyr::everything())

    # Combine into results df
    summ.all <- dplyr::full_join(summ.all, bam.summ, by = "libID")
  }

  #### Filtered paired alignment ####
  if(bam.filter){
    align.file2 <- list.files(data.dir, pattern="summary.align.filter|bam.filter.summary",
                              all.files=FALSE, full.names=TRUE, recursive=TRUE) %>%
      gsub("//", "/", .)

    align.lib2 <- readr::read_tsv(align.file2, col_names = FALSE,
                           col_types = readr::cols()) %>%
      dplyr::filter(grepl("bam",X1)) %>%
      dplyr::distinct(X1)

    # Align summary
    align.summ2 <- readr::read_tsv(align.file2, col_names = FALSE,
                            col_types = readr::cols()) %>%
      #Get sample name from filename
      dplyr::mutate(libID = factor(X1, levels=align.lib2$X1),
             libID = basename(as.character(libID)),
             libID = gsub("_filter_paired.bam|_filter.bam", "", libID)) %>%
      tidyr::fill(libID) %>%
      #Separate data to column
      tidyr::separate(X1, into=c("h","i"), sep=" \\+ 0 ", fill="right") %>%
      tidyr::drop_na(i) %>%
      #Recode data types (f)
      tidyr::separate(i, into=c("i"), sep="[(]", extra="drop") %>%
      dplyr::mutate(i = forcats::fct_recode(factor(i),
                            filtered.to.be.aligned="in total ",
                            filtered.secondary.align="secondary",
                            filtered.chimeric.align="supplementary",
                            filtered.PCR.dups="duplicates",
                            filtered.align="mapped ",
                            filtered.paired="paired in sequencing",
                            filtered.R1.paired="read1", filtered.R2.paired="read2",
                            filtered.align.paired="properly paired ",
                            filtered.both.align.paired= "with itself and mate mapped",
                            filtered.one.align.paired="singletons " ,
                            filtered.both.align.paired.diffCHR="with mate mapped to a different chr",
                            filtered.both.align.paired.diffCHR.mapq="with mate mapped to a different chr ")) %>%
      tidyr::pivot_wider(names_from = "i", values_from = "h") %>%
      dplyr::mutate_at(dplyr::vars(-libID), as.numeric)

    # Combine into results df
    summ.all <- dplyr::full_join(summ.all, align.summ2, by = "libID")
  }

  #### Gene counts ####
  if(count){
    count.file <- list.files(data.dir, pattern="featurecounts",
                             all.files=FALSE, full.names=TRUE, recursive=TRUE) %>%
      gsub("//", "/", .)

    count.summ <- readr::read_tsv(count.file, col_types = readr::cols()) %>%
      tidyr::pivot_longer(-Status) %>%
      #Get sample names from file
      dplyr::mutate(libID = basename(name),
             libID = gsub("_filter.paired.bam|_filter.bam", "", libID)) %>%
      dplyr::select(-name) %>%
      #Keep rows with non-zero values
      dplyr::filter(value !=0) %>%
      #Transpose
      tidyr::pivot_wider(names_from = Status, values_from = value)

    # Combine into results df
    summ.all <- dplyr::full_join(summ.all, count.summ, by = "libID")
  }

  #### Save ####
  #remove columns that are all blank or 0
  summ.all.filter <- summ.all %>%
    tidyr::drop_na(libID) %>%
    dplyr::mutate_at(dplyr::vars(-libID), as.numeric) %>%
    dplyr::select_if(~sum(!is.na(.)) > 0) %>%
    dplyr::select_if(~sum(.!= 0) > 0)

  return(summ.all.filter)
}
