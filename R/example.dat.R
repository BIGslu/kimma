#' kimma example DGEList.
#'
#' @format Formal class 'DGEList' [package "edgeR"] with 1 slot:
#' \enumerate{
#' \item \strong{counts} A matrix with 1000 rows and 12 columns
#' \describe{
#'   \item{rownames}{character. ENSEMBL gene ID.}
#'   \item{lib1}{integer. Counts in library 1.}
#'   \item{lib2}{integer. Counts in library 2.}
#'   \item{lib3}{integer. Counts in library 3.}
#'   \item{lib4}{integer. Counts in library 4.}
#'   \item{lib5}{integer. Counts in library 5.}
#'   \item{lib6}{integer. Counts in library 6.}
#'   \item{lib7}{integer. Counts in library 7.}
#'   \item{lib8}{integer. Counts in library 8.}
#'   \item{lib9}{integer. Counts in library 9.}
#'   \item{lib10}{integer. Counts in library 10.}
#'   \item{lib11}{integer. Counts in library 11.}
#'   \item{lib12}{integer. Counts in library 12.}
#'       }
#'
#' \item \strong{samples} A data frame with 12 rows and 9 columns
#' \describe{
#'   \item{group}{factor. No grouping was provided. All = 1.}
#'   \item{lib.size}{numeric. Total library size for this 1000 gene subset.}
#'   \item{norm.factors}{numeric. Normalizatin factors. No normalization was completed. All = 1.}
#'   \item{libID}{character. Unique library ID. Matches column names in counts.}
#'   \item{donorID}{character. Donor ID.}
#'   \item{median_cv_coverage}{numeric. Median coefficient of variation of coverage. Quality metric for sequencing libraries calculated from original full data set.}
#'   \item{virus}{Factor. Media samples with no virus ("none") vs virus-infected samples ("HRV").}
#'   \item{asthma}{Character. Asthma vs healthy.}
#'       }
#'
#' \item \strong{genes} A data frame with 1000 rows and 5 columns
#' \describe{
#'   \item{hgnc_symbol}{character. Current approved HGNC symbol.}
#'   \item{Previous symbols}{character. Previous HGNC symbols.}
#'   \item{Alias symbols}{character. Alias HGNC symbols.}
#'   \item{gene_biotype}{character. Gene product type. All = protein-coding.}
#'   \item{geneName}{character. ENSEMBL gene ID. Matches row names in counts.}
#'       }
#' }
#' @source \url{https://github.com/altman-lab/P259_pDC_public}
#' @references Dill-McFarland et al. 2022 Eosinophil-mediated suppression and Anti-IL-5 enhancement of plasmacytoid dendritic cell interferon responses in asthma. J Allergy Clin Immunol. 150(3):666-675. doi: 10.1016/j.jaci.2022.03.025
#' @description An edgeR DGEList data set containing unnormalized RNA-seq counts. RNA-seq of human dendritic cells cultured with and without virus. Samples from 3 donors and a random subset of 1000 genes were selected. Counts are unnormalized.
#' @docType data
#' @name example.dat
#' @keywords datasets
"example.dat"
