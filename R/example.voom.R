#' kimma example EList.
#'
#' @format Formal class 'EList' [package "limma"] with 1 slot:
#' \enumerate{
#' \item \strong{genes} A data frame with 1000 rows and 5 columns
#' \describe{
#'   \item{hgnc_symbol}{character. Current approved HGNC symbol.}
#'   \item{Previous symbols}{character. Previous HGNC symbols.}
#'   \item{Alias symbols}{character. Alias HGNC symbols.}
#'   \item{gene_biotype}{character. Gene product type. All = protein-coding.}
#'   \item{geneName}{character. ENSEMBL gene ID. Matches row names in E.}
#'       }
#'
#' \item \strong{targets} A data frame with 12 rows and 8 columns
#' \describe{
#'   \item{group}{factor. No grouping was provided. All = 1.}
#'   \item{lib.size}{numeric. Total library size for this 1000 gene subset.}
#'   \item{norm.factors}{numeric. TMM normalizatin factors.}
#'   \item{libID}{character. Unique library ID. Matches column names in E.}
#'   \item{donorID}{character. Donor ID.}
#'   \item{median_cv_coverage}{numeric. Median coefficient of variation of coverage. Quality metric for sequencing libraries calculated from original full data set.}
#'   \item{virus}{Factor. Media samples with no virus ("none") vs virus-infected samples ("HRV").}
#'   \item{asthma}{Factor. Asthma vs healthy.}
#'   }
#'
#' \item \strong{E} A matrix with 1000 rows and 12 columns
#' \describe{
#'   \item{rownames}{character. ENSEMBL gene ID.}
#'   \item{lib1}{numeric. log2 CPM in library 1.}
#'   \item{lib2}{numeric. log2 CPM in library 2.}
#'   \item{lib3}{numeric. log2 CPM in library 3.}
#'   \item{lib4}{numeric. log2 CPM in library 4.}
#'   \item{lib5}{numeric. log2 CPM in library 5.}
#'   \item{lib6}{numeric. log2 CPM in library 6.}
#'   \item{lib7}{numeric. log2 CPM in library 7.}
#'   \item{lib8}{numeric. log2 CPM in library 8.}
#'   \item{lib9}{numeric. log2 CPM in library 9.}
#'   \item{lib10}{numeric. log2 CPM in library 10.}
#'   \item{lib11}{numeric. log2 CPM in library 11.}
#'   \item{lib12}{numeric. log2 CPM in library 12.}
#'       }
#'
#' \item \strong{weights} A matrix with 1000 rows and 6 columns
#' \describe{
#'   \item{1}{numeric. limma gene weights for library 1.}
#'   \item{2}{numeric. limma gene weights for library 2.}
#'   \item{3}{numeric. limma gene weights for library 3.}
#'   \item{4}{numeric. limma gene weights for library 4.}
#'   \item{5}{numeric. limma gene weights for library 5.}
#'   \item{6}{numeric. limma gene weights for library 6.}
#'       }
#'
#' \item \strong{design} A matrix with 6 rows and 1 column
#' \describe{
#'   \item{GrandMean}{numeric. limma default design matrix.}
#'       }
#' }
#' @source \url{https://github.com/altman-lab/P259_pDC_public}
#' @references Dill-McFarland et al. 2022 Eosinophil-mediated suppression and Anti-IL-5 enhancement of plasmacytoid dendritic cell interferon responses in asthma. J Allergy Clin Immunol. 150(3):666-675. doi: 10.1016/j.jaci.2022.03.025
#' @description A limma EList data set containing normalized log2 RNA-seq counts. RNA-seq of human dendritic cells cultured with and without virus. Samples from 3 donors and a random subset of 1000 genes were selected. Counts are TMM normalized log2 counts per million (CPM).
#' @docType data
#' @name example.voom
#' @keywords datasets
"example.voom"
