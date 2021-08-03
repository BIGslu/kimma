#' kimma example kinship.
#'
#' @format A matrix with 3 rows and 3 variables:
#' \describe{
#'   \item{rowname}{Donor ID. Same as column names}
#'   \item{donor1}{numeric kinship (0-1) with donor 1}
#'   \item{donor2}{numeric kinship (0-1) with donor 2}
#'   \item{donor3}{numeric kinship (0-1) with donor 3}
#' }
#' @description Matrix of pairwise kinship values between donor 1,2,3. Values are dummy data with 1 for self comparison, 0.5 for siblings, and 0.1 for unrelated.
#' @docType data
#' @name example.kin
#' @keywords datasets
"example.kin"
