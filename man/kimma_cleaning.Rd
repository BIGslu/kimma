% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kimma_cleaning.R
\name{kimma_cleaning}
\alias{kimma_cleaning}
\title{kmFit data cleaning}
\usage{
kimma_cleaning(
  dat = NULL,
  kin = NULL,
  patientID = "ptID",
  libraryID = "libID",
  counts = NULL,
  meta = NULL,
  genes = NULL,
  subset.var = NULL,
  subset.lvl = NULL,
  subset.genes = NULL
)
}
\arguments{
\item{dat}{EList object output by voom( ). Contains counts (dat$E), meta (dat$targets), and genes (dat$genes).}

\item{kin}{Matrix with pairwise kinship values between individuals. Must be numeric with rownames.}

\item{patientID}{Character of variable name to match dat$targets to kinship row and column names.}

\item{libraryID}{Character of variable name to match dat$targets to dat$E colnames

Alternate data if not using EList object}

\item{counts}{Matrix of normalized expression. Rows are genes, columns are libraries.}

\item{meta}{Matrix or data frame of sample and individual metadata.}

\item{genes}{Matrix or data frame of gene metadata.

Subset data (optional)}

\item{subset.var}{Character list of variable name(s) to filter data by.}

\item{subset.lvl}{Character list of variable value(s) or level(s) to filter data to. Must match order of subset.var}

\item{subset.genes}{Character vector of genes to include in models.}
}
\value{
Data frame formatted for use in kmFit
}
\description{
Data
}
\keyword{internal}