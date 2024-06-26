% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dataProcessBayes.R
\name{dataProcessBayes}
\alias{dataProcessBayes}
\title{Process MS data using Bayesian ML: clean, normalize and summarize before differential analysis}
\usage{
dataProcessBayes(
  data,
  logTrans = 2,
  normalization = "equalizeMedians",
  nameStandards = NULL,
  fix_missing = NULL
)
}
\arguments{
\item{logTrans}{base of logarithm transformation: 2 (default) or 10.}

\item{normalization}{normalization to remove systematic bias between MS runs.
There are three different normalizations supported:
'equalizeMedians' (default) represents constant normalization (equalizing the medians)
based on reference signals is performed.
'quantile' represents quantile normalization based on reference signals
'globalStandards' represents normalization with global standards proteins.
If FALSE, no normalization is performed.}

\item{nameStandards}{optional vector of global standard peptide names.
Required only for normalization with global standard peptides.}

\item{fix_missing}{Optional, same as the `fix_missing` parameter in MSstatsConvert::MSstatsBalancedDesign function}

\item{raw}{name of the raw (input) data set.}
}
\description{
Process MS data using Bayesian ML: clean, normalize and summarize before differential analysis
}
\examples{
#TODO
NULL
}
