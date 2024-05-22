#' Process MS data using Bayesian ML: clean, normalize and summarize before differential analysis
#'
#' @param raw name of the raw (input) data set.
#' @param logTrans base of logarithm transformation: 2 (default) or 10.
#' @param normalization normalization to remove systematic bias between MS runs.
#' There are three different normalizations supported:
#' 'equalizeMedians' (default) represents constant normalization (equalizing the medians)
#' based on reference signals is performed.
#' 'quantile' represents quantile normalization based on reference signals
#' 'globalStandards' represents normalization with global standards proteins.
#' If FALSE, no normalization is performed.
#' @param nameStandards optional vector of global standard peptide names.
#' Required only for normalization with global standard peptides.
#' @param fix_missing Optional, same as the `fix_missing` parameter in MSstatsConvert::MSstatsBalancedDesign function
#'
#' @import MSstats
#'
#' @export
#'
#' @examples
#' #TODO
#' NULL
dataProcessBayes <- function(data,
                             logTrans=2,
                             normalization="equalizeMedians",
                             nameStandards=NULL,
                             fix_missing=NULL
                             ) {

  # Normalization
  peptides_dict = makePeptidesDictionary(as.data.table(unclass(data)), normalization)
  input = MSstatsPrepareForDataProcess(data, logTrans, fix_missing)
  input = MSstatsNormalize(input, normalization, peptides_dict, nameStandards)
  input = MSstatsMergeFractions(input)

  dpc_betas = calculateMNARCurve(input)

  summarized_results = MSstatsBayesSummarize(input, dpc_betas)

  MSstatsBayesSummarizationOutput

  return(summarized_results)

}
