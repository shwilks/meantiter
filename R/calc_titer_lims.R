
titer_fit_options <- function(
  min_titer_possible = -Inf,
  max_titer_possible = Inf
) {
  
  list(
    min_titer_possible = min_titer_possible,
    max_titer_possible = max_titer_possible
  )
  
}

#' Get titer limits
#'
#' Function for getting upper and lower limits of measured titers on the log scale.
#'
#' @param titers A numeric/character vector or matrix of titer measurements.
#' @param fit_opts A list of fitting options (see details).
#'
#' @return Returns a list of length two with values max_titers and min_titers, giving the
#' numeric vectors of the upper and lower bounds of the titers on the log scale.
#'
#' @details This function assumes that HI measurements were performed in 2-fold dilution
#' steps and converts them to the log scale using the formula:
#' \deqn{a + b}
#' Hence an HI titer of 20, which would convert to 1 via the transformation above, would be
#' assumed to have upper and lower limits of 1.5 and 0.5 respectively.
#'
#' In the case of non-detectable titers, such as <10, the lower bound of the measured value
#' is taken from the parameter \code{min_titer_possible}, defaulting to the value found from
#' a call to \code{get_lndscp_fit_defaults()}. For a greater than value, i.e. >1280, the
#' upper bound of the value is taken from the parameter \code{max_titer_possible}. You can
#' set different defaults by passing them as named arguments to the list, as shown in the
#' examples.
#'
#' @examples
#' # Calculate the titer limits of a set of HI titers
#' titer_lims <- get_titer_lims(titers = c("20", "320", "<10", ">1280"))
#'
#' # Calculate the titer limits assuming non-default upper and lower bounds for non-detectable
#' # and greater-than titers.
#' titer_lims <- get_titer_lims(titers = c("20", "320", "<10", ">1280"),
#'                              fit_opts = list(min_titer_possible = -Inf,
#'                                              max_titer_possible = 14))
#' @export
calc_titer_lims <- function(
  titers,
  dilution_stepsize,
  options = list()
) {
  
  # Get options
  options <- do.call(titer_fit_options, options)
  
  # Find less than and greater than titers and convert them to a numeric form
  lessthan_titers <- grepl(x = titers, pattern = "<")
  morethan_titers <- grepl(x = titers, pattern = ">")
  na_titers       <- is.na(titers) | grepl(x = titers, pattern = "\\*")
  
  numeric_titers <- titers
  numeric_titers[na_titers] <- NA
  numeric_titers <- as.numeric(gsub("(<|>)","",numeric_titers))
  
  # Convert titers to the log scale
  log_titers <- log2(numeric_titers / 10)
  log_titers[lessthan_titers] <- log_titers[lessthan_titers] - dilution_stepsize
  log_titers[morethan_titers] <- log_titers[morethan_titers] + dilution_stepsize
  max_titers <- log_titers + dilution_stepsize / 2
  min_titers <- log_titers - dilution_stepsize / 2
  min_titers[lessthan_titers] <- options$min_titer_possible
  max_titers[morethan_titers] <- options$max_titer_possible
  
  list(
    log_titers = log_titers,
    max_titers = max_titers,
    min_titers = min_titers
  )
  
}


#' @export
calc_titer_diff_lims <- function(
  titers1,
  titers2,
  dilution_stepsize,
  fit_opts = list()
) {
  
  # Get limits to titer measurements
  titer1_lims  <- calc_titer_lims(titers1, dilution_stepsize, fit_opts)
  titer2_lims <- calc_titer_lims(titers2, dilution_stepsize, fit_opts)
  
  # Compute maximum and mininmum differences
  max_diffs <- titer2_lims$max_titers - titer1_lims$min_titers
  min_diffs <- titer2_lims$min_titers - titer1_lims$max_titers
  
  logtiters1 <- titer1_lims$log_titers
  logtiters2 <- titer2_lims$log_titers
  
  # Return result
  list(
    max_diffs = max_diffs,
    min_diffs = min_diffs,
    logtiter_diffs = logtiters2 - logtiters1
  )
  
}



