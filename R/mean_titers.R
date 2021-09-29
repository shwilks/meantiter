
#' Calculate the geometric mean titer of a set of titers
#'
#' @param titers A vector of titers
#' @param method The method to use when dealing with censored titers (like
#'   `<10`), one of "replace_nd", "exclude_nd", "truncated_normal"
#' @param level The confidence level to use when calculating confidence intervals
#' @param sd By setting this you can fix the standard deviation assumed when finding 
#'   parameters for the normal distribution when using the "truncated_normal" approach
#' @param dilution_stepsize The dilution stepsize used in the assay, see `calc_titer_lims()`
#' @param options A named list of options to pass to `titer_fit_options()`
#'
#' @export
mean_titers <- function(
  titers, 
  method,
  level = 0.95,
  sd = NA,
  dilution_stepsize,
  options = list()
) {
  
  # Remove NA titers
  na_titers <- is.na(titers) | titers == "*"
  titers <- titers[!na_titers]
  
  # If length titers = 1 just return the titer or NA if thresholded
  if (length(titers) == 0) {
    return(
      list(
        mean = NA,
        sd = NA,
        mean_lower = NA,
        mean_upper = NA
      )
    )
  } else if (length(unique(titers)) == 1) {
    if (grepl("<|>|\\*", titers[1])) {
      return(
        list(
          mean = NA,
          sd = NA,
          mean_lower = NA,
          mean_upper = NA
        )
      )
    } else {
      return(
        list(
          mean = log2(as.numeric(titers[1]) / 10),
          sd = NA,
          mean_lower = NA,
          mean_upper = NA
        )
      )
    }
  }
  
  switch(
    method,
    "replace_nd" = mean_titers_replace_nd(
      titers = titers, 
      level = level, 
      dilution_stepsize = dilution_stepsize,
      sd = sd
    ),
    "exclude_nd" = mean_titers_exclude_nd(
      titers = titers, 
      level = level, 
      dilution_stepsize = dilution_stepsize,
      sd = sd
    ),
    "truncated_normal" = mean_titers_truncated_normal(
      titers = titers, 
      level = level, 
      dilution_stepsize = dilution_stepsize,
      sd = sd,
      options = options
    )
  )
  
}


mean_titers_replace_nd <- function(
  titers,
  level = 0.95,
  dilution_stepsize,
  sd = NA
) {
  
  lessthans <- substr(titers, 1, 1) == "<"
  morethans <- substr(titers, 1, 1) == ">"
  
  titers[lessthans | morethans] <- substr(
    x = titers[lessthans | morethans], 
    start = 2, 
    stop = nchar(titers[lessthans | morethans])
  )
  
  logtiters <- log2(as.numeric(titers) / 10)
  logtiters[lessthans] <- logtiters[lessthans] - dilution_stepsize
  logtiters[morethans] <- logtiters[morethans] + dilution_stepsize
  
  result <- Hmisc::smean.cl.normal(logtiters, conf.int = level)
  list(
    mean = unname(result["Mean"]),
    sd = sd(logtiters),
    mean_lower = unname(result["Lower"]),
    mean_upper = unname(result["Upper"])
  )
  
}



mean_titers_exclude_nd <- function(
  titers,
  level = 0.95,
  dilution_stepsize,
  sd = NA
) {
  
  lessthans <- substr(titers, 1, 1) == "<"
  morethans <- substr(titers, 1, 1) == ">"
  titers <- titers[!lessthans & !morethans]
  
  logtiters <- log2(as.numeric(titers) / 10)
  
  result <- Hmisc::smean.cl.normal(logtiters, conf.int = level)
  list(
    mean = result[["Mean"]],
    sd = sd(logtiters),
    mean_lower = result[["Lower"]],
    mean_upper = result[["Upper"]]
  )
  
}



mean_titers_truncated_normal <- function(
  titers,
  level = 0.95,
  dilution_stepsize,
  sd = NA,
  options = list()
) {
  
  # Get the titer limits
  titerlims <- calc_titer_lims(
    titers = titers,
    dilution_stepsize = dilution_stepsize,
    options = options
  )
  
  if (is.na(sd)) {
    
    start_pars <- list(
      mean = mean(titerlims$log_titers, na.rm = TRUE),
      sd = sd(titerlims$log_titers, na.rm = TRUE)
    )
    
    fixed_pars <- NULL
    
  } else {
    
    start_pars <- list(
      mean = mean(titerlims$log_titers, na.rm = TRUE)
    )
    
    fixed_pars <- list(
      sd = sd
    )
    
  }
  
  # Setup output
  output <- list(
    mean = NA,
    sd = NA,
    mean_lower = NA,
    mean_upper = NA
  )
  
  try({
    
    out <- capture.output({
      result <- fitdistrplus::fitdistcens(
        censdata = data.frame(
          left = titerlims$min_titers,
          right = titerlims$max_titers
        ),
        start = start_pars,
        fix.arg = fixed_pars,
        distr = "norm"
      )
      
      result_ci <- confint(result, level = level)
      
      output$mean <- result$estimate[["mean"]]
      output$mean_lower <- result_ci["mean", 1]
      output$mean_upper <- result_ci["mean", 2]
      
      if (is.na(sd)) {
        output$sd <- result$estimate[["sd"]]
      } else {
        output$sd <- sd
      }
      
    })
    
  }, silent = TRUE)
  
  # Return output
  output
  
}

