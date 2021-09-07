
#' @export
mean_titer_diffs <- function(
  titers1, 
  titers2,
  method,
  level = 0.95,
  dilution_stepsize
) {
  
  # Remove NA titers
  na_titers <- is.na(titers1) | titers1 == "*" | is.na(titers2) | titers2 == "*"
  titers1 <- titers1[!na_titers]
  titers2 <- titers2[!na_titers]
  
  switch(
    method,
    "maxlikelihood" = mean_titer_diffs_maxlikelihood(titers1, titers2, level, dilution_stepsize),
    "replace_nd" = mean_titer_diffs_replace_nd(titers1, titers2, level, dilution_stepsize),
    "exclude_nd" = mean_titer_diffs_exclude_nd(titers1, titers2, level, dilution_stepsize),
    "truncated_normal" = mean_titer_diffs_truncated_normal(titers1, titers2, level, dilution_stepsize)
  )
  
}



mean_titer_diffs_maxlikelihood <- function(
  titers1,
  titers2,
  level = 0.95,
  dilution_stepsize
) {
  
  # Get the titer limits
  titerlims <- calc_titer_diff_lims(titers1, titers2, dilution_stepsize)
  
  # Assign starting parameters
  start_mean <- mean(titerlims$logtiter_diffs)
  start_sd   <- sd(titerlims$logtiter_diffs)
  
  # Calculate the titer likelihood
  result <- nlminb(
    start = c(start_mean, start_sd),
    objective = calc_mean_titer_negll_by_par,
    max_titers = titerlims$max_diffs,
    min_titers = titerlims$min_diffs,
    titer_sd = NA,
    upper = c(Inf, Inf),
    lower = c(-Inf, 0.01)
  )
  
  if (is.na(level)) {
    
    list(
      mean = result$par[1],
      sd = result$par[2],
      mean_lower = NA,
      mean_upper = NA
    )
    
  } else {
    
    lower_ci_result <- nlminb(
      start = c(result$par[1] - 0.01, result$par[2]),
      objective = calc_mean_titer_ci_by_par,
      max_titers = titerlims$max_diffs,
      min_titers = titerlims$min_diffs,
      titer_sd = NA,
      upper = c(result$par[1], Inf),
      lower = c(-Inf, 0.01),
      target_negll = result$objective + qchisq(level, 1)/2
    )
    
    upper_ci_result <- nlminb(
      start = c(result$par[1] + 0.01, result$par[2]),
      objective = calc_mean_titer_ci_by_par,
      max_titers = titerlims$max_diffs,
      min_titers = titerlims$min_diffs,
      titer_sd = NA,
      upper = c(Inf, Inf),
      lower = c(result$par[1], 0.01),
      target_negll = result$objective + qchisq(level, 1)/2
    )
    
    list(
      mean_diff = result$par[1],
      sd = result$par[2],
      mean_diff_lower = lower_ci_result$par[1],
      mean_diff_upper = upper_ci_result$par[1]
    )
    
  }
  
  
}


mean_titer_diffs_replace_nd <- function(
  titers1,
  titers2,
  level = 0.95,
  dilution_stepsize
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
    mean_diff = unname(result["Mean"]),
    sd = sd(logtiters),
    mean_diff_lower = unname(result["Lower"]),
    mean_diff_upper = unname(result["Upper"])
  )
  
}



mean_titer_diffs_exclude_nd <- function(
  titers1,
  titers2,
  level = 0.95,
  dilution_stepsize
) {
  
  lessthans <- substr(titers, 1, 1) == "<"
  morethans <- substr(titers, 1, 1) == ">"
  titers <- titers[!lessthans & !morethans]
  
  logtiters <- log2(as.numeric(titers) / 10)
  
  result <- Hmisc::smean.cl.normal(logtiters, conf.int = level)
  list(
    mean_diff = result[["Mean"]],
    sd = sd(logtiters),
    mean_diff_lower = result[["Lower"]],
    mean_diff_upper = result[["Upper"]]
  )
  
}



mean_titer_diffs_truncated_normal <- function(
  titers1,
  titers2,
  level = 0.95,
  dilution_stepsize
) {
  
  # Get the titer limits
  titerlims <- calc_titer_diff_lims(
    titers1, 
    titers2, 
    dilution_stepsize
  )
  
  # Setup output
  output <- list(
    mean_diff = NA,
    sd = NA,
    mean_diff_lower = NA,
    mean_diff_upper = NA
  )
  
  try({
    
    out <- capture.output({
      result <- fitdistrplus::fitdistcens(
        censdata = data.frame(
          left = titerlims$min_diffs,
          right = titerlims$max_diffs
        ),
        start = list(
          mean = mean(titerlims$logtiter_diffs),
          sd = sd(titerlims$logtiter_diffs)
        ),
        distr = "norm",
        
      )
      
      result_ci <- confint(result, level = level)
      
      output <- list(
        mean_diff = result$estimate[["mean"]],
        sd = result$estimate[["sd"]],
        mean_diff_lower = result_ci["mean", 1],
        mean_diff_upper = result_ci["mean", 2]
      )
    })
    
  }, silent = TRUE)
  
  # Return output
  output
  
}

