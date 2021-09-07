
#' @export
mean_titers <- function(
  titers, 
  method,
  level = 0.95,
  sd = NA,
  dilution_stepsize
) {
  
  # Remove NA titers
  na_titers <- is.na(titers) | titers == "*"
  titers <- titers[!na_titers]
  
  switch(
    method,
    "maxlikelihood" = mean_titers_maxlikelihood(
      titers = titers, 
      level = level, 
      dilution_stepsize = dilution_stepsize,
      sd = sd
    ),
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
      sd = sd
    )
  )
  
}



mean_titers_maxlikelihood <- function(
  titers,
  level = 0.95,
  dilution_stepsize,
  sd = NA
) {
  
  # Get the titer limits
  titerlims <- calc_titer_lims(titers, dilution_stepsize)
  
  # Assign starting parameters
  start_mean <- mean(titerlims$log_titers)
  start_sd   <- sd(titerlims$log_titers)
  
  # Calculate the titer likelihood
  result <- nlminb(
    start = c(start_mean, start_sd),
    objective = calc_mean_titer_negll_by_par,
    max_titers = titerlims$max_titers,
    min_titers = titerlims$min_titers,
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
      max_titers = titerlims$max_titers,
      min_titers = titerlims$min_titers,
      titer_sd = NA,
      upper = c(result$par[1], Inf),
      lower = c(-Inf, 0.01),
      target_negll = result$objective + qchisq(level, 1)/2
    )

    upper_ci_result <- nlminb(
      start = c(result$par[1] + 0.01, result$par[2]),
      objective = calc_mean_titer_ci_by_par,
      max_titers = titerlims$max_titers,
      min_titers = titerlims$min_titers,
      titer_sd = NA,
      upper = c(Inf, Inf),
      lower = c(result$par[1], 0.01),
      target_negll = result$objective + qchisq(level, 1)/2
    )
    
    list(
      mean = result$par[1],
      sd = result$par[2],
      mean_lower = lower_ci_result$par[1],
      mean_upper = upper_ci_result$par[1]
    )
    
  }
  
  
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
  sd = NA
) {
  
  # Get the titer limits
  titerlims <- calc_titer_lims(titers, dilution_stepsize)
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

