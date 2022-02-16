
#' Calculate the mean difference between two paired sets of titers
#'
#' This function is useful for example for calculating differences where you have a set of
#' pre-vaccination and post-vaccination samples and you would like to know the mean response
#' size, accounting for non-detectable values.
#'
#' @param titers1 The first titerset from which to calculate the mean difference
#' @param titers2 The second titerset to which to calculate the mean difference
#' @param method The method to use when dealing with censored titers (like
#'   `<10`), one of "replace_nd", "exclude_nd", "truncated_normal"
#' @param level The confidence level to use when calculating confidence intervals
#' @param dilution_stepsize The dilution stepsize used in the assay, see `calc_titer_lims()`
#'
#' @export
mean_titer_diffs <- function(
  titers1,
  titers2,
  method,
  level = 0.95,
  dilution_stepsize
) {

  # Remove NA titers
  na_titers <- is.na(titers1) | titers1 == "*" | titers1 == "." | is.na(titers2) | titers2 == "*" | titers2 == "."
  titers1 <- titers1[!na_titers]
  titers2 <- titers2[!na_titers]

  # If length titers = 1 just return the titer or NA if thresholded
  if (length(titers1) == 0) {
    return(
      list(
        mean_diff = NA,
        sd = NA,
        mean_diff_lower = NA,
        mean_diff_upper = NA,
        mean_var = NA,
        sd_var = NA
      )
    )
  } else if (length(titers1) == 1 && length(titers2) == 1) {
    if (grepl("<|>", titers1) || titers1 == "*" || titers1 == "." || grepl("<|>", titers2) || titers2 == "*" || titers2 == ".") {
      return(
        list(
          mean_diff = NA,
          sd = NA,
          mean_diff_lower = NA,
          mean_diff_upper = NA,
          mean_var = NA,
          sd_var = NA
        )
      )
    } else {
      return(
        list(
          mean_diff = log2(as.numeric(titers2) / 10) - log2(as.numeric(titers1) / 10),
          sd = NA,
          mean_diff_lower = NA,
          mean_diff_upper = NA,
          mean_var = NA,
          sd_var = NA
        )
      )
    }
  }

  switch(
    method,
    "replace_nd" = mean_titer_diffs_replace_nd(titers1, titers2, level, dilution_stepsize),
    "exclude_nd" = mean_titer_diffs_exclude_nd(titers1, titers2, level, dilution_stepsize),
    "truncated_normal" = mean_titer_diffs_truncated_normal(titers1, titers2, level, dilution_stepsize)
  )

}


mean_titer_diffs_replace_nd <- function(
  titers1,
  titers2,
  level = 0.95,
  dilution_stepsize
) {

  lessthans1 <- substr(titers1, 1, 1) == "<"
  morethans1 <- substr(titers1, 1, 1) == ">"

  titers1[lessthans1 | morethans1] <- substr(
    x = titers1[lessthans1 | morethans1],
    start = 2,
    stop = nchar(titers1[lessthans1 | morethans1])
  )

  logtiters1 <- log2(as.numeric(titers1) / 10)
  logtiters1[lessthans1] <- logtiters1[lessthans1] - dilution_stepsize
  logtiters1[morethans1] <- logtiters1[morethans1] + dilution_stepsize

  lessthans2 <- substr(titers2, 1, 1) == "<"
  morethans2 <- substr(titers2, 1, 1) == ">"

  titers2[lessthans2 | morethans2] <- substr(
    x = titers2[lessthans2 | morethans2],
    start = 2,
    stop = nchar(titers2[lessthans2 | morethans2])
  )

  logtiters2 <- log2(as.numeric(titers2) / 10)
  logtiters2[lessthans2] <- logtiters2[lessthans2] - dilution_stepsize
  logtiters2[morethans2] <- logtiters2[morethans2] + dilution_stepsize

  logtiters <- logtiters1 - logtiters2

  result <- Hmisc::smean.cl.normal(logtiters, conf.int = level)
  list(
    mean_diff = unname(result["Mean"]),
    sd = sd(logtiters),
    mean_diff_lower = unname(result["Lower"]),
    mean_diff_upper = unname(result["Upper"]),
    mean_var = NA,
    sd_var = NA
  )

}


mean_titer_diffs_exclude_nd <- function(
  titers1,
  titers2,
  level = 0.95,
  dilution_stepsize
) {

  lessthans1 <- substr(titers1, 1, 1) == "<"
  morethans1 <- substr(titers1, 1, 1) == ">"
  titers1 <- titers1[!lessthans1 & !morethans1]

  logtiters1 <- log2(as.numeric(titers1) / 10)

  lessthans2 <- substr(titers2, 1, 1) == "<"
  morethans2 <- substr(titers2, 1, 1) == ">"
  titers2 <- titers2[!lessthans2 & !morethans2]

  logtiters2 <- log2(as.numeric(titers2) / 10)

  logtiters <- logtiters1 - logtiters2

  result <- Hmisc::smean.cl.normal(logtiters, conf.int = level)
  list(
    mean_diff = result[["Mean"]],
    sd = sd(logtiters),
    mean_diff_lower = result[["Lower"]],
    mean_diff_upper = result[["Upper"]],
    mean_var = NA,
    sd_var = NA
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
    mean_diff_upper = NA,
    mean_var = NA,
    sd_var = NA
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
        mean_diff_upper = result_ci["mean", 2],
        mean_var = result$vcov[1, 1],
        sd_var = result$vcov[2, 2]
      )
    })

  }, silent = TRUE)

  # Return output
  output

}

