
test_that("mean titer", {
  
  titers <- c("40", "20", "80")
  result <- mean_titers(
    titers = titers,
    method = "truncated_normal",
    dilution_stepsize = 0
  )
  
  expect_equal(
    round(result$mean, 3),
    2
  )
  
  expect_false(is.na(result$sd))
  
})


test_that("mean titer with fixed sd", {
  
  titers <- c("40", "20", "80")
  result <- mean_titers(
    titers = titers,
    method = "truncated_normal",
    dilution_stepsize = 0,
    sd = 2
  )
  
  expect_equal(
    result$mean,
    2
  )
  
  expect_equal(
    result$sd,
    2
  )
  
})
