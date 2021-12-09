test_that("fit basic model", {


  library(pbtkDB)
  data("data_total")

  modelData_original <- modelData_original(data_total, time_accumulation = 7)
  fit_original <- fitPBTK(modelData_original, chains = 4, iter = 1000)

  modelData <- modelData(data_total, time_accumulation = 7, jitter_t0 = 0.001)
  # fit_ode <- fitPBTK(modelData, chains = 1, iter = 10, algorithm = "Fixed_param")
  fit_ode <- fitPBTK(modelData, chains = 1, iter = 10, algorithm = "HMC")

})

