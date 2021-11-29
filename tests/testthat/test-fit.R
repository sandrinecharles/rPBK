test_that("fit basic model", {


  library(pbtkDB)
  data("data_total")

  modelData_original <- modelData_original(data_total, time_accumulation = 7)
  fit_original <- fitPBTK(modelData_original)

  modelData <- modelData(data_total, time_accumulation = 7)
  fit <- fitPBTK(modelData, jitter_t0 = -0.5)

})
