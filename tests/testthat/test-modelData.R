test_that("multiplication works", {

  library(pbtkDB)
  data("data_total")

  modelData_original <- modelData_original(data_total, time_accumulation = 7)

  modelData <- modelData(data_total, time_accumulation = 7)

})
