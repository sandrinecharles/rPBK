test_that("multiplication works", {

  library(pbtkDB)
  data("Male_Gammarus_Single")
  DS <- Male_Gammarus_Single
  # work only when replicate have the same length !!!
  dataset <- DS[DS$replicate == 1,]
  mData <- modelData(object = dataset,
                     col_time = "time",
                     col_replicate = "replicate",
                     col_exposure = "expw",
                     col_compartment = "conc",
                     time_accumulation = 4)

  fitPBTK_out <- fitPBTK_AD(mData)

  data("data_total")
  mData <- modelData(object = data_total,
                     col_time = "temps",
                     col_replicate = "replicat",
                     col_exposure = "condition",
                     col_compartment = c("intestin", "reste", "caecum", "cephalon"),
                     time_accumulation = 7)


  start_time <- Sys.time()
  fitPBTK_out <- fitPBTK_AD(mData, chains = 3, iter = 10000)
  end_time <- Sys.time()
  end_time - start_time

  plot_AD(fitPBTK_out)


  x = fitPBTK_out

  out_fit <- rstan::extract(x$stanfit)
  out_data <- x$stanPBTKdata

  df_ke <- out_fit$ke
  df_ku <- out_fit$U
  colnames(df_ke) <- out_data$col_compartment
  colnames(df_ku) <- out_data$col_compartment

  df_k <- out_fit$log10k

}

