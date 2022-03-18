test_that("sit ODE solver simple dataset", {

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
  # mData$C0_obs_comp = array(mData$C0_obs_comp, 1)
  fit_ode <- fitPBTK(mData, chains = 1, iter = 100)

  usethis::use_data(fit_ode)

  out <- rstan::extract(fit_ode$stanfit)
  # N_iter * N_time * N_comp
  dim(out$val_pred_comp)

  devtools::load_all()
  ls_out_init = lapply(1:mData$N_comp, function(i){
    df = .df_quant95(out$val_pred_comp[,,i])
    df$time = mData$time_obs_exp
    df$comp = mData$col_compartment[i]
    return(df)
  })

  df_out_init = do.call("rbind", ls_out_init)

  df_out =  df_out_init %>%
    tidyr::as_tibble()

  df_data = dataset %>%
    tidyr::pivot_longer(cols = mData$col_compartment)

  ggplot() +
    theme_minimal() +
    geom_ribbon(data = df_out,
                aes(x = time, ymin = qinf95, ymax=qsup95, fill = comp), alpha = 0.5) +
    geom_line(data = df_out,
              aes(x=time,y=q50,color=comp)) +
    geom_point(data = df_data,
               aes(x = time, y = value, color = name))

})


test_that("fit ODE solver", {

  library(pbtkDB)
  data("data_total")

  mData <- modelData(object = data_total,
            col_time = "temps",
            col_replicate = "replicat",
            col_exposure = "condition",
            col_compartment = c("intestin", "reste", "caecum", "cephalon"),
            time_accumulation = 7)

  fit_ode <- fitPBTK(mData, chains = 3, iter = 500)
  usethis::use_data(fit_ode)

  out <- rstan::extract(fit_ode$stanfit)
  # N_iter * N_time * N_comp
  dim(out$val_pred_comp)

  devtools::load_all()
  ls_out_init = lapply(1:mData$N_comp, function(i){
    df = .df_quant95(out$val_pred_comp[,,i])
    df$time = mData$time_obs_exp
    df$comp = mData$col_compartment[i]
    return(df)
  })

  df_out_init = do.call("rbind", ls_out_init)

  df_out =  df_out_init %>%
    tidyr::as_tibble()

  df_data = data_total %>%
    tidyr::pivot_longer(cols = mData$col_compartment)

  ggplot() +
    theme_minimal() +
    geom_ribbon(data = df_out,
                aes(x = time, ymin = qinf95, ymax=qsup95, fill = comp), alpha = 0.5) +
    geom_line(data = df_out,
              aes(x=time,y=q50,color=comp)) +
    geom_point(data = df_data,
               aes(x = temps, y = value, color = name))

  # A faire cf: morse / rbioacc
  ppc(fit_ode)
  prior_post(fit_ode)
  predict(fit_ode)

  # modelData <- modelData(data_total, time_accumulation = 7, jitter_t0 = 0)
  # modelData <- modelData(data_total, time_accumulation = 7)
  # # fit_ode <- fitPBTK(modelData, chains = 1, iter = 10, algorithm = "Fixed_param")
  # fit_ode <- fitPBTK(modelData, chains = 1, iter = 10, algorithm = "HMC")

})


