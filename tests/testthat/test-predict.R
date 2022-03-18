
test_that("multiplication works", {

  library(pbtkDB)
  ######### SIMPLE



  ######## Complex
  data("data_total")
  d <- modelData(object = data_total,
                 col_time = "temps",
                 col_replicate = "replicat",
                 col_exposure = "condition",
                 col_compartment = c("intestin", "caecum", "cephalon", "reste"),
                 time_accumulation = 7)

  # predict
  theta_matrix = matrix(c(
    c(1917,1571,91.1,135), # ku
    c(.506,.06,.053,.026), # ke
    rep(0,4*4)), # k
    byrow = FALSE, nrow = 4
  )

  N_samples = 2
  theta_array <- array(NA, dim=c(N_samples, 2+d$N_comp, d$N_comp))
  theta_array[1,,] <- t(theta_matrix)
  theta_array[2,,] <- t(theta_matrix)


  dataFit = list(
    N_comp = d$N_comp,
    N_exp = d$N_exp,

    N_obs_exp = d$N_obs_exp,
    time_obs_exp = d$time_obs_exp,
    val_obs_exp = d$val_obs_exp,

    N_pred_comp = 100,
    time_pred_comp = seq(0.01, 21, length.out = 100),
    C0_obs_comp = d$C0_obs_comp,

    t0 = 0,
    tacc = d$tacc,

    theta_array = array(t(theta_matrix), dim=c(N_samples, d$N_comp+2, d$N_comp)),

    N_samples = N_samples,
    # sigma = array(c(11,16,1.0,0.95), dim=c(N_samples, d$N_comp))
    sigma = array(c(1e-6,1e-6,1e-6,1e-6), dim=c(N_samples, d$N_comp))
  )
  devtools::load_all()
  stanfit <- rstan::sampling(stanmodels$PBTK_predict,
                             data = dataFit,
                             algorithm = "Fixed_param",
                             chains = 1, iter = N_samples, warmup = 0, verbose = TRUE)
  out = rstan::extract(stanfit)

  # N_chain * N_sample * N_pred_time * N_comp
  out$Cgen_comp[1,1,,]
  out$Cgen_comp[1,,,]


  dim(out$Cgen_comp[1,1,,])


  # N_sample * N_time
  m = matrix(1:9, ncol = 3)
  .df_quant95(m)
  # N_time * 3
  quantile(m[,1], c(0.025,0.5,0.975))



  .df_quant95()

  out$Cgen_comp[,,,3]

df_out =  out$Cgen_comp[1,,,] %>%
    tidyr::as_tibble() %>%
    dplyr::mutate(time = dataFit$time_pred_comp) %>%
    tidyr::pivot_longer(cols = -"time")

ggplot() +
  theme_minimal() +
  geom_line(data = df_out,
            aes(x=time,y=value,color=name))

})
