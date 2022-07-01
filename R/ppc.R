#' Posterior predictive check plot
#'
#' Plots posterior predictive check for \code{fitPBTK}
#'
#' @rdname PPC
#'
#' @param x an object used to select a method \code{ppc}
#' @param \dots Further arguments to be passed to generic methods
#'
#' @export
ppc <- function(x, ...){
  UseMethod("ppc")
}

#' @rdname PPC
#'
#' @export
#'
#'
ppc.fitPBTK <- function(x, ...){

  out_fit <- rstan::extract(x$stanfit)
  out_data <- x$stanPBTKdata

  # data.frame for observation
  ls_data <- lapply(1:out_data$N_comp, function(i_comp){
    df = data.frame(observation = c(out_data$val_obs_comp[,,i_comp]))
    # df$replicate = rep(1:out_data$N_rep, each = out_data$N_obs_comp)
    df$time = rep(out_data$time_obs_comp, out_data$N_rep)
    df$compartment = out_data$col_compartment[i_comp]
    return(df)
  })
  df_data = do.call("rbind", ls_data)

  # data.frame for prediction
  Cpred_quant <- lapply(1:out_data$N_comp, function(i_comp){
    df = .df_quant95(out_fit$Cpred_comp[,,i_comp])
    df$time = out_data$time_obs_comp
    df$compartment = out_data$col_compartment[i_comp]
    return(df)
  })
  df_fit <- do.call("rbind", Cpred_quant)

  df_ppc <- merge(x = df_data, y = df_fit, by = c("time", "compartment"), all = TRUE)
  df_ppc$color <- with(df_ppc, ifelse(qinf95 > observation | qsup95 < observation, "out", "in"))

  plt <- ggplot(data = df_ppc) +
    theme_classic() +
    labs(x = "Observation", y = "Prediction") +
    theme(legend.position="none") +
    scale_colour_manual(values = c("green", "red")) +
    geom_abline(slope = 1) +
    geom_linerange(
      aes(x = observation,ymin = qinf95, ymax = qsup95, color = color),
                   position = position_dodge(width=0.5)) +
    geom_point(aes(x = observation, y = q50),
               position = position_dodge(width=0.5))

  return(plt)
}
