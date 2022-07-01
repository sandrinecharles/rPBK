#' Plotting method for \code{fitPBTK} objects
#'
#' This is the generic \code{plot} S3 method for the
#' \code{fitTK}.  It plots the fit obtained for each
#' variable in the original dataset.
#'
#' @param x And object returned by \code{fitPBTK}
#' @param \dots Additional arguments
#'
#' @return a plot of class \code{ggplot}
#'
#' @export
#'
#' @import ggplot2
#'
plot.fitPBTK <- function(x, ...){

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

  # plot
  plt <- ggplot(data = df_fit) +
    theme_classic() +
    labs(x = "Time", y = "Concentration") +
    scale_y_continuous(limits = c(0,NA)) +
    geom_ribbon(
      aes_string(x = 'time', ymin = 'qinf95', ymax = 'qsup95'),
      fill = "grey80") +
    geom_line(
      aes_string(x = 'time', y = 'q50'),
      color = "orange") +
    geom_point(data = df_data,
               aes_string(x = 'time', y = 'observation' )) +
    facet_wrap(~compartment, scales = "free")

  return(plt)
}



############ INTERNAL


# @param x An object of class \code{fitPBTK}
#
# @importFrom stats quantile
#
# @export
#
.df_quant95 <- function(x,...){

  mat_quants = apply(x, 2, quantile, probs = c(0.025, 0.5, 0.975), ...)

  df = data.frame(
    q50 = mat_quants[2,],
    qinf95 = mat_quants[1,],
    qsup95 = mat_quants[3,]
  )
  return(df)
}
