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

  fit <- x

  df <- .df_for_plot(fit)

  # HACK TO BE > 0
  df$q50 <- ifelse(df$q50<0,0,df$q50)
  df$qinf95 <- ifelse(df$qinf95<0,0,df$qinf95)
  df$qsup95 <- ifelse(df$qsup95<0,0,df$qsup95)

  plt <- ggplot(data = df) +
    theme_classic() +
    labs(x = "Time", y = "Concentration") +
    # scale_y_continuous(limits = c(0,NA)) +
    geom_ribbon(
      aes_string(x = 'time', ymin = 'qinf95', ymax = 'qsup95'), fill = "grey80") +
    geom_line(aes_string(x = 'time', y = 'q50'), color = "orange") +
    geom_point(aes_string(x = 'time', y = 'observation' )) +
    facet_wrap(~variable, scales = "free")

  return(plt)
}

############ INTERNAL

#' @param x An object of class \code{fitPBTK}
#'
#' @importFrom stats quantile
#'
#' @export
#'
.df_quant95 <- function(x,...){

  mat_quants = apply(x, 2, quantile, probs = c(0.025, 0.5, 0.975), ...)

  df = data.frame(
    q50 = mat_quants[2,],
    qinf95 = mat_quants[1,],
    qsup95 = mat_quants[3,]
  )
  return(df)
}


.add_data = function(df_quant95,tp,data,id){
  if(is.vector(data)){
    df_quant95$time = tp
    df_quant95$observation = data
    df_quant95$replicate = 1
    df <- df_quant95
  } else{
    ls <- lapply(1:ncol(data),
                 function(i){
                   df_quant95$time = tp
                   df_quant95$observation = data[,i]
                   df_quant95$replicate = i
                   return(df_quant95)
                 })
    df <- do.call("rbind", ls)
  }
  df <- df[df$observation != Inf,]
  df$variable <-  id
  return(df)
}



.df_for_plot <- function(fit){
  fitMCMC <- rstan::extract(fit$stanfit)
  data <- fit$stanPBTKdata
  #
  ls_out <- list()
  ls_out$intestin  <- .add_data(
    df_quant95 = .df_quant95(fitMCMC$Cgen_intestin) ,
    tp = data$tp,
    data = data$Cobs_intestin ,
    id = "intestin "
  )
  ls_out$caecum <- .add_data(
    df_quant95 = .df_quant95(fitMCMC$Cgen_caecum),
    tp = data$tp,
    data = data$Cobs_caecum,
    id = "caecum"
  )
  ls_out$cephalon <- .add_data(
    df_quant95 = .df_quant95(fitMCMC$Cgen_cephalon),
    tp = data$tp,
    data = data$Cobs_cephalon,
    id = "cephalon"
  )
  ls_out$reste <- .add_data(
    df_quant95 = .df_quant95(fitMCMC$Cgen_reste),
    tp = data$tp,
    data = data$Cobs_reste,
    id = "reste"
  )

  df <- do.call("rbind", ls_out)

  return(df)
}

