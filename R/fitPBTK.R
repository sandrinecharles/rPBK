#' Bayesian inference of TK model with Stan
#'
#' @param stanPBTKdata List of Data require for computing
#' @param \dots Arguments passed to `rstan::sampling` (e.g. iter, chains).
#'
#' @return An object of class `fitPBTK` containing two object: \code{stanPBTKdata}
#' the data set used for inference and \code{stanfit}  returned by `rstan::sampling`
#'
#' @rdname fitPBTK
#'
#' @export
#'
fitPBTK <- function(stanPBTKdata, ...){
  UseMethod("fitPBTK")
}

#' @rdname fitPBTK
#'
#' @export
#'
fitPBTK.stanPBTKdata <- function(stanPBTKdata, ODE = FALSE, ...) {
  # if(ODE = TRUE){
  #   # remove additional variables
  #   dataFit <- stanPBTKdata
  #   dataFit$origin_data <- NULL
  #   dataFit$col_time <- NULL
  #   col_replicate <- NULL
  #   col_exposure  <- NULL
  #   col_compartment <- NULL
  #   stanfit <- rstan::sampling(stanmodels$PBTK, data = dataFit, ...)
  #   out <- list(stanPBTKdata = stanPBTKdata, stanfit = stanfit)
  #   class(out) <- append("fitPBTK", class(out))
  # } else{
    # clean variables
    x <- stanPBTKdata
    dataFit <- list(
      N_obs_comp = x$N_obs_comp,
      N_rep = x$N_rep,
      N_comp = x$N_comp,
      time_obs_comp = x$time_obs_comp,
      val_obs_comp = x$val_obs_comp,
      val_obs_exp = x$val_obs_exp[1],
      t0 = x$t0,
      tacc = x$tacc,
      ku_nest = x$ku_nest,
      ke_nest = x$ke_nest,
      k_nest = x$k_nest
    )
    stanfit <- rstan::sampling(stanmodels$PBTK_AD, data = dataFit, ...)
    out <- list(stanPBTKdata = stanPBTKdata, stanfit = stanfit)
    class(out) <- append("fitPBTK", class(out))
  # }

  return(out)
}



