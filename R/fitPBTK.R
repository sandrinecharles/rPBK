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
fitPBTK.stanPBTKdata <- function(stanPBTKdata, ...) {
  # remove additional variables
  dataFit <- stanPBTKdata
  dataFit$origin_data <- NULL
  dataFit$col_time <- NULL
  col_replicate <- NULL
  col_exposure  <- NULL
  col_compartment <- NULL
  stanfit <- rstan::sampling(stanmodels$PBTK, data = dataFit, ...)
  out <- list(stanPBTKdata = stanPBTKdata, stanfit = stanfit)
  class(out) <- append("fitPBTK", class(out))
  return(out)
}

#' @rdname fitPBTK
#'
#' @export
#'
fitPBTK_AD <- function(stanPBTKdata, ...) {
  # clean variables
  dataFit <- list(
    N_obs_comp = mData$N_obs_comp,
    N_rep = mData$N_rep,
    N_comp = mData$N_comp,
    time_obs_comp = mData$time_obs_comp,
    val_obs_comp = mData$val_obs_comp,
    t0 = mData$t0,
    tacc = mData$tacc,
    val_obs_exp = as.numeric(unique(mData$val_obs_exp))
  )
  stanfit <- rstan::sampling(stanmodels$PBTK_AD, data = dataFit, ...)
  out <- list(stanPBTKdata = stanPBTKdata, stanfit = stanfit)
  class(out) <- append("fitPBTK", class(out))
  return(out)
}


#' @rdname fitPBTK
#'
#' @export
#'
fitPBTK.stanPBTKoriginaldata <- function(stanPBTKdata, ...) {
  # remove additional variables
  dataFit <- stanPBTKdata
  dataFit$origin_data <- NULL
  dataFit$tp <- NULL
  stanfit <- rstan::sampling(stanmodels$PBTKoriginal, data = dataFit, ...)
  out <- list(stanPBTKdata = stanPBTKdata, stanfit = stanfit)
  class(out) <- append("fitPBTK", class(out))
  return(out)
}



