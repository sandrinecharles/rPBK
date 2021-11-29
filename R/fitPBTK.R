#' Posterior predictive check
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


#' Bayesian inference of TK model with Stan
#'
#' @rdname fitPBTK
#'
#' @export
#'
fitPBTK.stanPBTKdata <- function(stanPBTKdata, ...) {
  # remove additional variables
  dataFit <- stanPBTKdata
  dataFit$origin_data <- NULL
  stanfit <- rstan::sampling(stanmodels$PBTK, data = dataFit, ...)
  out <- list(stanPBTKdata = stanPBTKdata, stanfit = stanfit)
  class(out) <- append("fitPBTK", class(out))
  return(out)
}

#' Bayesian inference of TK model with Stan
#'
#' @rdname fitPBTK
#'
#' @export
#'
fitPBTK.stanPBTKoriginaldata <- function(stanPBTKdata, ...) {
  # remove additional variables
  dataFit <- stanPBTKdata
  dataFit$origin_data <- NULL
  stanfit <- rstan::sampling(stanmodels$PBTKoriginal, data = dataFit, ...)
  out <- list(stanPBTKdata = stanPBTKdata, stanfit = stanfit)
  class(out) <- append("fitPBTK", class(out))
  return(out)
}
