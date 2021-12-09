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
  stanfit <- rstan::sampling(stanmodels$PBTK, data = dataFit, ...)
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
  stanfit <- rstan::sampling(stanmodels$PBTKoriginal, data = dataFit, ...)
  out <- list(stanPBTKdata = stanPBTKdata, stanfit = stanfit)
  class(out) <- append("fitPBTK", class(out))
  return(out)
}


#' Prediction using Stan generated quantity simulator
#'
#' @param fitPBTK An object of class \code{fitPBTK}
#' @param \dots Supplementary arguments
#'
#' @return An object of class `predictPBTK`
#'
#' @rdname predictPBTK
#'
#' @export
#'
predictPBTK <- function(fitPBTK, ...){
  UseMethod("predictPBTK")
}


#' @rdname predictPBTK
#'
#' @export
#'
predictPBTK.fitPBTK <- function(fitPBTK, ...) {

  # remove additional variables
  dataFit <- .foo(fitPBTK)

  stanfit <- rstan::sampling(stanmodels$PBTK_predict, data = dataFit, algorithm = "Fixed_param", ...)
  out <- list(dataFit = dataFit, stanfit = stanfit)
  class(out) <- append("predictPBTK", class(out))
  return(out)
}

######## INTERNAL

.foo <- function(fitPBTK){
  return(ls= list())
}
