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


#' Interpolate function implemented in Stan only export for checking
#'
#' @export
export_interpolate <- function(x, xpt, ypt, chain = 1, iter =1, ...) {
  if(length(xpt) != length(ypt)) stop("length xpt and ypt mismatch")
  data = list(x=x,xpt=xpt,ypt=ypt,N=length(xpt))
  stanfit <- rstan::sampling(stanmodels$export_interpolate, data = data,
                             algorithm = "Fixed_param", chain = chain, iter = iter, ...)
  out <- rstan::extract(stanfit)
  return(out$y)
}
