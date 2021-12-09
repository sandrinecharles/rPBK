#' Create a list giving data and parameters to use in the model inference.
#'
#' @param object An object of class \code{data.frame}
#' @param time_accumulation A scalar givin accumulation time
#' @param elimination_rate A scalar for the elimination rate. Default is \code{NA}.
#' To remove elimination rate, set \code{elimination_rate = 0}.
#' @param \dots Further arguments to be passed to generic methods
#'
#' @export
#'
#' @return A \code{list} with data and parameters require for model inference.
#'
#'
modelData <- function(object, ...){
  UseMethod("modelData")
}


#' @rdname modelData
#'
#' @export
#' @importFrom stats na.omit
modelData.data.frame <- function(object, time_accumulation=7, jitter_t0 = 0.01,...){

  data_total <- object
  #Data sets (all data sets have the same columns names and the same number if lines)
  data_caecum <- dplyr::filter(data_total, name == "caecum")
  data_cephalon <- dplyr::filter(data_total, name == "cephalon")
  data_intestin <- dplyr::filter(data_total, name == "intestin")
  data_reste <- dplyr::filter(data_total, name == "reste")

  #Creation of objects containing the each type of data
  #time corresponding to each observation
  vttot<-data_caecum$temps
  tp <- unique(vttot)
  tp[1] <- tp[1]+jitter_t0

  replicate = unique(data_total$replicat)
  N_rep = length(replicate)
  # replicate
  Cobs_caecum <- do.call(
    "cbind",
    lapply(1:N_rep, function(i) data_caecum[data_caecum$replicat==i,]$concentration)
  )
  Cobs_cephalon <- do.call(
    "cbind",
    lapply(1:N_rep, function(i) data_caecum[data_caecum$replicat==i,]$concentration)
  )
  Cobs_reste <- do.call(
    "cbind",
     lapply(1:N_rep, function(i) data_caecum[data_caecum$replicat==i,]$concentration)
  )
  Cobs_intestin <- do.call(
    "cbind",
     lapply(1:N_rep, function(i) data_caecum[data_caecum$replicat==i,]$concentration)
  )

  N = nrow(Cobs_caecum)
  N_Cw = N

  #Value of each organ concentration at t=0 (mean of the 3 values)
  C0_caecum=mean(data_caecum[1:3,]$concentration)
  C0_cephalon=mean(data_cephalon[1:3,]$concentration)
  C0_reste=mean(data_reste[1:3,]$concentration)
  C0_intestin=mean(data_intestin[1:3,]$concentration)

  #Contaminant concentration in water
  Cw_sct=0.01108 #mean concentration in water during the accumulation phase
  Cw=rep(Cw_sct, N_Cw) # vector of concentration for each time of Temps vector (Cw during the accumultion phase and 0 during the depuration phase)

  rtrn_ls <- list(
    N = N,
    tp = tp,
    tp_eval = tp[-1],
    N_rep = N_rep,
    C0_caecum=C0_caecum,
    C0_cephalon=C0_cephalon,
    C0_intestin=C0_intestin,
    C0_reste=C0_reste,
    N_Cw = N_Cw,
    tp_Cw = tp,
    Cw=Cw,
    tacc=time_accumulation,
    Cobs_caecum=Cobs_caecum,
    Cobs_cephalon=Cobs_cephalon,
    Cobs_intestin=Cobs_intestin,
    Cobs_reste=Cobs_reste)

  class(rtrn_ls) <- append("stanPBTKdata", class(rtrn_ls))

  return(rtrn_ls)
}


##########
#' @rdname modelData
#'
#' @export
#' @importFrom stats na.omit
modelData_original <- function(object, time_accumulation=7, ...){

  data_total <- object
  #Data sets (all data sets have the same columns names and the same number if lines)
  data_caecum <- dplyr::filter(data_total, name == "caecum")
  data_cephalon <- dplyr::filter(data_total, name == "cephalon")
  data_intestin <- dplyr::filter(data_total, name == "intestin")
  data_reste <- dplyr::filter(data_total, name == "reste")

  #Creation of objects containing the each type of data
  #time corresponding to each observation
  vttot<-data_caecum$temps
  tp <- unique(vttot)
  tp[1] <- tp[1]

  replicate = unique(data_total$replicat)
  N_rep = length(replicate)
  # replicate
  Cobs_caecum <- do.call(
    "cbind",
    lapply(1:N_rep, function(i) data_caecum[data_caecum$replicat==i,]$concentration)
  )
  Cobs_cephalon <- do.call(
    "cbind",
    lapply(1:N_rep, function(i) data_caecum[data_caecum$replicat==i,]$concentration)
  )
  Cobs_reste <- do.call(
    "cbind",
    lapply(1:N_rep, function(i) data_caecum[data_caecum$replicat==i,]$concentration)
  )
  Cobs_intestin <- do.call(
    "cbind",
    lapply(1:N_rep, function(i) data_caecum[data_caecum$replicat==i,]$concentration)
  )
  N = nrow(Cobs_caecum)

  #Value of each organ concentration at t=0 (mean of the 3 values)
  C0_caecum=mean(data_caecum[1:3,]$concentration)
  C0_cephalon=mean(data_cephalon[1:3,]$concentration)
  C0_reste=mean(data_reste[1:3,]$concentration)
  C0_intestin=mean(data_intestin[1:3,]$concentration)

  # p time step for Euler discretisation of ODE
  p <- 0.5

  tacc <- time_accumulation #duration of the accumulation phase
  td <- max(tp) -  time_accumulation# duration of the depuration phase

  # vector of time for both accumulation and depuration phases
  tp_Cw <- seq(0,tacc+td,by=p)
  N_Cw <- length(tp_Cw)

  # position of tc in the vector called Temps
  posTc <- match(tacc,tp_Cw)

  # Contaminant concentration in water
  # mean concentration in water during the accumulation phase
  Cw_cst=0.01108
  # vector of concentration for each time of Temps vector (Cw during the accumultion phase and 0 during the depuration phase)
  Cw=c(rep(Cw_cst,posTc),rep(0,(N_Cw-posTc)))

  # Position of each time of observation in the vector Temps
  tp_position <- match(tp,tp_Cw)

  rtrn_ls <- list(
    N = N,
    N_rep = N_rep,
    C0_caecum = C0_caecum,
    C0_cephalon = C0_cephalon,
    C0_intestin = C0_intestin,
    C0_reste = C0_reste,
    N_Cw = length(tp_Cw),
    tp_Cw = tp_Cw,
    tp_position = tp_position,
    p = p,
    N = length(Cobs_caecum),
    Cobs_caecum = Cobs_caecum,
    Cobs_cephalon = Cobs_cephalon,
    Cobs_intestin = Cobs_intestin,
    Cobs_reste = Cobs_reste)

  class(rtrn_ls) <- append("stanPBTKoriginaldata", class(rtrn_ls))

  return(rtrn_ls)
}

