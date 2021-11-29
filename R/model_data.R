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
modelData_original <- function(object, time_accumulation, ...){


  data_total <- object
  #Data sets (all data sets have the same columns names and the same number if lines)
  data_caecum <- dplyr::filter(data_total, name == "caecum")
  data_cephalon <- dplyr::filter(data_total, name == "cephalon")
  data_intestin <- dplyr::filter(data_total, name == "intestin")
  data_reste <- dplyr::filter(data_total, name == "reste")

  #Creation of objects containing the each type of data
  #time corresponding to each observation
  vttot<-data_caecum$temps
  vt <- vttot[-(1:3)] #t=0 removed

  #observed data for each organ (C0 removed)
  Cobs_caecum <- data_caecum$concentration[-(1:3)]
  Cobs_cephalon <- data_cephalon$concentration[-(1:3)]
  Cobs_reste <- data_reste$concentration[-(1:3)]
  Cobs_intestin <- data_intestin$concentration[-(1:3)]

  #Value of each organ concentration at t=0 (mean of the 3 values)
  C0_caecum=mean(data_caecum[1:3,]$concentration)
  C0_cephalon=mean(data_cephalon[1:3,]$concentration)
  C0_reste=mean(data_reste[1:3,]$concentration)
  C0_intestin=mean(data_intestin[1:3,]$concentration)

  #p time step for Euler discretisation of ODE
  p=0.5
  tc=7 #duration of the accumulation phase
  td=14 # duration of the depuration phase
  Temps=seq(0,tc+td,by=p) #vector of time for both accumulation and depuration phases

  posTc=(tc/p)+1 #position of tc in the vector called Temps

  posTd=length(Temps)   #position of td in the vector called Temps (the last one)

  vtData=unique(vt) #time at which we have an observation

  #Contaminant concentration in water
  Cw=0.01108 #mean concentration in water during the accumulation phase
  vCw=c(rep(Cw,posTc),rep(0,(posTd-posTc))) # vector of concentration for each time of Temps vector (Cw during the accumultion phase and 0 during the depuration phase)

  #Position of each time of observation in the vector Temps
  Post2 <- vtData[1]/p+1
  Post5 <- vtData[2]/p+1
  Post7 <- vtData[3]/p+1
  Post9 <- vtData[4]/p+1
  Post15 <- vtData[5]/p+1
  Post21 <- vtData[6]/p+1

  rtrn_ls <- list(
    C0_caecum=C0_caecum,
    C0_cephalon=C0_cephalon,
    C0_intestin=C0_intestin,
    C0_reste=C0_reste,
    N_vCw = length(vCw),
    vCw=vCw,
    posTd=posTd,
    Post2=Post2,
    Post5=Post5,
    Post7=Post7,
    Post9=Post9,
    Post15=Post15,
    Post21=Post21,
    p=p,
    N=length(Cobs_caecum),
    Cobs_caecum=Cobs_caecum,
    Cobs_cephalon=Cobs_cephalon,
    Cobs_intestin=Cobs_intestin,
    Cobs_reste=Cobs_reste)

  class(rtrn_ls) <- append("stanPBTKoriginaldata", class(rtrn_ls))

  return(rtrn_ls)
}

