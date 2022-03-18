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
#'
#'
modelData.data.frame <- function(
  object,
  col_time = NA,
  col_replicate = NA,
  col_exposure = NA,
  col_compartment = NA,
  time_accumulation = NA,
  print_messages = FALSE,
  with_k = 1,
  ...){

  print_messages <- ifelse(print_messages == FALSE, 0, 1)

  data_total <- object

  # recover time objects
  if(!is.na(col_time)){
    data_time = data_total[[col_time]]
  } else{
    data_time = data_total$time
  }
  uniq_time <- sort(unique(data_time))
  N_time <- length(uniq_time)

  # recover replicate objects
  if(!is.na(col_replicate)){
    data_replicate = data_total[[col_replicate]]
  } else{
    data_replicate = data_total$replicate
  }
  uniq_replicate <- unique(data_replicate)
  N_replicate <- length(uniq_replicate)

  # recover exposure object
  N_exposure <- length(col_exposure)
  ls_exposure <- lapply(1:N_exposure, function(i_exp){
    r <- do.call(
      "cbind",
      lapply(1:N_replicate, function(i) data_total[data_replicate==uniq_replicate[1],col_exposure[i_exp]])
    )
  })
  Cobs_exposure <- array(unlist(ls_exposure), dim=c(N_time, N_exposure))

  # recover compartment object
  N_compartment <- length(col_compartment)
  ls_compartment <- lapply(1:N_compartment, function(i_comp){
    r <- do.call(
      "cbind",
      lapply(1:N_replicate, function(i) data_total[data_replicate==uniq_replicate[i],col_compartment[i_comp]])
    )
  })
  Cobs_compartment <- array(unlist(ls_compartment), dim=c(N_time, N_replicate, N_compartment))
  #Value of each organ concentration at t=0 (mean of the 3 values)
  Cobs_compartment_t0 <- sapply(1:N_compartment, function(i) mean(Cobs_compartment[1,1:N_replicate,i]))
  # to give dimension with a scalare
  Cobs_compartment_t0 <- array(Cobs_compartment_t0, N_compartment)

  # Creation of objects containing the each type of data
  # time corresponding to each observation
  tp <- uniq_time

  rtrn_ls <- list(
    origin_data = object,
    col_time = col_time,
    col_replicate = col_replicate,
    col_exposure = col_exposure,
    col_compartment = col_compartment,
    # --------
    N_obs_comp = N_time,
    time_obs_comp = tp,
    time_eval = tp[-1],
    N_rep = N_replicate,
    N_exp = N_exposure,
    N_comp = N_compartment,
    val_obs_comp = Cobs_compartment,
    C0_obs_comp = Cobs_compartment_t0,
    val_obs_exp = Cobs_exposure,
    N_obs_exp = N_time,
    time_obs_exp = tp,
    t0 = 0,
    tacc = time_accumulation,
    print_messages = print_messages,
    with_k = with_k)

  class(rtrn_ls) <- append("stanPBTKdata", class(rtrn_ls))

  return(rtrn_ls)
}


##########
#' @rdname modelData
#'
#' @export
#' @importFrom stats na.omit
modelData_original <- function(
  object,
  col_time = NA,
  col_replicate = NA,
  col_exposure = NA,
  col_compartment = NA,
  time_accumulation = NA, ...){

  data_total <- object

  # recover time objects
  if(!is.na(col_time)){
    data_time = data_total[[col_time]]
  } else{
    data_time = data_total$time
  }
  uniq_time <- sort(unique(data_time))
  N_time <- length(uniq_time)

  # recover replicate objects
  if(!is.na(col_replicate)){
    data_replicate = data_total[[col_replicate]]
  } else{
    data_replicate = data_total$replicate
  }
  uniq_replicate <- unique(data_replicate)
  N_replicate <- length(uniq_replicate)

  # recover exposure object
  N_exposure <- length(col_exposure)
  ls_exposure <- lapply(1:N_exposure, function(i_exp){
    r <- do.call(
      "cbind",
      lapply(1:N_replicate, function(i) data_total[data_replicate==uniq_replicate[1],col_exposure[i_exp]])
    )
  })
  Cobs_exposure <- array(unlist(ls_exposure), dim=c(N_time, N_exposure))

  # recover compartment object
  N_compartment <- length(col_compartment)
  ls_compartment <- lapply(1:N_compartment, function(i_comp){
    r <- do.call(
      "cbind",
      lapply(1:N_replicate, function(i) data_total[data_replicate==uniq_replicate[i],col_compartment[i_comp]])
    )
  })
  Cobs_compartment <- array(unlist(ls_compartment), dim=c(N_time, N_replicate, N_compartment))
  #Value of each organ concentration at t=0 (mean of the 3 values)
  Cobs_compartment_t0 <- sapply(1:N_compartment, function(i) mean(Cobs_compartment[1,1:N_replicate,i]))

  # p time step for Euler discretisation of ODE
  p <- 0.5

  tacc <- time_accumulation #duration of the accumulation phase
  td <- max(uniq_time) -  time_accumulation# duration of the depuration phase

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
  tp_position <- match(uniq_time,tp_Cw)

  rtrn_ls <- list(
    N_time = N_time,
    N_rep = N_replicate,
    N_comp = N_compartment,
    C0_comp = Cobs_compartment_t0,
    N_Cw = length(tp_Cw),
    tp_Cw = tp_Cw,
    Cw = Cw,
    tp_position = tp_position,
    tp = uniq_time,
    tacc = time_accumulation,
    p = p,
    Cobs_comp = Cobs_compartment)

  class(rtrn_ls) <- append("stanPBTKoriginaldata", class(rtrn_ls))

  return(rtrn_ls)
}

