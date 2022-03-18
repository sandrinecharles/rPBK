library(rstan)
library(dplyr)
library(tidyr)
library(ggplot2)

##################

model_ode <-
  '
  functions{

  int findfirst(real t, vector xt) {
  int i = 0 ;
  if(t == max(xt)){
    i = num_elements(xt) - 1 ;
    return i ;
  } else if(t < min(xt) || t > max(xt)){
    return i ;
  } else {
    while (t >= xt[i+1]){
      i = i+1 ;
    }
    return i ;
  }
}

real interpolate(real x, vector xpt, vector ypt){
  if(x >= min(xpt) && x <= max(xpt)){
    int idx = findfirst(x, xpt) ;
    return ypt[idx] + (x - xpt[idx]) * (ypt[idx+1] - ypt[idx]) / (xpt[idx+1] - xpt[idx]) ;
  } else{
    return 0.0 ;
  }
}


real[] ode_pbtk( real t,      // time
                 real[] y,    // variables
                 real[] theta,
                 real[] x_r,
                 int[] x_i) {

  // recover dimension
  int N_comp = x_i[1] ;
  int N_exp = x_i[2] ;
  int N_tp_exp = x_i[3] ;

  // recover parameters
  matrix[N_comp, N_comp+2] theta_matrix = to_matrix(theta, N_comp, N_comp+2) ;
  real ku[N_comp] = to_array_1d(theta_matrix[1:N_comp,1]) ;
  real ke[N_comp] = to_array_1d(theta_matrix[1:N_comp,2]) ;
  real k[N_comp, N_comp] = to_array_2d(theta_matrix[1:N_comp, 3:(N_comp+2)]) ;

  // recover exposure profiles
  matrix[N_tp_exp, 2+N_exp] x_r_matrix = to_matrix(x_r, N_tp_exp, 2+N_exp) ;
  real tacc = x_r_matrix[1,1] ;
  vector[N_tp_exp] tp_Cw = to_vector(x_r_matrix[1:N_tp_exp,2]) ;
  vector[N_tp_exp] Cw = to_vector(x_r_matrix[1:N_tp_exp,3]) ;

  // vector[1+n_met] dydt ;
  real dydt[N_comp] ;

  // Diagonals are 0
  for(i_comp in 1:N_comp){
    k[i_comp,i_comp] = 0 ;
  }

  if(t <= tacc){
    // Accumulation
    for(i in 1:N_comp){
       dydt[i] = ku[i] * interpolate(t, tp_Cw, Cw) - ke[i] * y[i] +
        to_row_vector(k[i,1:N_comp]) * to_vector(y[1:N_comp]) -
        sum(k[1:N_comp,i]) * y[i] ;
    }
  } else{
    // Depuration
    for(i in 1:N_comp){
       dydt[i] = - ke[i] * y[i] +
       to_row_vector(k[i,1:N_comp]) * to_vector(y[1:N_comp]) -
       sum(k[1:N_comp,i]) * y[i] ;
    }
  }
  return(dydt) ;
 }
real[,] ode_solver(real t0,      // time
                   real[] y0,    // variables
                   real[] time_eval,
                   real[] theta,
                   real[] x_r,
                   int[] x_i){
 return integrate_ode_rk45(ode_pbtk, y0, t0, time_eval, theta, x_r, x_i) ;
}
  }
'
expose_stan_functions(stanc(model_code = model_ode))

sol = ode_solver(
  0,
  c(.806,.647,0.067,0.048),
  seq(0, 21, length.out = 100),
  c(c(1917,1571,91.1,135), c(.506,.06,.053,.026), rep(0,16)),
  c(rep(7,7),c(0,2,5,7,9,15,21),rep(0.01108,7)),
  c(4,1,7))

df_out = do.call("rbind", sol) %>%
  tidyr::as_tibble() %>%
  dplyr::mutate(time = d$time_obs) %>%
  tidyr::pivot_longer(cols = -"time")

ggplot() +
  theme_minimal() +
  geom_line(data = df_out,
            aes(x=time,y=value,color=name))
