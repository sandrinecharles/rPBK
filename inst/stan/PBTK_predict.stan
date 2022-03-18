functions {

#include /include/interpolation.stan

}
data {
  // model structure
  int N_comp; // number of compartiment
  int N_exp; // number of exposure routes

  // Exposure
  int<lower=0> N_obs_exp;
  real time_obs_exp[N_obs_exp];
  real val_obs_exp[N_obs_exp, N_exp];

  // Compartiment
  int N_pred_comp ;
  real time_pred_comp[N_pred_comp] ;
  real C0_obs_comp[N_comp];

  // Time_accumulation ... TO BE REMOVED !!!
  real tacc;
  real t0;

  // DATA FOR SAMPLINGS
  int N_samples;
  real theta_array [N_samples, 2+N_comp, N_comp];
  matrix[N_samples, N_comp] sigma;
}
transformed data{

  real y0[N_comp] = C0_obs_comp;

  int<lower=0> x_int[3];

  real x_r_array[2+N_exp,N_obs_exp] ;
  real x_r[(2+N_exp)*N_obs_exp] ;

  real<lower=0> theta_matrix[N_samples, N_comp*(2+N_comp)];

  x_int[1] = N_comp ;
  x_int[2] = N_exp ;
  x_int[3] = N_obs_exp ;

  for(i in 1:N_obs_exp){
    x_r_array[1,i] = tacc ;
    x_r_array[2,i] = time_obs_exp[i] ;
    for(j in 1:N_exp){
      x_r_array[2+j,i] = val_obs_exp[i,j] ;
    }
  }
  x_r = to_array_1d(x_r_array) ;

  for(i in 1:N_samples){
    theta_matrix[i,1:N_comp*(2+N_comp)] = to_array_1d(theta_array[i,1:(2+N_comp),1:N_comp]) ;
  }

}
parameters {
}
model {
}
generated quantities {

  matrix[N_obs_exp, 2+N_exp] x_r_matrix = to_matrix(x_r, N_obs_exp, 2+N_exp) ;
  real tacc_test = x_r_matrix[1,1] ;
  vector[N_obs_exp] tp_Cw = to_vector(x_r_matrix[1:N_obs_exp, 2]) ;
  vector[N_obs_exp] Cw = to_vector(x_r_matrix[1:N_obs_exp, 3]) ;

  real y_sim[N_pred_comp, N_comp] ;
  real val_pred_comp[N_samples, N_pred_comp, N_comp];

  for(i in 1:N_samples) {

    y_sim = integrate_ode_rk45(ode_pbtk, y0, t0, time_pred_comp, theta_matrix[i,1:(N_comp*(2+N_comp))], x_r, x_int) ;

    for(i_comp in 1:N_comp){
      for(t in 1:N_pred_comp){
        val_pred_comp[i,t,i_comp] = normal_rng(y_sim[t,i_comp], sigma[i,i_comp]) ;
      }
    }
  }
}


