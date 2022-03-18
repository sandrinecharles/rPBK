// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//
functions {

#include /include/interpolation.stan

}
data {
  int print_messages;
  int with_k;

  int<lower=0> N_obs_comp;// number of data per replicate
  int<lower=0> N_rep; // Number of replicate
  int<lower=0> N_comp; // Number of compartiment

  real time_obs_comp[N_obs_comp];
  real time_eval[N_obs_comp-1];

  real val_obs_comp[N_obs_comp,N_rep,N_comp];
  real C0_obs_comp[N_comp];

  real t0;

  real tacc;

  int<lower=0> N_obs_exp;
  int<lower=0> N_exp;
  vector[N_obs_exp] time_obs_exp;
  real val_obs_exp[N_obs_exp, N_exp];

}
transformed data{

  real y0[N_comp] = C0_obs_comp;

  int<lower=0> x_int[3];

  real x_r_array[2+N_exp, N_obs_exp] ;
  real x_r[(2+N_exp)*N_obs_exp] ;

  x_int[1] = N_comp;
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
}
parameters {
  real log10ku[N_comp];
  real log10ke[N_comp];

  real log10k[N_comp,N_comp];

  real<lower=0> sigma[N_comp];
}
transformed parameters{

  real<lower=0> theta[N_comp*(2+N_comp)];

  real Cpred_comp[N_obs_comp, N_comp];
  real y_sim[N_obs_comp-1, N_comp] ;

  matrix[2+N_comp, N_comp] theta_matrix ;
  for(i in 1:N_comp){
    theta_matrix[1,i] = 10^log10ku[i] ;
    theta_matrix[2,i] = 10^log10ke[i] ;
    for(j in 1:N_comp){
      if(with_k == 1){
        theta_matrix[2+j,i] = 10^log10k[i,j] ;
      } else{
         theta_matrix[2+j,i] = 0 ;
      }
    }
  }
  theta = to_array_1d(theta_matrix) ;

  if(print_messages == 1){
    print("theta_matrix", theta_matrix);
    print("theta", theta);
    print("x_r", x_r);
    print("x_int", x_int);
  }

  y_sim = integrate_ode_rk45(ode_pbtk, y0, t0, time_eval, theta, x_r, x_int) ;

  for(i_comp in 1:N_comp){
    Cpred_comp[1,i_comp] = C0_obs_comp[i_comp];
    for(t in 1:N_obs_comp-1){
      Cpred_comp[t+1,i_comp] = y_sim[t,i_comp] ;
    }
  }
}
model {

  target += uniform_lpdf(log10ku | -5, 5);
  target += uniform_lpdf(log10ke | -5, 5);

  for(i_comp in 1:N_comp){
    target += uniform_lpdf(log10k[1:N_comp, i_comp] | -5, 5);
  }

  target += gamma_lpdf(sigma | 0.01,0.01);

  for(i_rep in 1:N_rep){
    for(i in 1:N_obs_comp){
      for(i_comp in 1:N_comp){
        target += normal_lpdf(val_obs_comp[i,i_rep,i_comp] | Cpred_comp[i,i_comp], sigma[i_comp]);
      }
    }
  }
}
generated quantities {
  // change Cgen_comp to val_gen_comp
  real  val_pred_comp[N_obs_comp, N_comp];

  for(i_comp in 1:N_comp){
    for(t in 1:N_obs_comp){
      val_pred_comp[t, i_comp] = normal_rng(Cpred_comp[t, i_comp], sigma[i_comp]) ;
    }
  }
}

