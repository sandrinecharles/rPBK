// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//
functions {

#include /include/interpolation.stan

}
data {
  int<lower=0> N_obs_comp;// number of data per replicate
  int<lower=0> N_rep; // Number of replicate
  int<lower=0> N_comp; // Number of compartiment

  vector[N_obs_comp] time_obs_comp;

  vector[N_comp] ke_nest;
  vector[N_comp] ku_nest;
  matrix[N_comp,N_comp] k_nest;

  real val_obs_comp[N_obs_comp,N_rep,N_comp];

  real t0;

  real tacc;

  real val_obs_exp;

}
parameters {
  real log10ku[N_comp];
  real log10ke[N_comp];

  real log10k[N_comp,N_comp];

  real<lower=0> sigma[N_comp];
}
transformed parameters{

  matrix[N_obs_comp, N_comp] Cpred_comp;

  vector[N_comp] ku ;
  vector[N_comp] ke ;
  matrix[N_comp, N_comp] k ;

  matrix[N_comp, N_comp] E ;
  matrix[N_comp, N_comp] I ;

  for(i in 1:N_comp){
    ku[i] = ku_nest[i] == 0 ? 0 :  10^log10ku[i] ;
    ke[i] = ke_nest[i] == 0 ? 0 :  10^log10ke[i] ;
    for(j in 1:N_comp){
      k[i,j] = k_nest[i,j] == 0 ? 0 :  10^log10k[i,j] ;
    }
  }

  E = matrix_E(ke,k,N_comp);
  I = matrix_I(N_comp);

  Cpred_comp = exact_AD_long(
    time_obs_comp, // vector time,
    N_obs_comp, // int N_time,
    tacc, // real tacc,
    E, // matrix E,
    I, // matrix I,
    ku, // vector U,
    N_comp, // int N_comp,
    val_obs_exp // real Cx
    ) ;

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
  real val_pred_comp[N_obs_comp, N_comp];
  real log_lik[N_obs_comp, N_rep, N_comp];

  for(i_comp in 1:N_comp){
    for(t in 1:N_obs_comp){
      val_pred_comp[t, i_comp] = normal_rng(Cpred_comp[t, i_comp], sigma[i_comp]) ;
    }
  }
  // Log likelihood
  for(i_rep in 1:N_rep){
    for(i in 1:N_obs_comp){
      for(i_comp in 1:N_comp){
        log_lik[i,i_rep,i_comp] = normal_lpdf(val_obs_comp[i,i_rep,i_comp] | Cpred_comp[i,i_comp], sigma[i_comp]);
      }
    }
  }
}


