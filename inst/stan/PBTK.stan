//
// This Stan program defines a simple model, with a
// vector of values 'y' modeled as normally distributed
// with mean 'mu' and standard deviation 'sigma'.
//
// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//

// The input data is a vector 'y' of length 'N'.
functions {

#include /include/interpolation.stan

}
data {
  int<lower=0> N;// number of data per replicate
  int<lower=0> N_rep; // Number of replicate

  real tp[N];
  real tp_eval[N-1];
  matrix[N,N_rep] Cobs_caecum;
  matrix[N,N_rep] Cobs_cephalon;
  matrix[N,N_rep] Cobs_intestin;
  matrix[N,N_rep] Cobs_reste;

  real C0_caecum;
  real C0_cephalon;
  real C0_reste;
  real C0_intestin;

  real tacc;

  int<lower=0> N_Cw;
  vector[N_Cw] tp_Cw;
  vector[N_Cw] Cw;

}
transformed data{

  real t0;
  real y0[4];

  int<lower=0> x_int[1];

  real x_r[1+N_Cw+N_Cw];

  t0 = 0;
  y0[1] = C0_intestin ;
  y0[2] = C0_caecum ;
  y0[3] = C0_cephalon ;
  y0[4] = C0_reste ;

  x_int[1] = N_Cw ;

  x_r[1] = tacc ;
  for(i in 1:N_Cw){
    x_r[1+i] = tp_Cw[i] ;
  }
  for(i in 1:N_Cw){
    x_r[i+N_Cw+1] = Cw[i] ;
  }

}
parameters {
  real log10ku[4];
  real log10ke[4];

  real log10k1[4];
  real log10k2[4];
  real log10k3[4];
  real log10k4[4];

  real<lower=0> sigma[4];
}
transformed parameters{

  real<lower=0> theta[24] ;

  vector[N] Cpred_caecum;
  vector[N] Cpred_cephalon;
  vector[N] Cpred_reste;
  vector[N] Cpred_intestin;

  real<lower=0> y_sim[N-1, 4] ;

  for(i in 1:4){
    theta[i] = 10^log10ku[i] ;
    theta[i+4] = 10^log10ke[i] ;
    theta[i+8] = 10^log10k1[i] ;
    theta[i+12] = 10^log10k2[i] ;
    theta[i+16] = 10^log10k3[i] ;
    theta[i+20] = 10^log10k4[i] ;
  }

  y_sim = integrate_ode_rk45(ode_pbtk, y0, t0, tp_eval, theta, x_r, x_int) ;

  Cpred_intestin[1] = C0_intestin ;
  Cpred_caecum[1] = C0_caecum ;
  Cpred_cephalon[1] = C0_cephalon ;
  Cpred_reste[1] = C0_reste ;

  for(t in 1:N-1){
      Cpred_intestin[t+1] = y_sim[t,1] ;
      Cpred_caecum[t+1] = y_sim[t,2] ;
      Cpred_cephalon[t+1] = y_sim[t,3] ;
      Cpred_reste[t+1] = y_sim[t,4] ;
  }
}
model {

  target += uniform_lpdf(log10ku | -5, 5);
  target += uniform_lpdf(log10ke | -5, 5);

  target += uniform_lpdf(log10k1 | -5, 5);
  target += uniform_lpdf(log10k2 | -5, 5);
  target += uniform_lpdf(log10k3 | -5, 5);
  target += uniform_lpdf(log10k4 | -5, 5);

  target += gamma_lpdf(sigma | 0.001,0.001);

  for(rep in 1:N_rep){
    for(i in 1:N){
      target += normal_lpdf(Cobs_intestin[i,rep] | Cpred_intestin[i], sigma[1]);
      target += normal_lpdf(Cobs_caecum[i,rep] | Cpred_caecum[i], sigma[2]) ;
      target += normal_lpdf(Cobs_cephalon[i,rep] | Cpred_cephalon[i], sigma[3]);
      target += normal_lpdf(Cobs_reste[i,rep] | Cpred_reste[i], sigma[4]);
    }
  }
}
generated quantities {

  vector[N] Cgen_caecum;
  vector[N] Cgen_cephalon;
  vector[N] Cgen_reste;
  vector[N] Cgen_intestin;

  for(t in 1:N){
    Cgen_intestin[t] = normal_rng(Cpred_intestin[t], sigma[1]) ;
    Cgen_caecum[t] = normal_rng(Cpred_caecum[t], sigma[2]) ;
    Cgen_cephalon[t] = normal_rng(Cpred_cephalon[t], sigma[3]) ;
    Cgen_reste[t] = normal_rng(Cpred_reste[t], sigma[4]) ;
  }
}

