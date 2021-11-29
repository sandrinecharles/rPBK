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

  x_int[1] = N_Cw ;

  t0 = 0;
  y0[1] = C0_intestin ;
  y0[2] = C0_caecum ;
  y0[3] = C0_cephalon ;
  y0[4] = C0_reste ;

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

  real theta[24] ;

  matrix[N,N_rep] Cpred_caecum;
  matrix[N,N_rep] Cpred_cephalon;
  matrix[N,N_rep] Cpred_reste;
  matrix[N,N_rep] Cpred_intestin;

  real y_sim[N, 4] ;

  for(i in 1:4){
    theta[i] = 10^log10ku[i] ;
    theta[i+4] = 10^log10ke[i] ;
    theta[i+8] = 10^log10k1[i] ;
    theta[i+12] = 10^log10k2[i] ;
    theta[i+16] = 10^log10k3[i] ;
    theta[i+20] = 10^log10k4[i] ;
  }

  for(rep in 1:N_rep){
    y_sim = integrate_ode_rk45(ode_pbtk, y0, t0, tp, theta, x_r, x_int) ;

    for(t in 1:N){
        Cpred_intestin[t,rep] = y_sim[t,1] ;
        Cpred_caecum[t,rep] = y_sim[t,2] ;
        Cpred_cephalon[t,rep] = y_sim[t,3] ;
        Cpred_reste[t,rep] = y_sim[t,4] ;
    }
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
      Cobs_intestin[i,rep] ~ normal(Cpred_intestin[i,rep], sigma[1]);
      Cobs_caecum[i,rep] ~ normal(Cpred_caecum[i,rep], sigma[2]);
      Cobs_cephalon[i,rep] ~ normal(Cpred_cephalon[i,rep], sigma[3]);
      Cobs_reste[i,rep] ~ normal(Cpred_reste[i,rep], sigma[4]);
    }
  }

}
generated quantities {

  matrix[N,N_rep] Cgen_caecum;
  matrix[N,N_rep] Cgen_cephalon;
  matrix[N,N_rep] Cgen_reste;
  matrix[N,N_rep] Cgen_intestin;
  for(rep in 1:N_rep){
    for(t in 1:N){
      Cgen_intestin[t,rep] = normal_rng(Cpred_intestin[t,rep], sigma[1]) ;
      Cgen_caecum[t,rep] = normal_rng(Cpred_caecum[t,rep], sigma[2]) ;
      Cgen_cephalon[t,rep] = normal_rng(Cpred_cephalon[t,rep], sigma[3]) ;
      Cgen_reste[t,rep] = normal_rng(Cpred_reste[t,rep], sigma[4]) ;
    }
  }
}

